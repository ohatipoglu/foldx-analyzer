import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
from plip.exchange.report import BindingSiteReport
from plip.structure.preparation import PDBComplex

from constants import (
    DEFAULT_NEIGHBOR_THRESHOLD,
    DEFAULT_TARGET_RESNUM,
    DEFAULT_TARGET_CHAIN,
    PYMOL_WIDTH,
    PYMOL_HEIGHT,
    PYMOL_EXECUTABLE_PATH
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# PLIP TXT Report Parsing (For Batch Processing)
# ---------------------------------------------------------------------------

def parse_plip_txt_file(filepath, target_residues=None, target_chain='A'):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    results = {'Hydrophobic': [], 'Hydrogen Bond': [], 'Salt Bridge': [], 'Water Bridge': [], 'pi-Cation': [], 'pi-Stacking': []}
    sections = [
        (r'\*\*Hydrophobic Interactions\*\*', 'Hydrophobic', 0, 1, 2, 3, 4, 5, 6),
        (r'\*\*Hydrogen Bonds\*\*', 'Hydrogen Bond', 0, 1, 2, 3, 4, 5, 7),
        (r'\*\*Salt Bridges\*\*', 'Salt Bridge', 0, 1, 2, 4, 5, 6, 7),
        (r'\*\*Water Bridges\*\*', 'Water Bridge', 0, 1, 2, 3, 4, 5, 6),
        (r'\*\*pi-Cation Interactions\*\*', 'pi-Cation', 0, 1, 2, 4, 5, 6, 7),
        (r'\*\*pi-Stacking Interactions\*\*', 'pi-Stacking', 0, 1, 2, 4, 5, 6, 7)
    ]
    for section_regex, key, resnr_idx, restype_idx, chain_idx, ligres_idx, ligtype_idx, ligchain_idx, dist_idx in sections:
        match = re.search(section_regex + r'.*?(?=\*\*|\Z)', content, re.DOTALL)
        if match:
            lines = match.group().split('\n')
            for line in lines:
                if line.startswith('|') and not line.startswith('| RESNR') and not line.startswith('+'):
                    cells = [c.strip() for c in line.split('|') if c.strip()]
                    try:
                        resnr = int(cells[resnr_idx])
                        if cells[chain_idx] == target_chain and (target_residues is None or resnr in target_residues):
                            results[key].append({
                                'prot_res': resnr, 'prot_type': cells[restype_idx], 'prot_chain': cells[chain_idx],
                                'lig_res': int(cells[ligres_idx]), 'lig_type': cells[ligtype_idx], 'lig_chain': cells[ligchain_idx],
                                'distance': cells[dist_idx]
                            })
                    except: continue
    return results

def get_mutation_name(filename):
    name = os.path.splitext(filename)[0]
    return name.replace('target_pa_', '').replace('_model.000.00', '')

# ---------------------------------------------------------------------------
# Word Reports
# ---------------------------------------------------------------------------

def create_plip_comparison_word_report(all_results, output_docx, base_residues=None):
    try:
        from docx import Document
        from docx.enum.text import WD_ALIGN_PARAGRAPH
    except ImportError: return False
    doc = Document()
    title = doc.add_heading('PLIP Interaction Comparison Report', 0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    if base_residues: doc.add_paragraph(f'Target residues (Chain A): {", ".join(map(str, sorted(base_residues)))} (±1)')
    
    itypes = ['Hydrophobic', 'Hydrogen Bond', 'Salt Bridge', 'Water Bridge', 'pi-Cation']
    mutations = sorted(all_results.keys())

    doc.add_heading('Interaction Counts', level=1)
    table = doc.add_table(rows=len(mutations)+1, cols=len(itypes)+1); table.style='Table Grid'
    for i, t in enumerate(['Mutation']+itypes): table.rows[0].cells[i].text = t
    for r, m in enumerate(mutations):
        table.rows[r+1].cells[0].text = m
        for c, t in enumerate(itypes): table.rows[r+1].cells[c+1].text = str(len(all_results[m].get(t, [])))
    
    doc.save(output_docx)
    return True

# ---------------------------------------------------------------------------
# Live PDB Analysis (Biopython & PLIP)
# ---------------------------------------------------------------------------

def load_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    return parser.get_structure("complex", pdb_file)

def list_all_cys(structure):
    res = []
    for m in structure:
        for c in m:
            for r in c:
                if r.get_resname() == "CYS": res.append((m.id, c.id, r.get_id()[1]))
    return res

def find_neighbors(structure, threshold, target_resnum, target_chain):
    target_atoms = []
    other_atoms = []
    for m in structure:
        for c in m:
            for r in c:
                if r.get_resname() == "CYS" and r.get_id()[1] == target_resnum and c.id == target_chain:
                    target_atoms.extend(list(r.get_atoms()))
                else: other_atoms.extend(list(r.get_atoms()))
    if not target_atoms: return []
    ns = NeighborSearch(other_atoms)
    found = set()
    for a in target_atoms:
        for r in ns.search(a.coord, threshold, level='R'):
            if r.get_parent().id != target_chain: found.add((r.get_resname(), r.get_id()[1], r.get_parent().id))
    return sorted(found)

def calculate_disulfide_distance(structure, c1_ch, c1_r, c2_ch, c2_r):
    sgs = {}
    for m in structure:
        for c in m:
            for r in c:
                if r.get_resname() == "CYS":
                    if c.id == c1_ch and r.get_id()[1] == c1_r: sgs['c1'] = next((a for a in r if a.get_name()=="SG"), None)
                    if c.id == c2_ch and r.get_id()[1] == c2_r: sgs['c2'] = next((a for a in r if a.get_name()=="SG"), None)
    if 'c1' in sgs and 'c2' in sgs and sgs['c1'] and sgs['c2']:
        return float(np.linalg.norm(sgs['c1'].coord - sgs['c2'].coord))
    return None

def run_live_plip(pdb_file, out_dir, ligand_sel="organic", interactive=False):
    try:
        os.makedirs(out_dir, exist_ok=True)
        mol = PDBComplex(); mol.load_pdb(pdb_file); mol.analyze()
        out = ["Live PLIP Analysis:\n"]; hotspots = set()
        for sid, site in mol.interaction_sets.items():
            report = BindingSiteReport(site)
            for attr in ['hydrophobic', 'hbond', 'saltbridge', 'pistacking', 'pication']:
                for inter in getattr(report, f"{attr}_interactions", []):
                    hotspots.add((str(inter.reschain), int(inter.resnr)))
                    out.append(f"  {attr.capitalize()}: {inter.restype}{inter.resnr}({inter.reschain}) dist:{getattr(inter, 'distance', 0):.2f}\n")
        
        img_path = os.path.join(out_dir, f"{Path(pdb_file).stem}_plip.png")
        logger.info(f"Generating PLIP image: {img_path}")
        script = write_pymol_script(pdb_file, img_path, ligand_sel, hotspots, interactive=interactive)
        logger.info(f"Using PyMOL script: {script}")
        exe = _resolve_pymol_executable()
        if interactive:
            # In interactive mode, we want to open the PyMOL window and not wait for it.
            # We don't use -c (headless) here.
            subprocess.Popen([exe, "-q", script])
            out.append("\nPyMOL has been launched in 3D interactive mode.\n")
        else:
            # Headless mode to just generate the image
            subprocess.run([exe, "-cq", script], check=True)
            
        return "".join(out), img_path
    except Exception as e: return f"Error: {e}", None

def write_pymol_script(pdb, img, lig_sel, hotspots, width=PYMOL_WIDTH, height=PYMOL_HEIGHT, interactive=False):
    # Use native absolute paths and escape backslashes for PyMOL's internal Python
    pdb_abs = str(Path(pdb).resolve())
    img_abs = str(Path(img).resolve())
    s_path = str(Path(img).resolve().with_suffix(".pml"))
    
    h_sel = " or ".join([f'(chain "{c}" and resi {r})' for c, r in hotspots])
    
    lines = [
        "reinitialize",
        "python",
        f'cmd.load(r"{pdb_abs}")',
        "python end",
        "hide everything",
        "show cartoon",
        "color gray70",
        "util.cnc",
        f"viewport {width}, {height}",
        "bg_color white"
    ]

    # Handle zoom and visualization
    if h_sel:
        lines += [f'select hp, ({h_sel})', 'show sticks, hp', 'color orange, hp', 'zoom hp, 10']
    elif lig_sel:
        # If no hotspots, try to zoom to ligand if it exists
        lines += [f'show sticks, ({lig_sel})', f'zoom ({lig_sel})']
    else:
        lines.append("zoom all")
        
    lines += [
        f"ray {width}, {height}",
        "python",
        f'cmd.png(r"{img_abs}", dpi=300)',
        "python end"
    ]
    
    if not interactive:
        lines.append("quit")
        
    with open(s_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return s_path

def _resolve_pymol_executable():
    if PYMOL_EXECUTABLE_PATH and shutil.which(PYMOL_EXECUTABLE_PATH): return PYMOL_EXECUTABLE_PATH
    for c in ["pymol", "pymol.exe"]:
        if shutil.which(c): return shutil.which(c)
    return "pymol"

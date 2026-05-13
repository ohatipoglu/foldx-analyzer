import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple

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

def parse_plip_txt_file(
    filepath: str,
    target_residues: Optional[Set[int]] = None,
    target_chain: str = 'A'
) -> Dict[str, List[Dict[str, Any]]]:
    """Parse a PLIP TXT report file and extract interaction data.

    Args:
        filepath: Path to the PLIP TXT report file
        target_residues: Optional set of target residue numbers to filter
        target_chain: Target chain identifier (default: 'A')

    Returns:
        Dictionary with interaction types as keys and lists of interaction details as values
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        logger.error("PLIP TXT file not found: %s", filepath)
        return {'Hydrophobic': [], 'Hydrogen Bond': [], 'Salt Bridge': [], 'Water Bridge': [], 'pi-Cation': [], 'pi-Stacking': []}
    except UnicodeDecodeError as e:
        logger.error("Failed to decode PLIP TXT file %s: %s", filepath, e)
        return {'Hydrophobic': [], 'Hydrogen Bond': [], 'Salt Bridge': [], 'Water Bridge': [], 'pi-Cation': [], 'pi-Stacking': []}
    
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
                    except (ValueError, IndexError, KeyError):
                        continue
    return results

def get_mutation_name(filename: str) -> str:
    """Extract mutation name from a filename.
    
    Args:
        filename: Filename with extension
        
    Returns:
        Mutation name with prefixes/suffixes removed
    """
    name = os.path.splitext(filename)[0]
    return name.replace('target_pa_', '').replace('_model.000.00', '')


def create_plip_comparison_word_report(
    all_results: Dict[str, Dict[str, List[Any]]],
    output_docx: str,
    base_residues: Optional[List[int]] = None
) -> bool:
    """Create a Word document comparing PLIP interaction counts across mutations.
    
    Args:
        all_results: Dictionary mapping mutation names to interaction results
        output_docx: Output path for the Word document
        base_residues: Optional list of target residue numbers to include in report
    
    Returns:
        True if report was created successfully, False otherwise
    """
    try:
        from docx import Document
        from docx.enum.text import WD_ALIGN_PARAGRAPH
    except ImportError:
        logger.error("python-docx not installed. Install with: pip install python-docx")
        return False
    
    try:
        doc = Document()
        title = doc.add_heading('PLIP Interaction Comparison Report', 0)
        title.alignment = WD_ALIGN_PARAGRAPH.CENTER
        if base_residues:
            doc.add_paragraph(f'Target residues (Chain A): {", ".join(map(str, sorted(base_residues)))} (±1)')

        itypes = ['Hydrophobic', 'Hydrogen Bond', 'Salt Bridge', 'Water Bridge', 'pi-Cation']
        mutations = sorted(all_results.keys())

        doc.add_heading('Interaction Counts', level=1)
        table = doc.add_table(rows=len(mutations)+1, cols=len(itypes)+1)
        table.style = 'Table Grid'
        for i, t in enumerate(['Mutation'] + itypes):
            table.rows[0].cells[i].text = t
        for r, m in enumerate(mutations):
            table.rows[r+1].cells[0].text = m
            for c, t in enumerate(itypes):
                table.rows[r+1].cells[c+1].text = str(len(all_results[m].get(t, [])))

        doc.save(output_docx)
        return True
    except Exception as e:
        logger.error("Failed to create Word report: %s", e)
        return False

# ---------------------------------------------------------------------------
# Live PDB Analysis (Biopython & PLIP)
# ---------------------------------------------------------------------------

def load_structure(pdb_file: str) -> Any:
    """Load a PDB structure using Biopython.

    Args:
        pdb_file: Path to the PDB file

    Returns:
        Biopython Structure object
    """
    parser = PDBParser(QUIET=True)
    return parser.get_structure("complex", pdb_file)


def list_all_cys(structure: Any) -> List[Tuple[Any, str, int]]:
    """List all cysteine residues in the structure.

    Args:
        structure: Biopython Structure object

    Returns:
        List of tuples (model_id, chain_id, residue_number) for each CYS residue
    """
    res: List[Tuple[Any, str, int]] = []
    for m in structure:
        for c in m:
            for r in c:
                if r.get_resname() == "CYS":
                    res.append((m.id, c.id, r.get_id()[1]))
    return res


def find_neighbors(
    structure: Any,
    threshold: float,
    target_resnum: int,
    target_chain: str
) -> List[Tuple[str, int, str]]:
    """Find neighboring residues near a target cysteine residue.

    Args:
        structure: Biopython Structure object
        threshold: Distance threshold in Angstroms
        target_resnum: Target residue number
        target_chain: Target chain identifier

    Returns:
        Sorted list of tuples (residue_name, residue_number, chain_id) for neighbors
    """
    target_atoms = []
    other_atoms = []
    for m in structure:
        for c in m:
            for r in c:
                if r.get_resname() == "CYS" and r.get_id()[1] == target_resnum and c.id == target_chain:
                    target_atoms.extend(list(r.get_atoms()))
                else:
                    other_atoms.extend(list(r.get_atoms()))
    if not target_atoms:
        logger.warning("Target CYS residue not found: chain=%s, resnum=%s", target_chain, target_resnum)
        return []
    ns = NeighborSearch(other_atoms)
    found = set()
    for a in target_atoms:
        for r in ns.search(a.coord, threshold, level='R'):
            if r.get_parent().id != target_chain:
                found.add((r.get_resname(), r.get_id()[1], r.get_parent().id))
    return sorted(found)


def calculate_disulfide_distance(
    structure: Any,
    c1_ch: str,
    c1_r: int,
    c2_ch: str,
    c2_r: int
) -> Optional[float]:
    """Calculate the distance between two cysteine SG atoms (disulfide bridge).

    Args:
        structure: Biopython Structure object
        c1_ch: First cysteine chain identifier
        c1_r: First cysteine residue number
        c2_ch: Second cysteine chain identifier
        c2_r: Second cysteine residue number

    Returns:
        Distance in Angstroms, or None if SG atoms not found
    """
    sgs = {}
    for m in structure:
        for c in m:
            for r in c:
                if r.get_resname() == "CYS":
                    if c.id == c1_ch and r.get_id()[1] == c1_r:
                        sgs['c1'] = next((a for a in r if a.get_name() == "SG"), None)
                    if c.id == c2_ch and r.get_id()[1] == c2_r:
                        sgs['c2'] = next((a for a in r if a.get_name() == "SG"), None)
    if 'c1' in sgs and 'c2' in sgs and sgs['c1'] and sgs['c2']:
        return float(np.linalg.norm(sgs['c1'].coord - sgs['c2'].coord))
    logger.warning(
        "Could not calculate disulfide distance: C1(%s,%s)-C2(%s,%s)",
        c1_ch, c1_r, c2_ch, c2_r
    )
    return None

def run_live_plip(
    pdb_file: str,
    out_dir: str,
    ligand_sel: str = "organic",
    interactive: bool = False
) -> Tuple[str, Optional[str]]:
    """Run PLIP analysis on a PDB file and generate visualization.

    Args:
        pdb_file: Path to the PDB file
        out_dir: Output directory for results
        ligand_sel: Ligand selection string for PyMOL
        interactive: If True, keep PyMOL open after analysis

    Returns:
        Tuple of (analysis_text, image_path) or (error_message, None)
    """
    try:
        os.makedirs(out_dir, exist_ok=True)
        mol = PDBComplex()
        mol.load_pdb(pdb_file)
        mol.analyze()
        out = ["Live PLIP Analysis:\n"]
        hotspots = set()
        for sid, site in mol.interaction_sets.items():
            report = BindingSiteReport(site)
            for attr in ['hydrophobic', 'hbond', 'saltbridge', 'pistacking', 'pication']:
                for inter in getattr(report, f"{attr}_interactions", []):
                    hotspots.add((str(inter.reschain), int(inter.resnr)))
                    out.append(f"  {attr.capitalize()}: {inter.restype}{inter.resnr}({inter.reschain}) dist:{getattr(inter, 'distance', 0):.2f}\n")

        img_path = os.path.join(out_dir, f"{Path(pdb_file).stem}_plip.png")
        logger.info("Generating PLIP image: %s", img_path)
        script = write_pymol_script(pdb_file, img_path, ligand_sel, hotspots, interactive=interactive)
        logger.info("Using PyMOL script: %s", script)
        exe = _resolve_pymol_executable()
        if interactive:
            subprocess.Popen([exe, "-q", script])
            out.append("\nPyMOL has been launched in 3D interactive mode.\n")
        else:
            subprocess.run([exe, "-cq", script], check=True)

        return "".join(out), img_path
    except FileNotFoundError as e:
        logger.error("File not found during PLIP analysis: %s", e)
        return f"Error: File not found - {e}", None
    except subprocess.CalledProcessError as e:
        logger.error("PyMOL execution failed: %s", e)
        return f"Error: PyMOL execution failed - {e}", None
    except Exception as e:
        logger.exception("Unexpected error during PLIP analysis")
        return f"Error: {e}", None

def write_pymol_script(
    pdb: str,
    img: str,
    lig_sel: str,
    hotspots: Set[Tuple[str, int]],
    width: int = PYMOL_WIDTH,
    height: int = PYMOL_HEIGHT,
    interactive: bool = False
) -> str:
    """Generate a PyMOL script for visualizing PLIP interactions.

    Args:
        pdb: Path to the PDB file
        img: Output image path
        lig_sel: Ligand selection string
        hotspots: Set of (chain, residue) tuples for hotspot residues
        width: Image width in pixels
        height: Image height in pixels
        interactive: If True, don't quit PyMOL after rendering

    Returns:
        Path to the generated .pml script file
    """
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

    if h_sel:
        lines += [f'select hp, ({h_sel})', 'show sticks, hp', 'color orange, hp', 'zoom hp, 10']
    elif lig_sel:
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


def _resolve_pymol_executable() -> str:
    """Resolve the PyMOL executable path from config or system PATH.

    Returns:
        PyMOL executable path or command
    """
    if PYMOL_EXECUTABLE_PATH and shutil.which(PYMOL_EXECUTABLE_PATH):
        return PYMOL_EXECUTABLE_PATH
    for c in ["pymol", "pymol.exe"]:
        if shutil.which(c):
            return shutil.which(c)
    return "pymol"

import logging
import os
import shutil
import subprocess
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
from PIL import Image, ImageTk
from plip.exchange.report import BindingSiteReport
from plip.structure.preparation import PDBComplex

from constants import (
    DEFAULT_NEIGHBOR_THRESHOLD,
    DEFAULT_DISULFIDE_MIN,
    DEFAULT_DISULFIDE_MAX,
    DEFAULT_TARGET_RESNUM,
    DEFAULT_TARGET_CHAIN,
    PLIP_IMAGE_PREVIEW_SIZE,
    PYMOL_WIDTH,
    PYMOL_HEIGHT,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# PyMOL script helpers
# ---------------------------------------------------------------------------

def _pymol_quote(path: str) -> str:
    s = str(path)
    s = s.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{s}"'


def _resolve_pymol_executable() -> str:
    for candidate in ("pymol", "pymol.exe", "PyMOL"):
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    return "pymol"


def _build_residue_selection(residues: set[tuple[str, int]]) -> str:
    if not residues:
        return ""
    parts: list[str] = []
    for chain, resi in sorted(residues, key=lambda x: (str(x[0]), int(x[1]))):
        chain_esc = str(chain).replace("\\", "\\\\").replace('"', '\\"')
        parts.append(f'(chain "{chain_esc}" and resi {int(resi)})')
    return " or ".join(parts)


def _write_pymol_script(
    pdb_file: str,
    image_file: str,
    ligand_selection: str = "organic",
    hotspot_residues: set[tuple[str, int]] | None = None,
    width: int = PYMOL_WIDTH,
    height: int = PYMOL_HEIGHT,
) -> str:
    pdb_path = str(Path(pdb_file).resolve())
    img_path = str(Path(image_file).resolve())
    script_path = str(Path(img_path).with_suffix(".pml"))

    ligand_sel = (ligand_selection or "organic").strip() or "organic"
    hotspot_sel = _build_residue_selection(hotspot_residues or set())

    script_lines = [
        "reinitialize",
        f"load {_pymol_quote(pdb_path)}",
        "hide everything",
        "show cartoon, polymer.protein",
        "color gray70, polymer.protein",
        f"show sticks, ({ligand_sel})",
        f"util.cnc ({ligand_sel})",
        "set stick_radius, 0.18",
        "bg_color white",
        "set antialias, 2",
        "set ray_opaque_background, off",
        f"viewport {width}, {height}",
    ]

    if hotspot_sel:
        script_lines += [
            f"select plip_hotspot, ({hotspot_sel}) and polymer.protein",
            "show sticks, plip_hotspot",
            "color tv_orange, plip_hotspot",
            "set stick_radius, 0.22, plip_hotspot",
        ]

    zoom_target = (f"({ligand_sel}) or plip_hotspot"
                   if hotspot_sel else f"({ligand_sel})")
    script_lines += [
        f"zoom {zoom_target}, 12",
        f"ray {width}, {height}",
        f"png {_pymol_quote(img_path)}, dpi=300",
        "quit",
    ]

    with open(script_path, "w", encoding="utf-8") as f:
        f.write("\n".join(script_lines) + "\n")

    return script_path


# ---------------------------------------------------------------------------
# PDB structure functions
# ---------------------------------------------------------------------------

def load_structure(pdb_file: str):
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", pdb_file)
        return structure
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: PDB file '{pdb_file}' not found.")
    except Exception as e:
        raise Exception(f"Error loading PDB file '{pdb_file}': {e}")


def list_all_cys(structure) -> tuple[str, list[tuple]]:
    cys_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "CYS":
                    res_id = residue.get_id()[1]
                    cys_list.append((model.id, chain.id, res_id))
    if not cys_list:
        return "Warning: No CYS residues found in the structure.\n", cys_list
    return "", cys_list


def find_neighbors(
    structure,
    threshold: float = DEFAULT_NEIGHBOR_THRESHOLD,
    target_resnum: int = DEFAULT_TARGET_RESNUM,
    target_chain: str = DEFAULT_TARGET_CHAIN,
) -> tuple[str, list[tuple]]:
    """Find residues within *threshold* Å of the target CYS residue in other chains."""
    target_atoms = []
    all_other_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if (residue.get_resname() == "CYS"
                        and residue.get_id()[1] == target_resnum
                        and chain.id == target_chain):
                    target_atoms.extend(list(residue.get_atoms()))
                else:
                    all_other_atoms.extend(list(residue.get_atoms()))

    if not target_atoms:
        raise ValueError(
            f"CYS {target_resnum} in chain {target_chain} not found.")

    ns = NeighborSearch(all_other_atoms)
    found_residues: set = set()
    for atom in target_atoms:
        for res in ns.search(atom.coord, threshold, level='R'):
            chain = res.get_parent()
            if chain.id != target_chain:
                found_residues.add((res.get_resname(), res.get_id()[1], chain.id))

    warning = ""
    if not found_residues:
        warning = (
            f"Warning: No neighboring residues found within {threshold} Å "
            f"of CYS {target_resnum} (chain {target_chain}) in other chains.\n"
        )
    return warning, sorted(found_residues)


def calculate_disulfide_distance(
    structure,
    cys1_chain: str = "A",
    cys1_res: int = 166,
    cys2_chain: str = "B",
    cys2_res: int = 56,
) -> float:
    cys1_sg = None
    cys2_sg = None

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "CYS":
                    if chain.id == cys1_chain and residue.get_id()[1] == cys1_res:
                        cys1_sg = next(
                            (a for a in residue.get_atoms() if a.get_name() == "SG"),
                            None,
                        )
                        if cys1_sg is None:
                            raise ValueError(
                                f"SG atom missing in CYS {cys1_res} (chain {cys1_chain}).")
                    elif chain.id == cys2_chain and residue.get_id()[1] == cys2_res:
                        cys2_sg = next(
                            (a for a in residue.get_atoms() if a.get_name() == "SG"),
                            None,
                        )
                        if cys2_sg is None:
                            raise ValueError(
                                f"SG atom missing in CYS {cys2_res} (chain {cys2_chain}).")

    if cys1_sg is None:
        raise ValueError(f"CYS {cys1_res} in chain {cys1_chain} not found.")
    if cys2_sg is None:
        raise ValueError(f"CYS {cys2_res} in chain {cys2_chain} not found.")

    # np.linalg.norm replaces manual math.sqrt(np.sum(diff * diff))
    return float(np.linalg.norm(cys1_sg.coord - cys2_sg.coord))


# ---------------------------------------------------------------------------
# PLIP analysis
# ---------------------------------------------------------------------------

def _collect_hotspot_residues(report: BindingSiteReport) -> set[tuple[str, int]]:
    residues: set[tuple[str, int]] = set()
    interaction_lists = [
        getattr(report, "hydrophobic_interactions", None),
        getattr(report, "hbond_interactions", None),
        getattr(report, "saltbridge_interactions", None),
        getattr(report, "pistacking_interactions", None),
        getattr(report, "pication_interactions", None),
    ]
    for lst in interaction_lists:
        if not lst:
            continue
        for inter in lst:
            chain = getattr(inter, "reschain", None)
            resnr = getattr(inter, "resnr", None)
            if chain is None or resnr is None:
                continue
            try:
                residues.add((str(chain), int(resnr)))
            except Exception:
                continue
    return residues


def run_plip_analysis(
    pdb_file: str,
    output_dir: str | None = None,
    ligand_selection: str = "organic",
) -> tuple[str, str | None]:
    try:
        mol = PDBComplex()
        mol.load_pdb(pdb_file)
        mol.analyze()

        output = ["PLIP Interaction Analysis Results:\n"]
        all_hotspots: set[tuple[str, int]] = set()

        for site_id, binding_site in mol.interaction_sets.items():
            output.append(f"\nBinding Site {site_id}:\n")
            report = BindingSiteReport(binding_site)
            all_hotspots |= _collect_hotspot_residues(report)

            if report.hydrophobic_interactions:
                output.append("Hydrophobic Interactions:\n")
                for inter in report.hydrophobic_interactions:
                    output.append(
                        f"  Residue {inter.resnr} ({inter.restype}, chain {inter.reschain}) "
                        f"with ligand atom {inter.ligatom} at distance {inter.distance:.2f} Å\n"
                    )

            if report.hbond_interactions:
                output.append("Hydrogen Bonds:\n")
                for inter in report.hbond_interactions:
                    output.append(
                        f"  Residue {inter.resnr} ({inter.restype}, chain {inter.reschain}) "
                        f"with ligand atom {inter.ligatom} at distance {inter.distance:.2f} Å\n"
                    )

            if report.saltbridge_interactions:
                output.append("Salt Bridges:\n")
                for inter in report.saltbridge_interactions:
                    output.append(
                        f"  Residue {inter.resnr} ({inter.restype}, chain {inter.reschain}) "
                        f"with ligand at distance {inter.distance:.2f} Å\n"
                    )

            if report.pistacking_interactions:
                output.append("Pi-Stacking Interactions:\n")
                for inter in report.pistacking_interactions:
                    output.append(
                        f"  Residue {inter.resnr} ({inter.restype}, chain {inter.reschain}) "
                        f"with ligand at distance {inter.distance:.2f} Å\n"
                    )

            if report.pication_interactions:
                output.append("Pi-Cation Interactions:\n")
                for inter in report.pication_interactions:
                    output.append(
                        f"  Residue {inter.resnr} ({inter.restype}, chain {inter.reschain}) "
                        f"with ligand at distance {inter.distance:.2f} Å\n"
                    )

            if not any([
                report.hydrophobic_interactions, report.hbond_interactions,
                report.saltbridge_interactions, report.pistacking_interactions,
                report.pication_interactions,
            ]):
                output.append("  No interactions detected for this binding site.\n")

        pdb_path = Path(pdb_file)
        out_dir = (Path(output_dir).expanduser().resolve()
                   if output_dir else pdb_path.parent.resolve())
        out_dir.mkdir(parents=True, exist_ok=True)

        image_file = str(out_dir / f"{pdb_path.stem}_plip.png")
        pymol_script = _write_pymol_script(
            pdb_file=pdb_file,
            image_file=image_file,
            ligand_selection=ligand_selection,
            hotspot_residues=all_hotspots,
        )

        pymol_exe = _resolve_pymol_executable()
        try:
            subprocess.run(
                [pymol_exe, "-cq", pymol_script],
                check=True, capture_output=True, text=True,
            )
        except FileNotFoundError as e:
            raise FileNotFoundError(
                "PyMOL executable not found. "
                "Install PyMOL and ensure 'pymol' is on PATH."
            ) from e
        except subprocess.CalledProcessError as e:
            stderr = (e.stderr or "").strip()
            stdout = (e.stdout or "").strip()
            details = "\n".join(s for s in [stdout, stderr] if s)
            raise RuntimeError(
                f"PyMOL failed while generating the image.\n{details}".strip()
            ) from e

        output.append(f"\nVisualization image generated: {image_file}. Displayed in GUI.\n")
        return "".join(output), image_file

    except Exception as e:
        return (
            f"Error in PLIP analysis: {e}. "
            "Ensure PyMOL is installed and PDB has a ligand.\n"
        ), None


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

class PDBAnalyzerApp:
    def __init__(self, root: tk.Misc) -> None:
        self.root = root
        self.root.title("PDB Analyzer")
        self.root.geometry("1200x900")

        self.pdb_file = tk.StringVar()
        self.output_dir = tk.StringVar(value=os.getcwd())
        self.ligand_selection = tk.StringVar(value="organic")
        self.threshold = tk.DoubleVar(value=DEFAULT_NEIGHBOR_THRESHOLD)
        self.disulfide_min = tk.DoubleVar(value=DEFAULT_DISULFIDE_MIN)
        self.disulfide_max = tk.DoubleVar(value=DEFAULT_DISULFIDE_MAX)

        # Shared PDB file selector — above the notebook, visible on both tabs
        self._build_shared_header()

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.biopython_frame = tk.Frame(self.notebook)
        self.notebook.add(self.biopython_frame, text="Biopython Analysis")

        self.plip_frame = tk.Frame(self.notebook)
        self.notebook.add(self.plip_frame, text="PLIP Analysis")

        self._build_biopython_tab()
        self._build_plip_tab()
        self._build_results_area()

        self.root.focus_force()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_shared_header(self) -> None:
        """Single PDB file selector shared across both tabs."""
        frame = tk.Frame(self.root)
        frame.pack(pady=10, fill=tk.X, padx=10)
        tk.Label(frame, text="PDB File:").pack(side=tk.LEFT, padx=5)
        tk.Entry(frame, textvariable=self.pdb_file,
                 width=50).pack(side=tk.LEFT, padx=5)
        tk.Button(frame, text="Browse",
                  command=self.browse_pdb).pack(side=tk.LEFT, padx=5)

    def _build_biopython_tab(self) -> None:
        frame = tk.Frame(self.biopython_frame)
        frame.pack(pady=10)
        tk.Label(frame, text="Neighbor Threshold (Å):").grid(
            row=0, column=0, padx=5, pady=5)
        tk.Entry(frame, textvariable=self.threshold,
                 width=10).grid(row=0, column=1, padx=5)
        tk.Label(frame, text="Disulfide Min (Å):").grid(
            row=1, column=0, padx=5, pady=5)
        tk.Entry(frame, textvariable=self.disulfide_min,
                 width=10).grid(row=1, column=1, padx=5)
        tk.Label(frame, text="Disulfide Max (Å):").grid(
            row=2, column=0, padx=5, pady=5)
        tk.Entry(frame, textvariable=self.disulfide_max,
                 width=10).grid(row=2, column=1, padx=5)
        tk.Button(frame, text="Run Biopython Analysis",
                  command=self.run_biopython_analysis).grid(
            row=3, columnspan=2, pady=10)

    def _build_plip_tab(self) -> None:
        frame_out = tk.Frame(self.plip_frame)
        frame_out.pack(pady=5, fill=tk.X)
        tk.Label(frame_out, text="Output Folder:").pack(side=tk.LEFT, padx=5)
        tk.Entry(frame_out, textvariable=self.output_dir,
                 width=50).pack(side=tk.LEFT, padx=5)
        tk.Button(frame_out, text="Browse",
                  command=self.browse_output_dir).pack(side=tk.LEFT, padx=5)

        frame_lig = tk.Frame(self.plip_frame)
        frame_lig.pack(pady=5, fill=tk.X)
        tk.Label(frame_lig, text="Ligand selection (PyMOL):").pack(
            side=tk.LEFT, padx=5)
        tk.Entry(frame_lig, textvariable=self.ligand_selection,
                 width=30).pack(side=tk.LEFT, padx=5)
        tk.Label(frame_lig,
                 text='Examples: organic | resn NAG | hetatm and not polymer'
                 ).pack(side=tk.LEFT, padx=5)

        frame_btn = tk.Frame(self.plip_frame)
        frame_btn.pack(pady=10)
        tk.Button(frame_btn, text="Run PLIP Analysis",
                  command=self._on_run_plip).pack(pady=10)

    def _build_results_area(self) -> None:
        self.biopython_results = scrolledtext.ScrolledText(
            self.biopython_frame, wrap=tk.WORD, width=80, height=20)
        self.biopython_results.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

        self.plip_results_frame = tk.Frame(self.plip_frame)
        self.plip_results_frame.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)
        self.plip_results = scrolledtext.ScrolledText(
            self.plip_results_frame, wrap=tk.WORD, width=80, height=10)
        self.plip_results.pack(side=tk.TOP, fill=tk.X)
        self.plip_image_label = tk.Label(self.plip_results_frame)
        self.plip_image_label.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

    # ------------------------------------------------------------------
    # File browser callbacks
    # ------------------------------------------------------------------

    def browse_pdb(self) -> None:
        file_path = filedialog.askopenfilename(
            filetypes=[("PDB Files", "*.pdb")])
        if file_path:
            self.pdb_file.set(file_path)
            try:
                self.output_dir.set(str(Path(file_path).resolve().parent))
            except Exception:
                pass

    def browse_output_dir(self) -> None:
        folder = filedialog.askdirectory()
        if folder:
            self.output_dir.set(folder)

    # ------------------------------------------------------------------
    # Analysis callbacks
    # ------------------------------------------------------------------

    def run_biopython_analysis(self) -> None:
        pdb_file = self.pdb_file.get()
        if not pdb_file:
            messagebox.showerror("Error", "Please select a PDB file.")
            return
        try:
            structure = load_structure(pdb_file)
            self.biopython_results.delete(1.0, tk.END)

            warning, cys_list = list_all_cys(structure)
            self.biopython_results.insert(tk.END, warning)
            self.biopython_results.insert(tk.END, "All CYS residues:\n")
            for model_id, chain_id, res_id in cys_list:
                self.biopython_results.insert(
                    tk.END, f"Model {model_id}: CYS {res_id} chain {chain_id}\n")

            try:
                neigh_warning, neighbors = find_neighbors(
                    structure, self.threshold.get())
                self.biopython_results.insert(tk.END, f"\n{neigh_warning}")
                self.biopython_results.insert(
                    tk.END,
                    f"CYS {DEFAULT_TARGET_RESNUM} (chain {DEFAULT_TARGET_CHAIN}) "
                    f"neighbors within {self.threshold.get()} Å (other chains):\n",
                )
                for resname, resnum, chain_id in neighbors:
                    self.biopython_results.insert(
                        tk.END, f"{resname} {resnum} chain {chain_id}\n")
            except ValueError as e:
                self.biopython_results.insert(
                    tk.END,
                    f"\nError in neighbor search: {e} Skipping neighbor analysis.\n",
                )

            try:
                dist = calculate_disulfide_distance(structure)
                self.biopython_results.insert(
                    tk.END,
                    f"\nCys166 (Ero1α, chain A) ↔ Cys56 (PDI, chain B) "
                    f"S–S distance: {dist:.2f} Å\n",
                )
                if self.disulfide_min.get() <= dist <= self.disulfide_max.get():
                    self.biopython_results.insert(
                        tk.END, "This distance is suitable for a disulfide bond.\n")
                else:
                    self.biopython_results.insert(
                        tk.END,
                        "This distance is not suitable for a disulfide bond.\n",
                    )
            except ValueError as e:
                self.biopython_results.insert(
                    tk.END,
                    f"\nError in distance calculation: {e} Skipping distance calculation.\n",
                )

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _on_run_plip(self) -> None:
        """Button handler for PLIP analysis (renamed to avoid collision with
        the module-level run_plip_analysis function)."""
        pdb_file = self.pdb_file.get()
        if not pdb_file:
            messagebox.showerror("Error", "Please select a PDB file.")
            return
        try:
            output, image_file = run_plip_analysis(
                pdb_file,
                output_dir=self.output_dir.get(),
                ligand_selection=self.ligand_selection.get(),
            )
            self.plip_results.delete(1.0, tk.END)
            self.plip_results.insert(tk.END, output)
            if image_file and os.path.exists(image_file):
                img = Image.open(image_file)
                img = img.resize(PLIP_IMAGE_PREVIEW_SIZE, Image.Resampling.LANCZOS)
                photo = ImageTk.PhotoImage(img)
                self.plip_image_label.config(image=photo)
                self.plip_image_label.image = photo  # prevent GC
            else:
                self.plip_image_label.config(image='')
                self.plip_image_label.image = None
        except Exception as e:
            messagebox.showerror("Error", f"PLIP analysis failed: {e}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.WARNING,
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
    )
    root = tk.Tk()
    PDBAnalyzerApp(root)
    root.mainloop()

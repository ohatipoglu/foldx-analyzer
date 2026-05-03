import logging
import os
import glob
import threading
from tkinter import filedialog, messagebox
from typing import Any, Dict, List, Optional

import customtkinter as ctk
from PIL import Image

# Import components
from foldx_analysis import FoldXAnalyzerGUI
import pdb_analyzer
from constants import (
    DEFAULT_NEIGHBOR_THRESHOLD,
    DEFAULT_DISULFIDE_MIN,
    DEFAULT_DISULFIDE_MAX,
    PLIP_IMAGE_PREVIEW_SIZE
)

logger = logging.getLogger(__name__)


class CombinedApp(ctk.CTk):
    """Main application window for Bioinformatics Analyzer Suite.
    
    Provides a navigation sidebar with three main modules:
    - FoldX Analyzer
    - PDB & PLIP (Live)
    - PLIP Batch (TXT)
    """
    
    def __init__(self) -> None:
        """Initialize the main application window."""
        super().__init__()

        self.title("Bioinformatics Analyzer Suite")
        self.geometry("1450x950")

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # Sidebar
        self.navigation_frame = ctk.CTkFrame(self, corner_radius=0)
        self.navigation_frame.grid(row=0, column=0, sticky="nsew")
        self.navigation_frame.grid_rowconfigure(5, weight=1)

        self.navigation_frame_label = ctk.CTkLabel(self.navigation_frame, text="  Main Menu",
                                                 font=ctk.CTkFont(size=15, weight="bold"))
        self.navigation_frame_label.grid(row=0, column=0, padx=20, pady=20)

        # Nav Buttons
        self.foldx_btn = self._add_nav_btn("FoldX Analyzer", 1, "foldx")
        self.pdb_live_btn = self._add_nav_btn("PDB & PLIP (Live)", 2, "pdb_live")
        self.plip_batch_btn = self._add_nav_btn("PLIP Batch (TXT)", 3, "plip_batch")

        # Appearance
        self.appearance_menu = ctk.CTkOptionMenu(self.navigation_frame, values=["System", "Dark", "Light"], command=ctk.set_appearance_mode)
        self.appearance_menu.grid(row=6, column=0, padx=20, pady=20)

        # Frames
        self.foldx_frame = ctk.CTkFrame(self, corner_radius=0, fg_color="transparent")
        self.foldx_gui = FoldXAnalyzerGUI(self.foldx_frame)

        self.pdb_live_frame = ctk.CTkFrame(self, corner_radius=0, fg_color="transparent")
        self._setup_pdb_live_ui()

        self.plip_batch_frame = ctk.CTkFrame(self, corner_radius=0, fg_color="transparent")
        self._setup_plip_batch_ui()

        self.select_frame_by_name("foldx")

    def _add_nav_btn(self, text: str, row: int, name: str) -> ctk.CTkButton:
        """Create a navigation button for the sidebar.
        
        Args:
            text: Button label text
            row: Grid row position
            name: Frame name to switch to when clicked
            
        Returns:
            Created button widget
        """
        btn = ctk.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=40,
            text=text,
            fg_color="transparent",
            text_color=("gray10", "gray90"),
            anchor="w",
            command=lambda: self.select_frame_by_name(name)
        )
        btn.grid(row=row, column=0, sticky="ew")
        return btn

    def select_frame_by_name(self, name: str) -> None:
        """Switch to the specified frame by name.
        
        Args:
            name: Name of the frame to display ('foldx', 'pdb_live', or 'plip_batch')
        """
        btns = {"foldx": self.foldx_btn, "pdb_live": self.pdb_live_btn, "plip_batch": self.plip_batch_btn}
        frames = {"foldx": self.foldx_frame, "pdb_live": self.pdb_live_frame, "plip_batch": self.plip_batch_frame}
        
        for k, b in btns.items(): b.configure(fg_color=("gray75", "gray25") if k == name else "transparent")
        for k, f in frames.items():
            if k == name:
                f.grid(row=0, column=1, sticky="nsew")
            else:
                f.grid_forget()

    # ---------------------------------------------------------------------------
    # 1. PDB & PLIP (Live) UI
    # ---------------------------------------------------------------------------
    
    def _setup_pdb_live_ui(self) -> None:
        """Set up the PDB & PLIP (Live) analysis user interface."""
        # Default path: data/PDB
        default_pdb_dir = os.path.join(os.getcwd(), "data", "PDB")
        self.live_pdb_var = ctk.StringVar(value=default_pdb_dir)
        self.live_out_dir = ctk.StringVar(value=os.path.join(default_pdb_dir, "output"))
        self.live_thresh = ctk.DoubleVar(value=DEFAULT_NEIGHBOR_THRESHOLD)
        self.live_ligand = ctk.StringVar(value="organic")
        self.live_interactive = ctk.BooleanVar(value=False)

        header = ctk.CTkFrame(self.pdb_live_frame, fg_color="transparent")
        header.pack(pady=20, fill='x', padx=30)
        ctk.CTkLabel(header, text="Live PDB & PLIP Analysis", font=ctk.CTkFont(size=20, weight="bold")).pack(side='left')

        f_row = ctk.CTkFrame(self.pdb_live_frame); f_row.pack(fill='x', padx=30, pady=10)
        ctk.CTkEntry(f_row, textvariable=self.live_pdb_var, width=500).pack(side='left', padx=10)
        ctk.CTkButton(f_row, text="Browse PDB", command=self._browse_live_pdb).pack(side='left')

        tabs = ctk.CTkTabview(self.pdb_live_frame); tabs.pack(fill='both', expand=True, padx=30, pady=10)
        t_bio = tabs.add("Biopython (CYS/Distance)"); t_plip = tabs.add("Live PLIP (3D)")

        # Biopython Tab
        ctk.CTkButton(t_bio, text="Run Structural Analysis", command=self._run_live_bio).pack(pady=10)
        self.live_bio_txt = ctk.CTkTextbox(t_bio, wrap='word'); self.live_bio_txt.pack(fill='both', expand=True, padx=10, pady=10)

        # PLIP Tab
        p_inputs = ctk.CTkFrame(t_plip); p_inputs.pack(fill='x', pady=5)
        ctk.CTkLabel(p_inputs, text="Ligand:").pack(side='left', padx=5)
        ctk.CTkEntry(p_inputs, textvariable=self.live_ligand).pack(side='left', padx=5)
        ctk.CTkCheckBox(p_inputs, text="Interactive (Keep PyMOL Open)", variable=self.live_interactive).pack(side='left', padx=10)
        ctk.CTkButton(t_plip, text="Run Live PLIP", command=self._run_live_plip).pack(pady=5)
        self.live_plip_txt = ctk.CTkTextbox(t_plip, height=150); self.live_plip_txt.pack(fill='x', padx=10)
        self.live_plip_img = ctk.CTkLabel(t_plip, text="Preview")
        self.live_plip_img.pack(fill='both', expand=True)

    def _browse_live_pdb(self) -> None:
        """Open file dialog to select a PDB file."""
        path = filedialog.askopenfilename(filetypes=[("PDB Files", "*.pdb")])
        if path:
            self.live_pdb_var.set(path)
            self.live_out_dir.set(os.path.join(os.path.dirname(path), "output"))

    def _run_live_bio(self) -> None:
        """Run Biopython structural analysis on the selected PDB file."""
        p = self.live_pdb_var.get()
        if not p or not os.path.exists(p):
            return
        s = pdb_analyzer.load_structure(p)
        self.live_bio_txt.delete("1.0", "end")
        cys = pdb_analyzer.list_all_cys(s)
        self.live_bio_txt.insert("end", f"Found {len(cys)} CYS residues.\n")
        dist = pdb_analyzer.calculate_disulfide_distance(s, "A", 166, "B", 56)
        if dist:
            self.live_bio_txt.insert("end", f"\nCys166(A)-Cys56(B) Distance: {dist:.2f} Å\n")

    def _run_live_plip(self) -> None:
        """Run live PLIP analysis on the selected PDB file in a background thread."""
        p = self.live_pdb_var.get()
        if not p or not os.path.exists(p):
            return
        out_dir = self.live_out_dir.get()
        os.makedirs(out_dir, exist_ok=True)
        self.live_plip_txt.delete("1.0", "end")
        self.live_plip_txt.insert("end", "Analyzing...")
        
        def _task() -> None:
            txt, img = pdb_analyzer.run_live_plip(
                p, out_dir, self.live_ligand.get(), interactive=self.live_interactive.get()
            )
            self.after(0, lambda: self._show_live_plip(txt, img))
        
        threading.Thread(target=_task, daemon=True).start()

    def _show_live_plip(self, txt: str, img: Optional[str]) -> None:
        """Display PLIP analysis results and preview image.
        
        Args:
            txt: Analysis text output
            img: Path to generated image (or None if failed)
        """
        self.live_plip_txt.delete("1.0", "end"); self.live_plip_txt.insert("end", txt)
        if img and os.path.exists(img):
            pi = Image.open(img); ci = ctk.CTkImage(pi, pi, size=PLIP_IMAGE_PREVIEW_SIZE)
            self.live_plip_img.configure(image=ci, text="")

    # ---------------------------------------------------------------------------
    # 2. PLIP Batch (TXT) UI
    # ---------------------------------------------------------------------------
    
    def _setup_plip_batch_ui(self) -> None:
        """Set up the PLIP Batch (TXT) report processing user interface."""
        # Default path: data/PLIP_TXT
        default_batch_dir = os.path.join(os.getcwd(), "data", "PLIP_TXT")
        self.batch_dir = ctk.StringVar(value=default_batch_dir)
        self.batch_out_dir = ctk.StringVar(value=os.path.join(default_batch_dir, "output"))
        self.batch_targets = ctk.StringVar(value="50, 53, 63, 90, 94, 95, 96, 98, 100, 117, 119, 125")
        self.batch_data = {}

        header = ctk.CTkFrame(self.plip_batch_frame, fg_color="transparent")
        header.pack(pady=20, fill='x', padx=30)
        ctk.CTkLabel(header, text="PLIP Batch Report Processing", font=ctk.CTkFont(size=20, weight="bold")).pack(side='left')

        f_row = ctk.CTkFrame(self.plip_batch_frame); f_row.pack(fill='x', padx=30, pady=10)
        ctk.CTkEntry(f_row, textvariable=self.batch_dir, width=500).pack(side='left', padx=10)
        ctk.CTkButton(f_row, text="Select Folder", command=self._browse_batch_folder).pack(side='left')

        c_row = ctk.CTkFrame(self.plip_batch_frame); c_row.pack(fill='x', padx=30, pady=5)
        ctk.CTkLabel(c_row, text="Target Residues (A):").pack(side='left', padx=5)
        ctk.CTkEntry(c_row, textvariable=self.batch_targets, width=300).pack(side='left', padx=5)
        ctk.CTkButton(c_row, text="Analyze TXTs", command=self._run_batch).pack(side='left', padx=10)
        self.batch_export_btn = ctk.CTkButton(c_row, text="Export Comparison", state="disabled", command=self._export_batch)
        self.batch_export_btn.pack(side='left', padx=5)

        self.batch_disp = ctk.CTkTextbox(self.plip_batch_frame, font=ctk.CTkFont(family="Courier", size=12))
        self.batch_disp.pack(fill='both', expand=True, padx=30, pady=20)

    def _browse_batch_folder(self) -> None:
        """Open directory dialog to select PLIP TXT folder."""
        path = filedialog.askdirectory()
        if path:
            self.batch_dir.set(path)
            self.batch_out_dir.set(os.path.join(path, "output"))

    def _run_batch(self) -> None:
        """Analyze all PLIP TXT files in the selected folder."""
        d = self.batch_dir.get()
        if not d or not os.path.exists(d):
            return
        files = glob.glob(os.path.join(d, "*.txt"))
        if not files:
            messagebox.showwarning("Warning", "No .txt files found.")
            return

        t_set: set[int] = set()
        b_res: list[int] = []
        for r in self.batch_targets.get().split(','):
            try:
                v = int(r.strip())
                b_res.append(v)
                t_set.add(v)
                t_set.add(v - 1)
                t_set.add(v + 1)
            except ValueError:
                continue

        self.batch_data: Dict[str, Dict[str, list]] = {}
        self.batch_disp.delete("1.0", "end")
        self.batch_disp.insert(
            "end",
            f"{'File':<25} | H-Phobic | H-Bond | Salt | Water | TOTAL\n" + "-" * 80 + "\n"
        )
        for f in sorted(files):
            name = pdb_analyzer.get_mutation_name(os.path.basename(f))
            res = pdb_analyzer.parse_plip_txt_file(f, t_set)
            self.batch_data[name] = res
            counts = [
                len(res.get(k, []))
                for k in ['Hydrophobic', 'Hydrogen Bond', 'Salt Bridge', 'Water Bridge']
            ]
            self.batch_disp.insert(
                "end",
                f"{name:<25} | {counts[0]:^8} | {counts[1]:^6} | {counts[2]:^4} | {counts[3]:^5} | {sum(counts)}\n"
            )

        if self.batch_data:
            self.batch_export_btn.configure(state="normal")

    def _export_batch(self) -> None:
        """Export batch analysis results to a Word document."""
        if not self.batch_data:
            return
        out_dir = self.batch_out_dir.get()
        os.makedirs(out_dir, exist_ok=True)

        default_name = "PLIP_Comparison_Report.docx"
        f = filedialog.asksaveasfilename(
            initialdir=out_dir,
            initialfile=default_name,
            defaultextension=".docx"
        )
        if f:
            success = pdb_analyzer.create_plip_comparison_word_report(self.batch_data, f)
            if success:
                messagebox.showinfo("Success", f"Report saved to:\n{f}")

if __name__ == "__main__":
    app = CombinedApp()
    app.mainloop()

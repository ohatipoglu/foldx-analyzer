import tkinter as tk

from foldx_analysis import FoldXAnalyzerGUI
from pdb_analyzer import PDBAnalyzerApp


def run_combined_gui() -> None:
    """
    Launch a combined GUI with two windows:
    - FoldX Analyzer
    - PDB / PLIP Analyzer

    Her iki uygulama da aynı Tk başlatıcısını (root) paylaşır, ayrı Toplevel pencerelerinde çalışır.
    """
    root = tk.Tk()
    root.withdraw()  # Ana pencereyi gizle, sadece alt pencereleri göster

    # FoldX penceresi
    foldx_win = tk.Toplevel(root)
    foldx_win.title("FoldX Analyzer")
    FoldXAnalyzerGUI(foldx_win)

    # PDB / PLIP penceresi
    pdb_win = tk.Toplevel(root)
    pdb_win.title("PDB / PLIP Analyzer")
    PDBAnalyzerApp(pdb_win)

    # Root'u görünmez ama event-loop sahibi olarak kullan
    root.mainloop()


if __name__ == "__main__":
    run_combined_gui()


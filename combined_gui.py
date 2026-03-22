import logging
import tkinter as tk

from foldx_analysis import FoldXAnalyzerGUI
from pdb_analyzer import PDBAnalyzerApp


def run_combined_gui() -> None:
    """Launch FoldX Analyzer and PDB/PLIP Analyzer in two coordinated windows.

    Both share a single hidden Tk root as the event-loop owner.  Closing
    either child window destroys the entire application cleanly.
    """
    root = tk.Tk()
    root.withdraw()

    foldx_win = tk.Toplevel(root)
    foldx_win.title("FoldX Analyzer")
    FoldXAnalyzerGUI(foldx_win)

    pdb_win = tk.Toplevel(root)
    pdb_win.title("PDB / PLIP Analyzer")
    PDBAnalyzerApp(pdb_win)

    def _on_close() -> None:
        root.destroy()

    foldx_win.protocol("WM_DELETE_WINDOW", _on_close)
    pdb_win.protocol("WM_DELETE_WINDOW", _on_close)
    root.protocol("WM_DELETE_WINDOW", _on_close)

    root.mainloop()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.WARNING,
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
    )
    run_combined_gui()

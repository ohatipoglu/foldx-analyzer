import logging
import tkinter as tk
from tkinter import messagebox

import customtkinter as ctk

from foldx_analysis import FoldXAnalyzerGUI
from pdb_analyzer import PDBAnalyzerApp


def run_combined_gui() -> None:
    """Launch FoldX Analyzer and PDB/PLIP Analyzer in two coordinated windows."""
    try:
        ctk.set_appearance_mode("System")
        ctk.set_default_color_theme("blue")

        root = ctk.CTk()
        root.withdraw()

        foldx_win = ctk.CTkToplevel(root)
        foldx_win.title("FoldX Analyzer")
        FoldXAnalyzerGUI(foldx_win)

        pdb_win = ctk.CTkToplevel(root)
        pdb_win.title("PDB / PLIP Analyzer")
        PDBAnalyzerApp(pdb_win)

        def _on_close() -> None:
            root.destroy()

        foldx_win.protocol("WM_DELETE_WINDOW", _on_close)
        pdb_win.protocol("WM_DELETE_WINDOW", _on_close)
        root.protocol("WM_DELETE_WINDOW", _on_close)

        root.mainloop()

    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}", exc_info=True)
        # Use a fallback Tkinter window for the error message if CTk fails
        root_tk = tk.Tk()
        root_tk.withdraw()
        messagebox.showerror(
            "Fatal Error",
            "An unexpected error occurred. Please check the logs for more details."
        )
        root_tk.destroy()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        handlers=[
            logging.FileHandler("app.log"),
            logging.StreamHandler()
        ]
    )
    run_combined_gui()

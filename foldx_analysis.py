"""
FoldX Analyzer GUI — UI layer only.

All data-processing logic lives in foldx_core.py.  This module is
responsible for:
  - Building and managing the Tkinter interface
  - Running file processing in a background thread
  - Rendering results as matplotlib charts embedded in tabs
  - Saving PNG and CSV outputs
"""

import logging
import os
import threading

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from constants import (
    ALL_COMMANDS, AMINO_ACIDS, RNA_BASES,
    BG_COLOR, STATUS_BG, STATUS_FG,
    FIGURE_SIZE, FIGURE_DPI, SAVE_DPI,
    BUILDMODEL_MAX_N, ANALYSECOMPLEX_MAX_N,
    ENERGY_THRESHOLD_GREEN, ENERGY_THRESHOLD_RED,
    CMD_UNKNOWN,
)
from foldx_core import (
    collect_files, find_header_and_read,
    detect_command, PROCESSORS,
)

logger = logging.getLogger(__name__)


class FoldXAnalyzerGUI:
    def __init__(self, root: tk.Misc) -> None:
        self.root = root
        self.root.title("FoldX Analyzer v8")
        self.root.geometry("1200x900")
        self.root.configure(bg=BG_COLOR)

        ttk.Style().theme_use('clam')

        self._build_ui()

        # tab_name → FigureCanvasTkAgg; _tab_frames tracks the Frame widgets
        # so they can be properly destroyed on re-run.
        self.tabs: dict[str, FigureCanvasTkAgg] = {}
        self._tab_frames: list[ttk.Frame] = []

        self._renderers = {
            "PositionScan":   self._render_positionscan,
            "RepairPDB":      self._render_repairpdb,
            "BuildModel":     self._render_buildmodel,
            "AnalyseComplex": self._render_analysecomplex,
            "Stability":      self._render_stability,
            "Pssm":           self._render_pssm,
            "RnaScan":        self._render_rnascan,
        }

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self) -> None:
        top = tk.Frame(self.root, bg=BG_COLOR)
        top.pack(pady=10, fill='x', padx=20)

        # File / folder selector
        tk.Label(top, text="fxout Dosyası/Klasörü:",
                 font=("Helvetica", 10), bg=BG_COLOR).pack(anchor='w')
        file_row = tk.Frame(top, bg=BG_COLOR)
        file_row.pack(fill='x', pady=2)
        self.file_path = tk.StringVar()
        tk.Entry(file_row, textvariable=self.file_path,
                 width=70, state='readonly').pack(side='left', padx=(0, 5))
        ttk.Button(file_row, text="Tek Dosya",
                   command=self.browse_file).pack(side='left', padx=2)
        ttk.Button(file_row, text="Toplu Klasör",
                   command=self.browse_folder).pack(side='left')

        # Mode
        mode_row = tk.Frame(top, bg=BG_COLOR)
        mode_row.pack(fill='x', pady=8)
        tk.Label(mode_row, text="Mod:",
                 font=("Helvetica", 10), bg=BG_COLOR).pack(side='left')
        self.mode_var = tk.StringVar(value="auto")
        ttk.Radiobutton(mode_row, text="Otomatik", variable=self.mode_var,
                        value="auto",
                        command=self._toggle_manual).pack(side='left', padx=(10, 5))
        ttk.Radiobutton(mode_row, text="Manuel", variable=self.mode_var,
                        value="manual",
                        command=self._toggle_manual).pack(side='left', padx=5)
        self.manual_combo = ttk.Combobox(mode_row, values=ALL_COMMANDS,
                                         state="disabled", width=15)
        self.manual_combo.pack(side='left', padx=5)
        self.manual_combo.set(ALL_COMMANDS[0])

        # Graph title
        tk.Label(top, text="Grafik Başlığı:",
                 font=("Helvetica", 10), bg=BG_COLOR).pack(anchor='w', pady=(10, 2))
        self.title_var = tk.StringVar(value="FoldX Analizi")
        tk.Entry(top, textvariable=self.title_var, width=70).pack(pady=2, fill='x')

        # Output directory
        out_row = tk.Frame(top, bg=BG_COLOR)
        out_row.pack(fill='x', pady=8)
        tk.Label(out_row, text="Çıktı Klasörü:",
                 font=("Helvetica", 10), bg=BG_COLOR).pack(side='left')
        self.out_dir = tk.StringVar(value=os.path.join(os.getcwd(), "output"))
        tk.Entry(out_row, textvariable=self.out_dir,
                 width=60).pack(side='left', padx=5)
        ttk.Button(out_row, text="Değiştir",
                   command=self.browse_output).pack(side='left')

        # Run button
        btn_row = tk.Frame(top, bg=BG_COLOR)
        btn_row.pack(pady=15)
        self.analyze_btn = ttk.Button(btn_row, text="Analizi Başlat",
                                      command=self.run_analysis)
        self.analyze_btn.pack()

        # Progress bar
        self.progress = ttk.Progressbar(top, mode='determinate', maximum=100)
        self.progress.pack(fill='x', pady=5)

        # Status bar
        self.status = tk.Label(self.root, text="Hazır.", relief='sunken',
                               anchor='w', bg=STATUS_BG, fg=STATUS_FG)
        self.status.pack(side='bottom', fill='x')

        # Graph tabs
        self.tab_control = ttk.Notebook(self.root)
        self.tab_control.pack(fill='both', expand=True, padx=20, pady=10)

    # ------------------------------------------------------------------
    # UI helpers
    # ------------------------------------------------------------------

    def _toggle_manual(self) -> None:
        state = "readonly" if self.mode_var.get() == "manual" else "disabled"
        self.manual_combo.config(state=state)

    def browse_file(self) -> None:
        path = filedialog.askopenfilename(
            filetypes=[("FoldX Output", "*.fxout"), ("All", "*.*")])
        if path:
            self.file_path.set(path)

    def browse_folder(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.file_path.set(path)

    def browse_output(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.out_dir.set(path)

    def _set_status(self, text: str) -> None:
        self.status.config(text=text)

    def _update_progress(self, value: int) -> None:
        self.progress['value'] = value

    def _clear_tabs(self) -> None:
        for tab_id in self.tab_control.tabs():
            self.tab_control.forget(tab_id)
        for frame in self._tab_frames:
            frame.destroy()
        self._tab_frames.clear()
        self.tabs.clear()

    # ------------------------------------------------------------------
    # Graph renderers — one method per FoldX command
    # ------------------------------------------------------------------

    def _render_positionscan(
        self, ax: plt.Axes, results: list[dict], aux: list[str] | None,
    ) -> None:
        amino_acids = aux or AMINO_ACIDS
        stats = (pd.DataFrame(results)
                 .groupby('amino_acid')['energy']
                 .agg(['mean', 'std', 'count'])
                 .reset_index())
        stats['sem'] = stats['std'] / np.sqrt(stats['count'])
        order = {aa: i for i, aa in enumerate(amino_acids)}
        stats = stats.sort_values(by='amino_acid', key=lambda x: x.map(order))
        bars = ax.bar(stats['amino_acid'], stats['mean'],
                      yerr=stats['sem'], capsize=5,
                      color='lightgray', edgecolor='black')
        for bar, val in zip(bars, stats['mean']):
            bar.set_color(
                'green' if val < ENERGY_THRESHOLD_GREEN
                else 'red' if val > ENERGY_THRESHOLD_RED
                else 'gray'
            )
        ax.axhline(0, color='black', lw=0.8)
        ax.set_ylabel('ΔΔG (kcal/mol)')

    def _render_repairpdb(
        self, ax: plt.Axes, results: list[dict], aux,
    ) -> None:
        df_res = pd.DataFrame(results)
        ax.plot(df_res['model'], df_res['energy'], 'o-', color='blue')
        ax.set_xlabel('Model')
        ax.set_ylabel('Total Energy')
        ax.grid(True, alpha=0.3)

    def _render_buildmodel(
        self, ax: plt.Axes, results: list[dict], aux,
    ) -> None:
        df_res = pd.DataFrame(results)
        if len(df_res) > BUILDMODEL_MAX_N:
            df_res = df_res.iloc[:BUILDMODEL_MAX_N]
            ax.set_title(ax.get_title() + f" (ilk {BUILDMODEL_MAX_N} satır)")
        ax.plot(range(1, len(df_res) + 1), df_res['energy'], 'o-', color='purple')
        ax.set_xlabel('Sıra')
        ax.set_ylabel('Energy')
        ax.grid(True, alpha=0.3)

    def _render_analysecomplex(
        self, ax: plt.Axes, results: list[dict], aux,
    ) -> None:
        df_res = pd.DataFrame(results)
        if len(df_res) > ANALYSECOMPLEX_MAX_N:
            df_res = df_res.iloc[:ANALYSECOMPLEX_MAX_N]
            ax.set_title(ax.get_title() + f" (ilk {ANALYSECOMPLEX_MAX_N} satır)")
        x = np.arange(len(df_res))
        bottom = np.zeros(len(df_res))
        for comp, color in [('backbone', '#1f77b4'),
                             ('sidechain', '#ff7f0e'),
                             ('interaction', '#2ca02c')]:
            if comp in df_res.columns:
                vals = pd.to_numeric(df_res[comp], errors='coerce').fillna(0).to_numpy()
                ax.bar(x, vals, bottom=bottom, label=comp, color=color, alpha=0.85)
                bottom += vals
        labels = (df_res['pdb'].astype(str).tolist()
                  if 'pdb' in df_res.columns else [str(i) for i in x])
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=90, fontsize=8)
        ax.set_ylabel('Energy components (stacked)')
        ax.legend()
        ax.grid(True, axis='y', alpha=0.2)

    def _render_stability(
        self, ax: plt.Axes, results: list[dict], aux,
    ) -> None:
        df_res = pd.DataFrame(results)
        ax.bar(['Stability'], df_res['energy'].astype(float), color='teal')
        ax.set_ylabel('Total Energy')
        ax.grid(True, axis='y', alpha=0.3)

    def _render_pssm(
        self, ax: plt.Axes, results: list[dict], aux,
    ) -> None:
        df_all = pd.DataFrame(results)
        num = df_all.apply(pd.to_numeric, errors='coerce')
        keep_cols = [c for c in num.columns
                     if num[c].notna().sum() >= max(1, int(0.3 * len(num)))]
        num = num[keep_cols].dropna(how='all')
        if num.empty or num.shape[1] < 2:
            ax.text(0.01, 0.99,
                    "PSSM grafik üretilemedi (sayısal veri bulunamadı).",
                    transform=ax.transAxes, va='top', ha='left')
            ax.axis('off')
        else:
            im = ax.imshow(num.to_numpy(), aspect='auto',
                           interpolation='nearest', cmap='viridis')
            ax.set_xlabel('Kolonlar')
            ax.set_ylabel('Satırlar')
            ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    def _render_rnascan(
        self, ax: plt.Axes, results: list[dict], aux: list[str] | None,
    ) -> None:
        bases = aux or RNA_BASES
        df_res = (pd.DataFrame(results)
                  .groupby('base')['energy']
                  .agg(['mean', 'std', 'count'])
                  .reset_index())
        df_res['sem'] = df_res['std'] / np.sqrt(df_res['count'])
        order = {b: i for i, b in enumerate(bases)}
        df_res = df_res.sort_values(by='base', key=lambda x: x.map(order))
        ax.bar(df_res['base'], df_res['mean'],
               yerr=df_res['sem'], capsize=5,
               color='lightgray', edgecolor='black')
        ax.axhline(0, color='black', lw=0.8)
        ax.set_ylabel('ΔΔG (kcal/mol)')

    # ------------------------------------------------------------------
    # Analysis orchestration
    # ------------------------------------------------------------------

    def run_analysis(self) -> None:
        input_path = self.file_path.get()
        if not input_path:
            messagebox.showerror("Hata", "Lütfen dosya/klasör seçin!")
            return

        output_dir = self.out_dir.get()
        os.makedirs(output_dir, exist_ok=True)

        files = collect_files(input_path)
        if not files:
            if not os.path.exists(input_path):
                messagebox.showerror("Hata", "Geçersiz yol!")
            else:
                messagebox.showwarning(
                    "Uyarı", "Seçilen klasörde .fxout dosyası bulunamadı.")
            return

        # Capture all UI state before the worker thread starts.
        mode = self.mode_var.get()
        manual_cmd = self.manual_combo.get()
        graph_title = self.title_var.get() or "FoldX Analizi"

        self._clear_tabs()
        self.analyze_btn.config(state='disabled')
        self.progress.config(maximum=len(files), value=0)

        threading.Thread(
            target=self._analysis_worker,
            args=(files, output_dir, mode, manual_cmd, graph_title),
            daemon=True,
        ).start()

    def _analysis_worker(
        self,
        files: list[str],
        output_dir: str,
        mode: str,
        manual_cmd: str,
        graph_title: str,
    ) -> None:
        """Run in a background thread; posts all UI updates via root.after."""
        total = len(files)
        error_files: list[str] = []

        for idx, file in enumerate(files):
            basename = os.path.basename(file)
            self.root.after(
                0, self._set_status,
                f"İşleniyor: {basename} ({idx + 1}/{total})",
            )
            try:
                df, _, _, header_lines = find_header_and_read(file)
                detected = detect_command(df, header_lines, basename)

                if detected == CMD_UNKNOWN and mode == "auto":
                    self.root.after(
                        0, self._set_status,
                        f"Uyarı: {basename} için komut tespit edilemedi, atlandı.",
                    )
                else:
                    command = manual_cmd if mode == "manual" else detected
                    if command not in PROCESSORS:
                        logger.error("Unknown command '%s' for file %s", command, basename)
                    else:
                        results, aux = PROCESSORS[command](df)
                        if results:
                            self.root.after(
                                0, self._add_result,
                                file, idx, command, results, aux,
                                graph_title, output_dir,
                            )
                        else:
                            self.root.after(
                                0, self._set_status,
                                f"{basename} için geçerli sonuç elde edilemedi.",
                            )
            except Exception:
                logger.exception("Error processing %s", file)
                error_files.append(basename)
                self.root.after(0, self._set_status, f"Hata: {basename}")
            finally:
                self.root.after(0, self._update_progress, idx + 1)

        self.root.after(0, self._analysis_complete, total, output_dir, error_files)

    def _add_result(
        self,
        file: str,
        idx: int,
        command: str,
        results: list[dict],
        aux,
        graph_title: str,
        output_dir: str,
    ) -> None:
        """Called on the main thread: create a tab, render the chart, save files."""
        basename = os.path.basename(file)
        tab_name = basename[:20]
        if tab_name in self.tabs:
            tab_name += f"_{idx}"

        tab = ttk.Frame(self.tab_control)
        self._tab_frames.append(tab)
        self.tab_control.add(tab, text=tab_name)

        fig = plt.Figure(figsize=FIGURE_SIZE, dpi=FIGURE_DPI)
        ax = fig.add_subplot(111)
        ax.set_title(f"{graph_title} - {basename}")

        renderer = self._renderers.get(command)
        if renderer:
            renderer(ax, results, aux)

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)
        self.tabs[tab_name] = canvas

        # Save outputs
        stem = os.path.splitext(basename)[0]
        file_out_dir = os.path.join(output_dir, stem)
        os.makedirs(file_out_dir, exist_ok=True)
        fig.savefig(
            os.path.join(file_out_dir, f"{command}_graph.png"),
            dpi=SAVE_DPI, format='png', bbox_inches='tight',
        )
        pd.DataFrame(results).to_csv(
            os.path.join(file_out_dir, f"{command}_data.csv"), index=False,
        )

    def _analysis_complete(
        self,
        total: int,
        output_dir: str,
        error_files: list[str],
    ) -> None:
        self.analyze_btn.config(state='normal')
        if error_files:
            failed = "\n".join(error_files)
            self._set_status(
                f"Tamamlandı ({total} dosya). {len(error_files)} hata var.")
            messagebox.showerror(
                "Analiz Hataları",
                f"Şu dosyalar işlenemedi:\n{failed}",
            )
        else:
            self._set_status(f"Toplam {total} dosya analiz edildi!")
            messagebox.showinfo(
                "Başarılı!",
                f"{total} dosya analiz edildi!\nÇıktılar: {output_dir}",
            )


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.WARNING,
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
    )
    root = tk.Tk()
    FoldXAnalyzerGUI(root)
    root.mainloop()

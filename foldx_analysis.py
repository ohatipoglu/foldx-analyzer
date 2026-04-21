"""
FoldX Analyzer GUI — UI layer only.

All data-processing logic lives in foldx_core.py.  This module is
responsible for:
  - Building and managing the Tkinter/CustomTkinter interface
  - Running file processing in a background thread
  - Rendering results as matplotlib charts embedded in tabs
  - Saving PNG and CSV outputs
"""

import logging
import os
import threading
from tkinter import filedialog, messagebox

import customtkinter as ctk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from constants import (
    ALL_COMMANDS, AMINO_ACIDS, RNA_BASES,
    CMD_UNKNOWN,
    FIGURE_SIZE, FIGURE_DPI, SAVE_DPI,
    BUILDMODEL_MAX_N, ANALYSECOMPLEX_MAX_N,
    ENERGY_THRESHOLD_GREEN, ENERGY_THRESHOLD_RED,
)
from foldx_core import (
    collect_files, find_header_and_read,
    detect_command, PROCESSORS,
)

logger = logging.getLogger(__name__)


class FoldXAnalyzerGUI:
    def __init__(self, root: ctk.CTkFrame | ctk.CTkToplevel) -> None:
        self.root = root
        # Base settings - only apply if root is a window
        if isinstance(self.root, (ctk.CTk, ctk.CTkToplevel)):
            self.root.geometry("1200x900")
        
        self._build_ui()

        # tab_name → FigureCanvasTkAgg; _tab_frames tracks the Frame widgets
        # so they can be properly destroyed on re-run.
        self.tabs: dict[str, FigureCanvasTkAgg] = {}
        self._tab_frames: list[ctk.CTkFrame] = []

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
        top_frame = ctk.CTkFrame(self.root)
        top_frame.pack(pady=10, fill='x', padx=20)

        # File / folder selector
        ctk.CTkLabel(top_frame, text="fxout Dosyası/Klasörü:", font=ctk.CTkFont(weight="bold")).pack(anchor='w', padx=10, pady=(10, 0))
        file_row = ctk.CTkFrame(top_frame, fg_color="transparent")
        file_row.pack(fill='x', pady=5, padx=10)
        self.file_path = ctk.StringVar(value=os.path.join(os.getcwd(), "data", "fxout"))
        file_entry = ctk.CTkEntry(file_row, textvariable=self.file_path, width=500, state='readonly')
        file_entry.pack(side='left', padx=(0, 10))
        ctk.CTkButton(file_row, text="Tek Dosya", command=self.browse_file, width=100).pack(side='left', padx=5)
        ctk.CTkButton(file_row, text="Toplu Klasör", command=self.browse_folder, width=100).pack(side='left')

        # Mode
        mode_row = ctk.CTkFrame(top_frame, fg_color="transparent")
        mode_row.pack(fill='x', pady=5, padx=10)
        ctk.CTkLabel(mode_row, text="Mod:", font=ctk.CTkFont(weight="bold")).pack(side='left', padx=(0, 10))
        self.mode_var = ctk.StringVar(value="auto")
        ctk.CTkRadioButton(mode_row, text="Otomatik", variable=self.mode_var, value="auto", command=self._toggle_manual).pack(side='left', padx=10)
        ctk.CTkRadioButton(mode_row, text="Manuel", variable=self.mode_var, value="manual", command=self._toggle_manual).pack(side='left', padx=10)
        self.manual_combo = ctk.CTkComboBox(mode_row, values=ALL_COMMANDS, state="disabled", width=150)
        self.manual_combo.pack(side='left', padx=10)
        self.manual_combo.set(ALL_COMMANDS[0])

        # Graph title
        ctk.CTkLabel(top_frame, text="Grafik Başlığı:", font=ctk.CTkFont(weight="bold")).pack(anchor='w', padx=10, pady=(5, 0))
        self.title_var = ctk.StringVar(value="FoldX Analizi")
        ctk.CTkEntry(top_frame, textvariable=self.title_var, width=500).pack(pady=5, padx=10, anchor='w')

        # Output directory
        out_row = ctk.CTkFrame(top_frame, fg_color="transparent")
        out_row.pack(fill='x', pady=5, padx=10)
        ctk.CTkLabel(out_row, text="Çıktı Klasörü:", font=ctk.CTkFont(weight="bold")).pack(side='left', padx=(0, 10))
        self.out_dir = ctk.StringVar(value=os.path.join(os.getcwd(), "data", "fxout", "output"))
        ctk.CTkEntry(out_row, textvariable=self.out_dir, width=400).pack(side='left', padx=(0, 10))
        ctk.CTkButton(out_row, text="Değiştir", command=self.browse_output, width=100).pack(side='left')

        # Run button
        btn_row = ctk.CTkFrame(top_frame, fg_color="transparent")
        btn_row.pack(pady=15)
        self.analyze_btn = ctk.CTkButton(btn_row, text="Analizi Başlat", command=self.run_analysis, font=ctk.CTkFont(weight="bold", size=14))
        self.analyze_btn.pack()

        # Progress bar
        self.progress = ctk.CTkProgressBar(top_frame, mode='determinate')
        self.progress.pack(fill='x', pady=10, padx=20)
        self.progress.set(0)

        # Status bar
        self.status = ctk.CTkLabel(self.root, text="Hazır.", anchor='w', fg_color=("gray85", "gray25"), padx=10)
        self.status.pack(side='bottom', fill='x')

        # Graph tabs (Using CTkTabview)
        self.tab_control = ctk.CTkTabview(self.root)
        self.tab_control.pack(fill='both', expand=True, padx=20, pady=10)

    # ------------------------------------------------------------------
    # UI helpers
    # ------------------------------------------------------------------

    def _update_output_path(self, input_path: str) -> None:
        if os.path.isfile(input_path):
            new_out = os.path.join(os.path.dirname(input_path), "output")
        else:
            new_out = os.path.join(input_path, "output")
        self.out_dir.set(new_out)

    def _toggle_manual(self) -> None:
        state = "normal" if self.mode_var.get() == "manual" else "disabled"
        self.manual_combo.configure(state=state)

    def browse_file(self) -> None:
        path = filedialog.askopenfilename(
            filetypes=[("FoldX Output", "*.fxout"), ("All", "*.*")])
        if path:
            self.file_path.set(path)
            self._update_output_path(path)

    def browse_folder(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.file_path.set(path)
            self._update_output_path(path)

    def browse_output(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.out_dir.set(path)

    def _set_status(self, text: str) -> None:
        self.status.configure(text=text)

    def _update_progress(self, value: float) -> None:
        self.progress.set(value)

    def _clear_tabs(self) -> None:
        for tab_name in list(self.tabs.keys()):
            try:
                self.tab_control.delete(tab_name)
            except ValueError:
                pass
        self.tabs.clear()
        self._tab_frames.clear()

    # ------------------------------------------------------------------
    # Graph renderers — one method per FoldX command
    # ------------------------------------------------------------------

    def _render_positionscan(
        self, ax: plt.Axes, stats_df: pd.DataFrame, aux: list[str] | None,
    ) -> None:
        amino_acids = aux or AMINO_ACIDS
        if stats_df.empty:
            ax.text(0.5, 0.5, "Analiz için veri bulunamadı.",
                    ha='center', va='center', transform=ax.transAxes)
            return

        order = {aa: i for i, aa in enumerate(amino_acids)}
        stats_df = stats_df.sort_values(by='amino_acid', key=lambda x: x.map(order))
        bars = ax.bar(stats_df['amino_acid'], stats_df['mean'],
                      yerr=stats_df['sem'], capsize=5,
                      color='lightgray', edgecolor='black')
        for bar, val in zip(bars, stats_df['mean']):
            bar.set_color(
                'green' if val < ENERGY_THRESHOLD_GREEN
                else 'red' if val > ENERGY_THRESHOLD_RED
                else 'gray'
            )
        ax.axhline(0, color='black', lw=0.8)
        ax.set_ylabel('ΔΔG (kcal/mol)')

    def _render_repairpdb(
        self, ax: plt.Axes, df_res: pd.DataFrame, aux,
    ) -> None:
        ax.plot(df_res['model'], df_res['energy'], 'o-', color='blue')
        ax.set_xlabel('Model')
        ax.set_ylabel('Total Energy')
        ax.grid(True, alpha=0.3)

    def _render_buildmodel(
        self, ax: plt.Axes, df_res: pd.DataFrame, aux,
    ) -> None:
        if len(df_res) > BUILDMODEL_MAX_N:
            df_res = df_res.iloc[:BUILDMODEL_MAX_N]
            ax.set_title(ax.get_title() + f" (ilk {BUILDMODEL_MAX_N} satır)")
        ax.plot(range(1, len(df_res) + 1), df_res['energy'], 'o-', color='purple')
        ax.set_xlabel('Sıra')
        ax.set_ylabel('Energy')
        ax.grid(True, alpha=0.3)

    def _render_analysecomplex(
        self, ax: plt.Axes, df_res: pd.DataFrame, aux,
    ) -> None:
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
        self, ax: plt.Axes, df_res: pd.DataFrame, aux,
    ) -> None:
        ax.bar(['Stability'], df_res['energy'].astype(float), color='teal')
        ax.set_ylabel('Total Energy')
        ax.grid(True, axis='y', alpha=0.3)

    def _render_pssm(
        self, ax: plt.Axes, df_res: pd.DataFrame, aux,
    ) -> None:
        num = df_res.apply(pd.to_numeric, errors='coerce')
        keep_cols = [c for c in num.columns
                     if c in num and num[c].notna().sum() >= max(1, int(0.3 * len(num)))]
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
        self, ax: plt.Axes, stats_df: pd.DataFrame, aux: list[str] | None,
    ) -> None:
        bases = aux or RNA_BASES
        if stats_df.empty:
            ax.text(0.5, 0.5, "Analiz için veri bulunamadı.",
                    ha='center', va='center', transform=ax.transAxes)
            return

        order = {b: i for i, b in enumerate(bases)}
        stats_df = stats_df.sort_values(by='base', key=lambda x: x.map(order))
        ax.bar(stats_df['base'], stats_df['mean'],
               yerr=stats_df['sem'], capsize=5,
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
        try:
            os.makedirs(output_dir, exist_ok=True)
        except Exception as e:
            messagebox.showerror("Hata", f"Çıktı klasörü oluşturulamadı:\n{e}")
            return

        files = collect_files(input_path)
        if not files:
            if not os.path.exists(input_path):
                messagebox.showerror("Hata", "Geçersiz yol!")
            else:
                messagebox.showwarning(
                    "Uyarı", "Seçilen klasörde .fxout dosyası bulunamadı.")
            return

        mode = self.mode_var.get()
        manual_cmd = self.manual_combo.get()
        graph_title = self.title_var.get() or "FoldX Analizi"

        self._clear_tabs()
        self.analyze_btn.configure(state='disabled')
        self.progress.set(0)

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
                        if not results.empty:
                            self.root.after(
                                0, self._add_result,
                                file, idx, command, results, aux,
                                graph_title, output_dir,
                            )
                        else:
                            self.root.after(
                                0, self._set_status,
                                f"{basename} için geçerli sonuç elde edilemedi veya ayrıştırılamadı.",
                            )

            except Exception as e:
                logger.exception("Error processing %s", file)
                error_files.append(basename)
                self.root.after(0, self._set_status, f"Hata: {basename}")
            finally:
                progress_val = (idx + 1) / total
                self.root.after(0, self._update_progress, progress_val)

        self.root.after(0, self._analysis_complete, total, output_dir, error_files)

    def _add_result(
        self,
        file: str,
        idx: int,
        command: str,
        results_df: pd.DataFrame,
        aux,
        graph_title: str,
        output_dir: str,
    ) -> None:
        basename = os.path.basename(file)
        tab_name = basename[:20]
        if tab_name in self.tabs:
            tab_name += f"_{idx}"

        try:
             tab_frame = self.tab_control.add(tab_name)
        except ValueError:
             tab_name += f"_v{idx}"
             tab_frame = self.tab_control.add(tab_name)

        self._tab_frames.append(tab_frame)

        fig = plt.Figure(figsize=FIGURE_SIZE, dpi=FIGURE_DPI)
        ax = fig.add_subplot(111)
        ax.set_title(f"{graph_title} - {basename}")

        renderer = self._renderers.get(command)
        if renderer:
            try:
                renderer(ax, results_df, aux)
            except Exception as e:
                logger.error(f"Error rendering chart for {basename}: {e}")
                ax.text(0.5, 0.5, f"Grafik çizim hatası:\n{e}", ha='center', va='center')

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=tab_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)
        self.tabs[tab_name] = canvas

        try:
            stem = os.path.splitext(basename)[0]
            file_out_dir = os.path.join(output_dir, stem)
            os.makedirs(file_out_dir, exist_ok=True)
            
            fig.savefig(
                os.path.join(file_out_dir, f"{command}_graph.png"),
                dpi=SAVE_DPI, format='png', bbox_inches='tight',
            )
            results_df.to_csv(
                os.path.join(file_out_dir, f"{command}_data.csv"), index=False,
            )
        except Exception as e:
             logger.error(f"Failed to save outputs for {basename}: {e}")

    def _analysis_complete(
        self,
        total: int,
        output_dir: str,
        error_files: list[str],
    ) -> None:
        self.analyze_btn.configure(state='normal')
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
    ctk.set_appearance_mode("System")
    ctk.set_default_color_theme("blue")
    root = ctk.CTk()
    FoldXAnalyzerGUI(root)
    root.mainloop()

# foldx_gui_analyzer_v8.py
# FoldX v8: TOPLU ANALİZ + UI TAM DÜZELTME
# Sekmeli Grafik, Otomatik Algılama, PNG Kaydet

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
from datetime import datetime
import re

class FoldXAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("FoldX Analyzer v8 - Toplu + UI Sabit")
        self.root.geometry("1200x900")
        self.root.configure(bg='#f4f6f9')

        style = ttk.Style()
        style.theme_use('clam')

        # === ÜST KISIM ===
        top_frame = tk.Frame(root, bg='#f4f6f9')
        top_frame.pack(pady=10, fill='x', padx=20)

        # Dosya / Klasör
        tk.Label(top_frame, text="fxout Dosyası/Klasörü:", font=("Helvetica", 10), bg='#f4f6f9').pack(anchor='w')
        file_sub = tk.Frame(top_frame, bg='#f4f6f9')
        file_sub.pack(fill='x', pady=2)
        self.file_path = tk.StringVar()
        tk.Entry(file_sub, textvariable=self.file_path, width=70, state='readonly').pack(side='left', padx=(0,5))
        ttk.Button(file_sub, text="Tek Dosya", command=self.browse_file).pack(side='left', padx=2)
        ttk.Button(file_sub, text="Toplu Klasör", command=self.browse_folder).pack(side='left')

        # Mod
        mode_frame = tk.Frame(top_frame, bg='#f4f6f9')
        mode_frame.pack(fill='x', pady=8)
        tk.Label(mode_frame, text="Mod:", font=("Helvetica", 10), bg='#f4f6f9').pack(side='left')
        self.mode_var = tk.StringVar(value="auto")
        ttk.Radiobutton(mode_frame, text="Otomatik", variable=self.mode_var, value="auto", command=self.toggle_manual).pack(side='left', padx=(10,5))
        ttk.Radiobutton(mode_frame, text="Manuel", variable=self.mode_var, value="manual", command=self.toggle_manual).pack(side='left', padx=5)
        self.manual_combo = ttk.Combobox(mode_frame, values=["PositionScan", "RepairPDB", "BuildModel", "AnalyseComplex", "Stability", "Pssm", "RnaScan"], state="disabled", width=15)
        self.manual_combo.pack(side='left', padx=5)
        self.manual_combo.set("PositionScan")

        # Başlık
        tk.Label(top_frame, text="Grafik Başlığı:", font=("Helvetica", 10), bg='#f4f6f9').pack(anchor='w', pady=(10,2))
        self.title_var = tk.StringVar(value="FoldX Analizi")
        tk.Entry(top_frame, textvariable=self.title_var, width=70).pack(pady=2, fill='x')

        # Çıktı
        out_sub = tk.Frame(top_frame, bg='#f4f6f9')
        out_sub.pack(fill='x', pady=8)
        tk.Label(out_sub, text="Çıktı Klasörü:", font=("Helvetica", 10), bg='#f4f6f9').pack(side='left')
        self.out_dir = tk.StringVar(value=os.path.join(os.getcwd(), "output"))
        tk.Entry(out_sub, textvariable=self.out_dir, width=60).pack(side='left', padx=5)
        ttk.Button(out_sub, text="Değiştir", command=self.browse_output).pack(side='left')

        # Butonlar
        btn_frame = tk.Frame(top_frame, bg='#f4f6f9')
        btn_frame.pack(pady=15)
        self.analyze_btn = ttk.Button(btn_frame, text="Analizi Başlat", command=self.run_analysis)
        self.analyze_btn.pack()

        # İlerleme
        self.progress = ttk.Progressbar(top_frame, mode='determinate', maximum=100)
        self.progress.pack(fill='x', pady=5)

        # Durum
        self.status = tk.Label(root, text="Hazır.", relief='sunken', anchor='w', bg='#ecf0f1', fg='#2c3e50')
        self.status.pack(side='bottom', fill='x')

        # === GRAFİK SEKMELERİ ===
        self.tab_control = ttk.Notebook(root)
        self.tab_control.pack(fill='both', expand=True, padx=20, pady=10)

        self.tabs = {}  # {tab_name: canvas}

    def toggle_manual(self):
        self.manual_combo.config(state="readonly" if self.mode_var.get() == "manual" else "disabled")

    def browse_file(self):
        file = filedialog.askopenfilename(filetypes=[("FoldX Output", "*.fxout"), ("All", "*.*")])
        if file:
            self.file_path.set(file)

    def browse_folder(self):
        folder = filedialog.askdirectory()
        if folder:
            self.file_path.set(folder)

    def browse_output(self):
        dir_path = filedialog.askdirectory()
        if dir_path:
            self.out_dir.set(dir_path)

    def find_header_and_read(self, file_path):
        encodings = ['utf-8', 'latin-1', 'cp1254']
        last_error = None
        for enc in encodings:
            try:
                with open(file_path, 'r', encoding=enc) as f:
                    lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.strip().startswith('Pdb\t'):
                        df = pd.read_csv(file_path, sep='\t', skiprows=i, encoding=enc)
                        return df, enc, i, lines[:i]
            except Exception as e:
                # Konsola yaz, GUI'yi bozmadan devam et
                print(f"[find_header_and_read] {file_path} için {enc} kodlamasında hata: {e}")
                last_error = e
                continue
        raise ValueError(f"Dosya okunamadı: {file_path}. Son hata: {last_error}")

    def detect_command(self, df, header_lines):
        cols = set(df.columns.str.lower())
        pdb_sample = str(df['Pdb'].iloc[0]).lower()

        if 'interaction energy' in cols and any('wt_' in str(p).lower() for p in df['Pdb']):
            return "PositionScan"
        if 'total energy' in cols and re.search(r'repair_\d+', pdb_sample):
            return "RepairPDB"
        if 'average' in cols or 'raw' in os.path.basename(self.file_path.get()).lower():
            return "BuildModel"
        if all(x in cols for x in ['backbone', 'sidechain', 'interaction']):
            return "AnalyseComplex"
        if 'total energy' in cols and len(df) == 1:
            return "Stability"
        if 'pssm' in ' '.join(header_lines).lower():
            return "Pssm"
        if any(base in pdb_sample for base in ['a', 'u', 'g', 'c']) and 'rna' in ' '.join(header_lines).lower():
            return "RnaScan"
        return "Bilinmiyor"

    def _first_existing_column(self, df, candidates):
        cols_lower = {c.lower(): c for c in df.columns}
        for cand in candidates:
            if cand.lower() in cols_lower:
                return cols_lower[cand.lower()]
        return None

    def _to_numeric_safe(self, series):
        return pd.to_numeric(series, errors='coerce')

    def _make_text_tab(self, tab_name, text):
        tab = ttk.Frame(self.tab_control)
        self.tab_control.add(tab, text=tab_name)
        lbl = tk.Label(tab, text=text, justify='left', anchor='nw')
        lbl.pack(fill='both', expand=True, padx=12, pady=12)
        return tab

    def _process_positionscan(self, df):
        results = []
        df = df.copy()
        df['run'] = np.nan
        df['replicate'] = np.nan
        df['is_wt'] = df['Pdb'].astype(str).str.contains('WT_', na=False)
        prev_run = prev_rep = None
        for i, pdb in enumerate(df['Pdb'].astype(str)):
            if df.at[i, 'is_wt']:
                df.at[i, 'run'], df.at[i, 'replicate'] = prev_run, prev_rep
            else:
                try:
                    parts = pdb.split('_')
                    run = int(parts[-2])
                    rep = int(parts[-1].split('.')[0])
                    df.at[i, 'run'], df.at[i, 'replicate'] = run, rep
                    prev_run, prev_rep = run, rep
                except Exception as e:
                    print(f"[process_positionscan] Pdb adı parse edilemedi ({pdb}): {e}")
                    continue

        energy_col = self._first_existing_column(df, ['Interaction Energy', 'interaction energy'])
        if not energy_col:
            print("[process_positionscan] 'Interaction Energy' kolonu bulunamadı.")
            return results, None

        amino_acids = ['A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        for run in df['run'].dropna().unique():
            for rep in df['replicate'].dropna().unique():
                mut = df[(df['run'] == run) & (df['replicate'] == rep) & (~df['is_wt'])]
                wt = df[(df['run'] == run) & (df['replicate'] == rep) & (df['is_wt'])]
                if not mut.empty and not wt.empty:
                    try:
                        run_int = int(run)
                        if 1 <= run_int <= len(amino_acids):
                            diff = self._to_numeric_safe(mut[energy_col]).iloc[0] - self._to_numeric_safe(wt[energy_col]).iloc[0]
                            aa = amino_acids[run_int - 1]
                            results.append({'amino_acid': aa, 'energy': float(diff)})
                        else:
                            print(f"[process_positionscan] Geçersiz run indexi: {run_int}")
                    except Exception as e:
                        print(f"[process_positionscan] run/aa eşleme hatası (run={run}): {e}")
        return results, amino_acids

    def _process_repairpdb(self, df):
        results = []
        df = df.copy()
        model_col = 'model'
        df[model_col] = df['Pdb'].astype(str).str.extract(r'_(\d+)').astype(float)
        energy_col = self._first_existing_column(df, ['total energy', 'Total Energy'])
        if not energy_col:
            print("[process_repairpdb] 'total energy' kolonu bulunamadı.")
            return results
        for _, row in df.iterrows():
            if pd.isna(row[model_col]):
                print(f"[process_repairpdb] Model numarası bulunamadı: {row.get('Pdb')}")
                continue
            energy_val = self._to_numeric_safe(pd.Series([row[energy_col]])).iloc[0]
            if pd.isna(energy_val):
                continue
            results.append({'model': int(row[model_col]), 'energy': float(energy_val)})
        return results

    def _process_buildmodel(self, df):
        # BuildModel/average outputs vary; try to extract a sensible energy series.
        df = df.copy()
        energy_col = self._first_existing_column(df, ['average', 'Average', 'total energy', 'Total Energy', 'Energy', 'energy'])
        if not energy_col:
            print("[process_buildmodel] Enerji kolonu bulunamadı (average/total energy/energy).")
            return []
        energies = self._to_numeric_safe(df[energy_col])
        results = []
        for pdb, val in zip(df.get('Pdb', pd.Series(range(len(df)))), energies):
            if pd.isna(val):
                continue
            results.append({'pdb': str(pdb), 'energy': float(val)})
        return results

    def _process_analysecomplex(self, df):
        # Expect backbone/sidechain/interaction-like columns; create per-Pdb component rows.
        df = df.copy()
        back_col = self._first_existing_column(df, ['backbone', 'Backbone'])
        side_col = self._first_existing_column(df, ['sidechain', 'Sidechain', 'side chain', 'Side Chain'])
        inter_col = self._first_existing_column(df, ['interaction', 'Interaction'])
        if not (back_col or side_col or inter_col):
            print("[process_analysecomplex] backbone/sidechain/interaction kolonları bulunamadı.")
            return []
        results = []
        for _, row in df.iterrows():
            rec = {'pdb': str(row.get('Pdb', ''))}
            if back_col:
                rec['backbone'] = float(self._to_numeric_safe(pd.Series([row[back_col]])).iloc[0]) if pd.notna(row[back_col]) else np.nan
            if side_col:
                rec['sidechain'] = float(self._to_numeric_safe(pd.Series([row[side_col]])).iloc[0]) if pd.notna(row[side_col]) else np.nan
            if inter_col:
                rec['interaction'] = float(self._to_numeric_safe(pd.Series([row[inter_col]])).iloc[0]) if pd.notna(row[inter_col]) else np.nan
            results.append(rec)
        # Drop rows that are all-nan for components
        df_res = pd.DataFrame(results)
        comp_cols = [c for c in ['backbone', 'sidechain', 'interaction'] if c in df_res.columns]
        df_res = df_res.dropna(subset=comp_cols, how='all')
        return df_res.to_dict('records')

    def _process_stability(self, df):
        df = df.copy()
        energy_col = self._first_existing_column(df, ['total energy', 'Total Energy', 'Energy', 'energy'])
        if not energy_col or len(df) == 0:
            return []
        val = self._to_numeric_safe(df[energy_col]).iloc[0]
        if pd.isna(val):
            return []
        pdb_val = str(df.get('Pdb', pd.Series([''])).iloc[0]) if 'Pdb' in df.columns else ''
        return [{'pdb': pdb_val, 'energy': float(val)}]

    def _process_pssm(self, df, header_lines):
        # PSSM formats differ; best effort:
        # - If it has 'Position' and 20 AA columns => heatmap
        # - Else store the table as-is for CSV export and show basic line plot of numeric columns.
        df = df.copy()
        df.columns = [str(c) for c in df.columns]
        numeric_df = df.apply(pd.to_numeric, errors='coerce')
        if numeric_df.notna().sum().sum() == 0:
            return []
        # Save whole numeric table (plus non-numeric) as records; plotting will decide later.
        out = df.fillna('').to_dict('records')
        return out

    def _process_rnascan(self, df):
        # Similar to PositionScan but bases A,U,G,C (1..4). Use same energy diff logic.
        results = []
        df = df.copy()
        df['run'] = np.nan
        df['replicate'] = np.nan
        df['is_wt'] = df['Pdb'].astype(str).str.contains('WT_', na=False)
        prev_run = prev_rep = None
        for i, pdb in enumerate(df['Pdb'].astype(str)):
            if df.at[i, 'is_wt']:
                df.at[i, 'run'], df.at[i, 'replicate'] = prev_run, prev_rep
            else:
                try:
                    parts = pdb.split('_')
                    run = int(parts[-2])
                    rep = int(parts[-1].split('.')[0])
                    df.at[i, 'run'], df.at[i, 'replicate'] = run, rep
                    prev_run, prev_rep = run, rep
                except Exception as e:
                    print(f"[process_rnascan] Pdb adı parse edilemedi ({pdb}): {e}")
                    continue

        energy_col = self._first_existing_column(df, ['Interaction Energy', 'interaction energy'])
        if not energy_col:
            print("[process_rnascan] 'Interaction Energy' kolonu bulunamadı.")
            return results, None

        bases = ['A', 'U', 'G', 'C']
        for run in df['run'].dropna().unique():
            for rep in df['replicate'].dropna().unique():
                mut = df[(df['run'] == run) & (df['replicate'] == rep) & (~df['is_wt'])]
                wt = df[(df['run'] == run) & (df['replicate'] == rep) & (df['is_wt'])]
                if not mut.empty and not wt.empty:
                    try:
                        run_int = int(run)
                        if 1 <= run_int <= len(bases):
                            diff = self._to_numeric_safe(mut[energy_col]).iloc[0] - self._to_numeric_safe(wt[energy_col]).iloc[0]
                            base = bases[run_int - 1]
                            results.append({'base': base, 'energy': float(diff)})
                        else:
                            print(f"[process_rnascan] Geçersiz run indexi: {run_int}")
                    except Exception as e:
                        print(f"[process_rnascan] run/base eşleme hatası (run={run}): {e}")
        return results, bases

    def run_analysis(self):
        input_path = self.file_path.get()
        output_dir = self.out_dir.get()
        os.makedirs(output_dir, exist_ok=True)

        if not input_path:
            messagebox.showerror("Hata", "Lütfen dosya/klasör seçin!")
            return

        self.analyze_btn.config(state='disabled')
        self.progress['value'] = 0

        files = []
        if os.path.isfile(input_path):
            files = [input_path]
        elif os.path.isdir(input_path):
            files = [os.path.join(input_path, f) for f in os.listdir(input_path) if f.endswith('.fxout')]
        else:
            messagebox.showerror("Hata", "Geçersiz yol!")
            return

        total = len(files)
        if total == 0:
            messagebox.showwarning("Uyarı", "Seçilen klasörde .fxout dosyası bulunamadı.")
            self.analyze_btn.config(state='normal')
            return
        self.progress['maximum'] = total

        for idx, file in enumerate(files):
            self.status.config(text=f"İşleniyor: {os.path.basename(file)} ({idx+1}/{total})")
            self.root.update_idletasks()

            try:
                df, _, _, header_lines = self.find_header_and_read(file)
                detected = self.detect_command(df, header_lines)
                if detected == "Bilinmiyor" and self.mode_var.get() == "auto":
                    messagebox.showwarning(
                        "Komut Tespit Edilemedi",
                        f"{os.path.basename(file)} dosyası için FoldX komutu otomatik tespit edilemedi.\n"
                        f"Lütfen manuel modu kullanarak uygun komutu seçin."
                    )
                    continue
                mode = self.mode_var.get()
                if mode == "manual":
                    detected = self.manual_combo.get()

                # === SONUÇ İŞLE ===
                results = []
                aux = None
                if detected == "PositionScan":
                    results, aux = self._process_positionscan(df)

                elif detected == "RepairPDB":
                    results = self._process_repairpdb(df)

                elif detected == "BuildModel":
                    results = self._process_buildmodel(df)

                elif detected == "AnalyseComplex":
                    results = self._process_analysecomplex(df)

                elif detected == "Stability":
                    results = self._process_stability(df)

                elif detected == "Pssm":
                    results = self._process_pssm(df, header_lines)

                elif detected == "RnaScan":
                    results, aux = self._process_rnascan(df)

                if not results:
                    self.status.config(text=f"{os.path.basename(file)} için geçerli sonuç elde edilemedi.")
                    continue

                # === GRAFİK TAB ===
                tab_name = os.path.basename(file)[:20]
                if tab_name in self.tabs:
                    tab_name += f"_{idx}"
                tab = ttk.Frame(self.tab_control)
                self.tab_control.add(tab, text=tab_name)

                fig = plt.Figure(figsize=(10, 6), dpi=100)
                ax = fig.add_subplot(111)

                if detected == "PositionScan" and results:
                    stats = pd.DataFrame(results).groupby('amino_acid')['energy'].agg(['mean','std','count']).reset_index()
                    stats['sem'] = stats['std']/np.sqrt(stats['count'])
                    amino_acids = aux if aux else ['A','R','D','N','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
                    order = {aa:i for i,aa in enumerate(amino_acids)}
                    stats = stats.sort_values(by='amino_acid', key=lambda x: x.map(order))
                    bars = ax.bar(stats['amino_acid'], stats['mean'], yerr=stats['sem'], capsize=5, color='lightgray', edgecolor='black')
                    for bar, val in zip(bars, stats['mean']):
                        bar.set_color('green' if val < -0.5 else 'red' if val > 0.5 else 'gray')
                    ax.axhline(0, color='black', lw=0.8)
                    ax.set_ylabel('ΔΔG (kcal/mol)')

                elif detected == "RepairPDB" and results:
                    df_res = pd.DataFrame(results)
                    ax.plot(df_res['model'], df_res['energy'], 'o-', color='blue')
                    ax.set_xlabel('Model'); ax.set_ylabel('Total Energy'); ax.grid(True, alpha=0.3)

                elif detected == "BuildModel" and results:
                    df_res = pd.DataFrame(results)
                    # Use order as provided; show first N if too large
                    max_n = 60
                    if len(df_res) > max_n:
                        df_res = df_res.iloc[:max_n].copy()
                        ax.set_title(ax.get_title() + f" (ilk {max_n} satır)")
                    ax.plot(range(1, len(df_res) + 1), df_res['energy'], 'o-', color='purple')
                    ax.set_xlabel('Sıra')
                    ax.set_ylabel('Energy')
                    ax.grid(True, alpha=0.3)

                elif detected == "AnalyseComplex" and results:
                    df_res = pd.DataFrame(results)
                    max_n = 40
                    if len(df_res) > max_n:
                        df_res = df_res.iloc[:max_n].copy()
                        ax.set_title(ax.get_title() + f" (ilk {max_n} satır)")
                    x = np.arange(len(df_res))
                    bottom = np.zeros(len(df_res))
                    for comp, color in [('backbone', '#1f77b4'), ('sidechain', '#ff7f0e'), ('interaction', '#2ca02c')]:
                        if comp in df_res.columns:
                            vals = self._to_numeric_safe(df_res[comp]).fillna(0).to_numpy()
                            ax.bar(x, vals, bottom=bottom, label=comp, color=color, alpha=0.85)
                            bottom = bottom + vals
                    labels = df_res['pdb'].astype(str).tolist() if 'pdb' in df_res.columns else [str(i) for i in x]
                    ax.set_xticks(x)
                    ax.set_xticklabels(labels, rotation=90, fontsize=8)
                    ax.set_ylabel('Energy components (stacked)')
                    ax.legend()
                    ax.grid(True, axis='y', alpha=0.2)

                elif detected == "Stability" and results:
                    df_res = pd.DataFrame(results)
                    ax.bar(['Stability'], df_res['energy'].astype(float), color='teal')
                    ax.set_ylabel('Total Energy')
                    ax.grid(True, axis='y', alpha=0.3)

                elif detected == "Pssm" and results:
                    # Best-effort plot: heatmap of numeric matrix (drop non-numeric columns)
                    df_all = pd.DataFrame(results)
                    num = df_all.apply(pd.to_numeric, errors='coerce')
                    # Drop columns with too many NaNs
                    keep_cols = [c for c in num.columns if num[c].notna().sum() >= max(1, int(0.3 * len(num)))]
                    num = num[keep_cols].dropna(how='all')
                    if num.empty or num.shape[1] < 2:
                        ax.text(0.01, 0.99, "PSSM grafik üretilemedi (sayısal veri bulunamadı).", transform=ax.transAxes,
                                va='top', ha='left')
                        ax.axis('off')
                    else:
                        im = ax.imshow(num.to_numpy(), aspect='auto', interpolation='nearest', cmap='viridis')
                        ax.set_xlabel('Kolonlar')
                        ax.set_ylabel('Satırlar')
                        ax.set_title('PSSM (heatmap)')
                        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

                elif detected == "RnaScan" and results:
                    df_res = pd.DataFrame(results).groupby('base')['energy'].agg(['mean', 'std', 'count']).reset_index()
                    df_res['sem'] = df_res['std'] / np.sqrt(df_res['count'])
                    bases = aux if aux else ['A', 'U', 'G', 'C']
                    order = {b: i for i, b in enumerate(bases)}
                    df_res = df_res.sort_values(by='base', key=lambda x: x.map(order))
                    ax.bar(df_res['base'], df_res['mean'], yerr=df_res['sem'], capsize=5, color='lightgray', edgecolor='black')
                    ax.axhline(0, color='black', lw=0.8)
                    ax.set_ylabel('ΔΔG (kcal/mol)')

                title = self.title_var.get() or f"{detected}"
                ax.set_title(f"{title} - {os.path.basename(file)}")

                fig.tight_layout()
                canvas = FigureCanvasTkAgg(fig, master=tab)
                canvas.draw()
                canvas.get_tk_widget().pack(fill='both', expand=True)
                self.tabs[tab_name] = canvas

                # === KAYDET ===
                file_out_dir = os.path.join(output_dir, os.path.splitext(os.path.basename(file))[0])
                os.makedirs(file_out_dir, exist_ok=True)
                png_path = os.path.join(file_out_dir, f"{detected}_graph.png")
                fig.savefig(png_path, dpi=300, format='png', bbox_inches='tight')
                csv_path = os.path.join(file_out_dir, f"{detected}_data.csv")
                pd.DataFrame(results).to_csv(csv_path, index=False)

            except Exception as e:
                print(f"[run_analysis] {file} işlenirken hata: {e}")
                self.status.config(text=f"Hata: {os.path.basename(file)}")
                messagebox.showerror("Analiz Hatası",
                                     f"{os.path.basename(file)} dosyası işlenirken hata oluştu:\n{e}")
                continue
            finally:
                self.progress['value'] = idx + 1
                self.root.update_idletasks()

        self.status.config(text=f"Toplam {total} dosya analiz edildi!")
        self.analyze_btn.config(state='normal')
        messagebox.showinfo("Başarılı!", f"{total} dosya analiz edildi!\nÇıktılar: {output_dir}")

# === ÇALIŞTIR ===
if __name__ == "__main__":
    root = tk.Tk()
    app = FoldXAnalyzerGUI(root)
    root.mainloop()
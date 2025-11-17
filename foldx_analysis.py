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
        for enc in encodings:
            try:
                with open(file_path, 'r', encoding=enc) as f:
                    lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.strip().startswith('Pdb\t'):
                        df = pd.read_csv(file_path, sep='\t', skiprows=i, encoding=enc)
                        return df, enc, i, lines[:i]
            except:
                continue
        raise ValueError("Dosya okunamadı.")

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
        self.progress['maximum'] = total

        for idx, file in enumerate(files):
            self.status.config(text=f"İşleniyor: {os.path.basename(file)} ({idx+1}/{total})")
            self.root.update_idletasks()

            try:
                df, _, _, header_lines = self.find_header_and_read(file)
                detected = self.detect_command(df, header_lines)
                mode = self.mode_var.get()
                if mode == "manual":
                    detected = self.manual_combo.get()

                # === SONUÇ İŞLE ===
                results = []
                if detected == "PositionScan":
                    df['run'] = np.nan; df['replicate'] = np.nan; df['is_wt'] = df['Pdb'].str.contains('WT_', na=False)
                    prev_run = prev_rep = None
                    for i, pdb in enumerate(df['Pdb']):
                        if df.at[i, 'is_wt']:
                            df.at[i, 'run'], df.at[i, 'replicate'] = prev_run, prev_rep
                        else:
                            try:
                                parts = pdb.split('_')
                                run = int(parts[-2]); rep = int(parts[-1].split('.')[0])
                                df.at[i, 'run'], df.at[i, 'replicate'] = run, rep
                                prev_run, prev_rep = run, rep
                            except: pass
                    amino_acids = ['A','R','D','N','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
                    for run in df['run'].dropna().unique():
                        for rep in df['replicate'].dropna().unique():
                            mut = df[(df['run']==run)&(df['replicate']==rep)&(~df['is_wt'])]
                            wt = df[(df['run']==run)&(df['replicate']==rep)&(df['is_wt'])]
                            if not mut.empty and not wt.empty:
                                diff = mut['Interaction Energy'].iloc[0] - wt['Interaction Energy'].iloc[0]
                                aa = amino_acids[int(run)-1]
                                results.append({'amino_acid': aa, 'energy': diff})

                elif detected == "RepairPDB":
                    df['model'] = df['Pdb'].str.extract(r'_(\d+)').astype(float)
                    for _, row in df.iterrows():
                        results.append({'model': int(row['model']), 'energy': row['total energy']})

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
                self.status.config(text=f"Hata: {os.path.basename(file)}")
                continue

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
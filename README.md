# FoldX Analyzer

**Türkçe** | [English](#english)

---

## Türkçe

FoldX çıktı dosyalarını (`.fxout`) otomatik olarak analiz eden, grafiklerini GUI içinde gösteren ve PNG + CSV olarak kaydeden masaüstü uygulaması.

### Özellikler

- 7 FoldX komutunu **otomatik algılar** (veya manuel seçim)
- **Toplu analiz** — klasör seçerek yüzlerce `.fxout` dosyasını aynı anda işler
- Her dosya ayrı **sekmede grafik** olarak gösterilir
- ΔΔG bar plot: yeşil = stabilizasyon, kırmızı = destabilizasyon
- **300 DPI PNG + CSV** çıktılar otomatik isimlendirilip kaydedilir
- Analiz arka planda çalışır, arayüz donmaz
- PDB yapı analizi: CYS komşu arama, disülfür mesafesi hesabı (Biopython)
- Protein-ligand etkileşim analizi ve PyMOL görselleştirme (PLIP)
- %100 Python, açık kaynak, ücretsiz

### Desteklenen FoldX Komutları

| Komut | Grafik Türü | Açıklama |
|---|---|---|
| PositionScan | Bar plot | Amino asit mutasyon enerji etkileri (ΔΔG) |
| RepairPDB | Çizgi grafik | Model başına total enerji |
| BuildModel | Çizgi grafik | Sıralı model enerjileri |
| AnalyseComplex | Yığılmış bar | Backbone / sidechain / interaction bileşenleri |
| Stability | Tek bar | Tek stabilite ölçümü |
| Pssm | Heatmap | Pozisyon özgül skor matrisi |
| RnaScan | Bar plot | RNA baz mutasyon etkileri (A/U/G/C) |

### Kurulum

```bash
git clone https://github.com/ohatipoglu/foldx-analyzer.git
cd foldx-analyzer
python -m venv venv
venv\Scripts\activate        # Windows
# source venv/bin/activate   # macOS / Linux
pip install -r requirements.txt
```

> **Not:** PLIP görselleştirme için sisteme PyMOL kurulu olmalı ve `pymol` komutu PATH'te erişilebilir olmalıdır.

### Kullanım

```bash
python combined_gui.py
```

İki pencere açılır:

1. **FoldX Analyzer** — `.fxout` dosyası veya klasörü seçin → "Analizi Başlat"
2. **PDB / PLIP Analyzer** — `.pdb` dosyası seçin → Biopython veya PLIP analizini çalıştırın

Çıktılar `output/<dosya_adı>/` klasörüne kaydedilir.

### Proje Yapısı

```
foldx-analyzer/
├── combined_gui.py      # Giriş noktası — her iki pencereyi başlatır
├── foldx_analysis.py    # FoldX Analyzer GUI
├── pdb_analyzer.py      # PDB / PLIP Analyzer GUI
├── foldx_core.py        # Saf veri işleme katmanı (UI bağımsız)
├── constants.py         # Uygulama geneli sabitler
├── data/                # Örnek .fxout dosyaları
└── output/              # Üretilen grafik ve CSV çıktıları
```

### Gereksinimler

- Python ≥ 3.12
- pandas, numpy, matplotlib
- biopython, plip (PDB analizi için)
- PyMOL (PLIP görselleştirme için, opsiyonel)

---

## English

<a name="english"></a>

A desktop GUI application that automatically analyzes FoldX output files (`.fxout`), displays charts inside the interface, and saves PNG + CSV outputs.

### Features

- **Auto-detects** all 7 FoldX commands (or manual override)
- **Batch processing** — select a folder to analyze hundreds of `.fxout` files at once
- Each file is shown as a chart in its own **tab**
- ΔΔG bar plots: green = stabilizing, red = destabilizing
- **300 DPI PNG + CSV** outputs saved and named automatically
- Analysis runs in a background thread — UI stays responsive
- PDB structure analysis: CYS neighbor search, disulfide distance calculation (Biopython)
- Protein–ligand interaction analysis and PyMOL visualization (PLIP)
- 100% Python, open source, free

### Supported FoldX Commands

| Command | Chart Type | Description |
|---|---|---|
| PositionScan | Bar plot | Amino acid mutation energy effects (ΔΔG) |
| RepairPDB | Line chart | Total energy per model |
| BuildModel | Line chart | Sequential model energies |
| AnalyseComplex | Stacked bar | Backbone / sidechain / interaction components |
| Stability | Single bar | Single stability measurement |
| Pssm | Heatmap | Position-specific scoring matrix |
| RnaScan | Bar plot | RNA base mutation effects (A/U/G/C) |

### Installation

```bash
git clone https://github.com/ohatipoglu/foldx-analyzer.git
cd foldx-analyzer
python -m venv venv
venv\Scripts\activate        # Windows
# source venv/bin/activate   # macOS / Linux
pip install -r requirements.txt
```

> **Note:** PLIP visualization requires PyMOL to be installed and the `pymol` command available on PATH.

### Usage

```bash
python combined_gui.py
```

Two windows open:

1. **FoldX Analyzer** — select a `.fxout` file or folder → click "Analizi Başlat"
2. **PDB / PLIP Analyzer** — select a `.pdb` file → run Biopython or PLIP analysis

Outputs are saved to `output/<filename>/`.

### Project Structure

```
foldx-analyzer/
├── combined_gui.py      # Entry point — launches both windows
├── foldx_analysis.py    # FoldX Analyzer GUI
├── pdb_analyzer.py      # PDB / PLIP Analyzer GUI
├── foldx_core.py        # Pure data-processing layer (no UI dependency)
├── constants.py         # Application-wide constants
├── data/                # Sample .fxout files
└── output/              # Generated charts and CSV outputs
```

### Requirements

- Python ≥ 3.12
- pandas, numpy, matplotlib
- biopython, plip (for PDB analysis)
- PyMOL (for PLIP visualization, optional)

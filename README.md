# Bioinformatics Analyzer Suite (FoldX & PLIP)

**Türkçe** | [English](#english)

> **📢 Latest:** Version 0.2.0 includes comprehensive code quality improvements, full English language support, type hints, and a detailed [macOS Setup Guide](MACOS_SETUP_GUIDE.md).

---

## Türkçe

**Bioinformatics Analyzer Suite**, protein yapılarını ve enerjilerini analiz etmek için geliştirilmiş kapsamlı bir masaüstü uygulamasıdır. FoldX çıktılarını otomatik olarak işler, protein-ligand etkileşimlerini (PLIP) analiz eder ve yapısal biyoinformatik hesaplamaları yapar. Modern bir arayüz (CustomTkinter) ile donatılmış olan uygulama, hem bireysel hem de toplu işlemleri destekler.

### Temel Özellikler

#### 1. FoldX Analizörü
- **Otomatik Komut Algılama:** 7 farklı FoldX komutunun (`PositionScan`, `RepairPDB`, `BuildModel`, `AnalyseComplex`, `Stability`, `PSSM`, `RnaScan`) çıktılarını otomatik olarak tanır.
- **Toplu İşlem:** Bir klasör dolusu `.fxout` dosyasını tek seferde analiz eder ve sonuçları sekmeli görünümde sunar.
- **Görselleştirme ve Çıktı:** Analiz sonuçlarını grafiklere döker, yüksek çözünürlüklü PNG ve CSV formatında `output` klasörüne kaydeder.

#### 2. PDB & Yapısal Analiz (Biopython)
- **CYS Arama:** PDB dosyasındaki tüm Sistein (CYS) kalıntılarını listeler.
- **Mesafe Hesaplama:** Belirlenen kalıntılar arasındaki (örneğin disülfür köprüleri) mesafeyi otomatik hesaplar.
- **Komşu Analizi:** Hedef kalıntıların çevresindeki etkileşim alanlarını belirler.

#### 3. PLIP (Protein-Ligand Interaction Profiler) Entegrasyonu
- **Canlı Analiz (3D):** Bir PDB dosyasını PLIP ile analiz eder ve sonuçları doğrudan interaktif **3D PyMOL** üzerinde görselleştirir.
- **Batch TXT İşleme:** PLIP tarafından üretilen metin raporlarını (`.txt`) toplu olarak okur ve analiz eder.
- **Word Raporu Dışa Aktarma:** Farklı mutasyonların veya dosyaların etkileşim sayılarını (Hidrojen bağı, Hidrofobik etkileşim, Tuz köprüsü vb.) karşılaştıran profesyonel bir **Microsoft Word (.docx)** raporu oluşturur.

#### 4. Kullanıcı Deneyimi ve Altyapı
- **Modern Arayüz:** CustomTkinter ile yan menü (sidebar) tabanlı, şık ve hızlı bir kullanım.
- **Yapılandırma (config.ini):** PyMOL yolu, grafik çözünürlükleri ve enerji eşik değerleri gibi ayarları kalıcı olarak kaydeder.
- **Arka Plan İşlemleri:** Uzun süren analizler sırasında arayüzün donmasını engelleyen multi-threading yapısı.
- **Detaylı Loglama:** Hata takibi için `app.log` dosyası.
- **Tip Güvenliği:** Tüm fonksiyonlar type hints ile dokümante edilmiştir.
- **Test Coverage:** 20+ unit test ile kritik fonksiyonlar güvence altındadır.

### Kurulum

**1. Gerekli Programlar:**
- Python ≥ 3.12
- **PyMOL:** Görselleştirmeler için gereklidir. Sistem `PATH`'ine eklenmiş olmalı veya `config.ini` dosyasında yolu belirtilmelidir.

**2. Bağımlılıkların Kurulumu:**

```bash
git clone https://github.com/ohatipoglu/foldx_analyzer.git
cd foldx_analyzer

# Sanal ortam oluşturup aktif ettikten sonra:
pip install .
```

**⚠️ ÖNEMLİ (Windows - PLIP Kurulumu):**
Windows üzerinde PLIP ve OpenBabel kurulumu `pip` ile hata verebilir. Conda kullanmanız önerilir:
```bash
conda install -c conda-forge openbabel plip
```

### 🍎 macOS Kullanıcıları

macOS üzerinde kurulum için detaylı talimatlar: **[MACOS_SETUP_GUIDE.md](MACOS_SETUP_GUIDE.md)**

**Hızlı Başlangıç (Anaconda ile):**
```bash
conda create -n foldx-env python=3.12 -y
conda activate foldx-env
conda install -c conda-forge pandas numpy matplotlib biopython openbabel plip pymol pytest -y
pip install customtkinter python-docx
pip install -e .
python foldx_analyzer_core/combined_gui.py
```

### Kullanım

Uygulamayı başlatmak için:
```bash
python foldx_analyzer_core/combined_gui.py
```

### Yapılandırma (`config.ini`)

Uygulama ilk açılışta varsayılan ayarlarla bir `config.ini` oluşturur. PyMOL yolunu buradan güncelleyebilirsiniz:
```ini
[Settings]
pymol_executable_path = C:/Program Files/PyMOL/PyMOL.exe
```

---

## English
<a name="english"></a>

**Bioinformatics Analyzer Suite** is a comprehensive desktop application developed for analyzing protein structures and energies. It automatically processes FoldX outputs, analyzes protein-ligand interactions (PLIP), and performs structural bioinformatics calculations. Equipped with a modern interface (CustomTkinter), the application supports both individual and batch processing.

### Key Features

#### 1. FoldX Analyzer
- **Auto Command Detection:** Automatically recognizes outputs from 7 different FoldX commands (`PositionScan`, `RepairPDB`, `BuildModel`, `AnalyseComplex`, `Stability`, `PSSM`, `RnaScan`).
- **Batch Processing:** Analyzes an entire folder of `.fxout` files at once and presents results in a tabbed view.
- **Visualization & Export:** Converts analysis results into charts and saves them in high-resolution PNG and CSV formats to the `output` folder.

#### 2. PDB & Structural Analysis (Biopython)
- **CYS Search:** Lists all Cysteine (CYS) residues in the PDB file.
- **Distance Calculation:** Automatically calculates distances between specified residues (e.g., disulfide bridges).
- **Neighbor Analysis:** Identifies interaction zones around target residues.

#### 3. PLIP (Protein-Ligand Interaction Profiler) Integration
- **Live Analysis (3D):** Analyzes a PDB file with PLIP and visualizes the results directly in interactive **3D PyMOL**.
- **Batch TXT Processing:** Batch reads and analyzes text reports (`.txt`) generated by PLIP.
- **Word Report Export:** Generates a professional **Microsoft Word (.docx)** report comparing interaction counts (Hydrogen bonds, Hydrophobic interactions, Salt bridges, etc.) across different mutations or files.

#### 4. User Experience & Infrastructure
- **Modern UI:** A sleek and fast sidebar-based interface built with CustomTkinter.
- **Configuration (config.ini):** Permanently saves settings such as PyMOL path, graph resolutions, and energy thresholds.
- **Background Processing:** Multi-threading structure prevents UI freezing during long-running analyses.
- **Detailed Logging:** `app.log` file for error tracking.
- **Type Safety:** All functions are documented with type hints for better IDE support.
- **Test Coverage:** 20+ unit tests cover critical functionality.

### Installation

**1. Prerequisites:**
- Python ≥ 3.12
- **PyMOL:** Required for visualizations. Must be added to system `PATH` or specified in `config.ini`.

**2. Install Dependencies:**

```bash
git clone https://github.com/ohatipoglu/foldx_analyzer.git
cd foldx_analyzer

# After creating and activating a virtual environment:
pip install .
```

**⚠️ IMPORTANT (Windows - PLIP Installation):**
Installing PLIP and OpenBabel on Windows via `pip` may fail. Conda is recommended:
```bash
conda install -c conda-forge openbabel plip
```

### 🍎 macOS Users

For detailed macOS installation instructions: **[MACOS_SETUP_GUIDE.md](MACOS_SETUP_GUIDE.md)**

**Quick Start (with Anaconda):**
```bash
conda create -n foldx-env python=3.12 -y
conda activate foldx-env
conda install -c conda-forge pandas numpy matplotlib biopython openbabel plip pymol pytest -y
pip install customtkinter python-docx
pip install -e .

> 💡 **Important:** First, install scientific packages with `conda`, then the rest with `pip`. Reversing the order can lead to environment instability.

python foldx_analyzer_core/combined_gui.py
```

### Usage

To launch the application:
```bash
python foldx_analyzer_core/combined_gui.py
```

### Configuration (`config.ini`)

The application creates a `config.ini` with default settings on the first run. You can update the PyMOL path here:
```ini
[Settings]
pymol_executable_path = C:/Program Files/PyMOL/PyMOL.exe
```

### Running Tests

```bash
python -m pytest tests/
```

**Test Coverage:**
- `test_foldx_core.py` - 20 tests for FoldX data processing
- `test_pdb_analyzer.py` - 15+ tests for PLIP parsing and structural analysis (requires Biopython)

---

## Development

### Code Quality Features (v0.2.0+)

- **Type Hints:** All public functions include type annotations for better IDE support
- **Docstrings:** Comprehensive documentation following Google style
- **Error Handling:** Specific exception types with detailed logging
- **Config Validation:** Automatic validation and correction of invalid settings

### Project Structure

```
foldx_analyzer/
├── combined_gui.py      # Main application entry point
├── foldx_analysis.py    # FoldX GUI module
├── foldx_core.py        # FoldX business logic (tested)
├── pdb_analyzer.py      # PLIP & Biopython module (tested)
├── constants.py         # Configuration and constants
├── tests/
│   ├── test_foldx_core.py
│   └── test_pdb_analyzer.py
├── MACOS_SETUP_GUIDE.md  # macOS installation guide
└── README.md
```

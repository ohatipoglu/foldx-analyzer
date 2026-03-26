# FoldX Analyzer

**Türkçe** | [English](#english)

---

## Türkçe

FoldX çıktı dosyalarını (`.fxout`) otomatik olarak analiz eden, grafiklerini modern bir arayüz (CustomTkinter) içinde gösteren ve PNG + CSV olarak kaydeden masaüstü uygulaması. 

Ayrıca PDB (Protein Data Bank) yapı analizi ve PLIP (Protein-Ligand Interaction Profiler) entegrasyonu ile protein-ligand etkileşimlerini PyMOL üzerinde görselleştirme yeteneklerine sahiptir.

### Özellikler

- **Modern Arayüz:** CustomTkinter ile oluşturulmuş, kullanıcı dostu ve şık bir görünüm.
- **7 FoldX Komutunu Algılama:** PositionScan, RepairPDB, BuildModel, AnalyseComplex vb. komutların çıktılarını otomatik tanır veya manuel seçime izin verir.
- **Toplu İşlem (Batch Processing):** Bir klasör seçerek yüzlerce `.fxout` dosyasını tek seferde analiz etme imkanı.
- **Özelleştirilebilir Yapılandırma:** `config.ini` dosyası üzerinden PyMOL yolu, grafik çözünürlükleri (DPI) ve enerji eşik değerlerini kolayca ayarlama.
- **Sekmeli Görünüm:** Her bir analiz sonucu ayrı bir sekmede, detaylı grafiklerle gösterilir.
- **Otomatik Çıktı Kaydı:** Yüksek çözünürlüklü (300 DPI) PNG grafikler ve işlenmiş CSV verileri `output` klasörüne otomatik kaydedilir.
- **Arka Planda Çalışma:** Uzun süren analizler arka planda yürütülür, böylece arayüz donmaz.
- **Hata Yönetimi ve Loglama:** Beklenmeyen hatalar için kullanıcı dostu uyarılar ve detaylı loglama (`app.log`).
- **Birim Testleri (Unit Tests):** Veri işleme mantığının doğruluğunu garanti eden test altyapısı.
- **Biopython & PLIP Entegrasyonu:** CYS komşu araması, disülfür mesafesi hesaplama ve protein-ligand etkileşim analizleri.

### Kurulum

**1. Gerekli Programların Yüklenmesi:**
- Python ≥ 3.12 yüklü olmalıdır.
- (İsteğe bağlı) PLIP görselleştirmeleri için PyMOL yüklü olmalı ve sistem `PATH`'ine eklenmiş olmalı veya `config.ini`'de yolu belirtilmelidir.

**2. Projenin İndirilmesi ve Bağımlılıkların Kurulumu:**

```bash
git clone https://github.com/ohatipoglu/foldx-analyzer.git
cd foldx-analyzer

# Sanal ortam oluşturma ve aktifleştirme
python -m venv venv
# Windows için:
venv\Scripts\activate
# macOS / Linux için:
# source venv/bin/activate

# Bağımlılıkları yükleme
pip install -r requirements.txt
# VEYA modern arayüz için (eğer requirements.txt'de yoksa):
pip install customtkinter pandas numpy matplotlib biopython plip pytest
```

### Yapılandırma (`config.ini`)

İlk çalıştırmada otomatik olarak bir `config.ini` dosyası oluşturulacaktır. Bu dosyayı düzenleyerek uygulama ayarlarını kendinize göre özelleştirebilirsiniz:

```ini
[Settings]
pymol_executable_path = pymol  # Veya C:/Program Files/PyMOL/PyMOL.exe

[Graph]
figure_dpi = 100
save_dpi = 300
stabilizing_threshold = -0.5
destabilizing_threshold = 0.5

[PDB]
neighbor_threshold = 6.0
disulfide_min_dist = 2.0
disulfide_max_dist = 2.2
```

### Kullanım

```bash
python combined_gui.py
```

Uygulama iki pencere açacaktır:
1. **FoldX Analyzer:** Klasör/dosya seçimi ve grafik analizleri için ana pencere.
2. **PDB / PLIP Analyzer:** PDB analizleri ve PyMOL entegrasyonu için yardımcı pencere.

### Testleri Çalıştırma

Kodun doğruluğunu kontrol etmek için:
```bash
python -m pytest tests/
```

---

## English
<a name="english"></a>

A desktop application that automatically analyzes FoldX output files (`.fxout`), displays charts within a modern interface (CustomTkinter), and saves them as PNG + CSV.

It also features PDB (Protein Data Bank) structure analysis and PLIP (Protein-Ligand Interaction Profiler) integration to visualize protein-ligand interactions in PyMOL.

### Features

- **Modern UI:** User-friendly and sleek interface built with CustomTkinter.
- **Auto-detects 7 FoldX Commands:** Automatically recognizes outputs from PositionScan, RepairPDB, BuildModel, AnalyseComplex, etc., or allows manual selection.
- **Batch Processing:** Analyze hundreds of `.fxout` files at once by selecting a folder.
- **Customizable Configuration:** Easily adjust PyMOL paths, graph resolutions (DPI), and energy thresholds via the `config.ini` file.
- **Tabbed View:** Each analysis result is displayed in a separate tab with detailed charts.
- **Automatic Output Saving:** High-resolution (300 DPI) PNG charts and processed CSV data are automatically saved to the `output` folder.
- **Background Processing:** Long-running analyses are executed in the background, keeping the UI responsive.
- **Error Handling & Logging:** User-friendly alerts for unexpected errors and detailed logging (`app.log`).
- **Unit Tests:** Test infrastructure ensuring the correctness of the data processing logic.
- **Biopython & PLIP Integration:** CYS neighbor search, disulfide distance calculation, and protein-ligand interaction analysis.

### Installation

**1. Prerequisites:**
- Python ≥ 3.12 must be installed.
- (Optional) PyMOL must be installed for PLIP visualizations and added to the system `PATH`, or its path must be specified in `config.ini`.

**2. Clone and Install Dependencies:**

```bash
git clone https://github.com/ohatipoglu/foldx-analyzer.git
cd foldx-analyzer

# Create and activate a virtual environment
python -m venv venv
# For Windows:
venv\Scripts\activate
# For macOS / Linux:
# source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
# OR for the modern UI (if not in requirements.txt):
pip install customtkinter pandas numpy matplotlib biopython plip pytest
```

### Configuration (`config.ini`)

A `config.ini` file will be created automatically on the first run. You can customize the application settings by editing this file:

```ini
[Settings]
pymol_executable_path = pymol  # Or C:/Program Files/PyMOL/PyMOL.exe

[Graph]
figure_dpi = 100
save_dpi = 300
stabilizing_threshold = -0.5
destabilizing_threshold = 0.5

[PDB]
neighbor_threshold = 6.0
disulfide_min_dist = 2.0
disulfide_max_dist = 2.2
```

### Usage

```bash
python combined_gui.py
```

The application will launch two windows:
1. **FoldX Analyzer:** Main window for folder/file selection and chart analysis.
2. **PDB / PLIP Analyzer:** Auxiliary window for PDB analysis and PyMOL integration.

### Running Tests

To verify the correctness of the code:
```bash
python -m pytest tests/
```

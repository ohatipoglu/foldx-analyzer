# macOS Kurulum ve Hata Ayıklama (Troubleshooting) Günlüğü

Bu doküman, FoldX ve PLIP entegrasyonu içeren `foldx_analyzer` projesinin macOS üzerinde (venv kullanılarak) kurulmaya çalışılırken karşılaşılan sorunları ve nihai Conda çözümünü belgelemektedir.

## Karşılaşılan Sorunlar ve Nedenleri

### Sorun 1: Poetry "ModuleOrPackageNotFoundError" Hatası
**Durum:** `pip install -e .` komutu çalıştırıldığında Poetry paket klasörünü bulamadı.
**Neden:** Python standartlarında proje ana dizini ile paket klasörü aynı isimde (tireli) olamaz. Ayrıca proje dosyaları ana dizine dağınıktı.
**Çözüm:** 1. `foldx_analyzer_core` adında alt tireli bir klasör oluşturuldu.
2. İçine `__init__.py` eklendi ve tüm `.py` dosyaları buraya taşındı.
3. `pyproject.toml` dosyası `packages = [{ include = "foldx_analyzer_core" }]` şeklinde güncellendi.

### Sorun 2: OpenBabel ve PLIP Derleme (SWIG) Hatası
**Durum:** Sanal ortamda (`venv`) `pip install plip` komutu çalıştırıldığında `Error: SWIG failed. Is Open Babel installed?` hatası alındı.
**Neden:** PLIP, OpenBabel'in sistemde C++ seviyesinde derlenmiş olmasını (SWIG ile bağlanmasını) bekler. macOS'te çip mimarisi (Intel vs Apple Silicon) ve Homebrew yollarının karmaşıklığı nedeniyle pip üzerinden derleme işlemi sürekli path (LDFLAGS, CPPFLAGS) hataları verdi.
**Çözüm:** Biyoinformatik araçları için standart `venv` kullanımından vazgeçilip, önceden derlenmiş (pre-built) binary'ler sunan **Conda** ekosistemine geçiş yapıldı.

### Sorun 3: Conda Hizmet Şartları (ToS) Engeli
**Durum:** `conda create` komutu "Terms of Service have not been accepted" hatası verdi.
**Neden:** Anaconda'nın yeni lisans/kullanım politikası gereği main ve r kanallarının ToS onayı istemesi.
**Çözüm:** Terminalden manuel olarak onay verildi:
`conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main`

### Sorun 4: PyMOL'un Conda-Forge'da Bulunamaması
**Durum:** `conda install -c conda-forge pymol` komutu "PackagesNotFoundError" hatası verdi.
**Neden:** PyMOL'un açık kaynaklı versiyonu conda-forge üzerinde `pymol` değil, `pymol-open-source` adıyla barındırılmaktadır.
**Çözüm:** Kurulum komutu `conda install -c conda-forge pymol-open-source` olarak güncellendi.

---

## Nihai ve Sorunsuz Çözüm (Conda)

Tüm sorunları aşan ve projeyi 2 dakika içinde ayağa kaldıran nihai komut dizisi şu şekildedir:

```bash
# 1. Miniconda Kurulumu (Eğer yoksa)
brew install --cask miniconda
conda init zsh

# 2. ToS Onayları (Yeni kurulumlarda gerekebilir)
conda tos accept --override-channels --channel [https://repo.anaconda.com/pkgs/main](https://repo.anaconda.com/pkgs/main)
conda tos accept --override-channels --channel [https://repo.anaconda.com/pkgs/r](https://repo.anaconda.com/pkgs/r)

# 3. Ortam Yaratma ve Aktivasyon
conda create -n foldx_env python=3.12 -y
conda activate foldx_env

# 4. Kritik Biyoinformatik Paketlerinin Kurulumu
conda install -c conda-forge openbabel plip pymol-open-source -y

# 5. Projenin Kurulumu
pip install -e .
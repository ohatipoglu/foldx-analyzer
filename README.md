# Otomatik .fxout Analiz GUI

**FoldX’in tüm komutlarını tek tıkla analiz eden, tamamen Türkçe GUI araç!**

PositionScan • RepairPDB • BuildModel • AnalyseComplex • Stability • Pssm • RnaScan  
→ Otomatik algılar → Grafik GUI içinde gösterir → PNG + CSV otomatik kaydeder

![demo](https://github.com/ohatipoglu/foldx-analyzer/raw/main/demo.gif)  


### Özellikler
- 7+ FoldX komutunu **otomatik algılar**  
- **Toplu analiz** (klasör seç, yüzlerce .fxout aynı anda)  
- Her dosya ayrı **sekmede grafik** gösterilir  
- ΔΔG bar plot (yeşil = stabilizasyon, kırmızı = destabilizasyon)  
- 300 DPI PNG + CSV çıktılar (otomatik isimlendirilir)  
- %100 Python, açık kaynak, ücretsiz

### Kurulum (30 saniye)

```bash
# 1. Projeyi indir
git clone https://github.com/ohatipoglu/foldx-analyzer.git

# 2. Klasöre gir
cd foldx-analyzer

# 3. Sanal ortam kur ve çalıştır (tek seferlik)
python -m venv venv
venv\Scripts\activate        # Windows
# source venv/bin/activate   # Mac/Linux
pip install -r requirements.txt


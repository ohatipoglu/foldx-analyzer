# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FoldX Analyzer is a Python desktop GUI application that automatically analyzes FoldX protein engineering output files (`.fxout` format). It detects which FoldX command generated the file, processes the data, and produces matplotlib visualizations and CSV exports. The UI text is in Turkish.

## Running the Application

```bash
python combined_gui.py      # Launches both FoldX Analyzer and PDB Analyzer windows
python foldx_analysis.py    # Only FoldX Analyzer GUI
python pdb_analyzer.py      # Only PDB/PLIP Analyzer GUI
```

## Setup

```bash
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

There are no automated tests in this project.

## Architecture

### Entry Points
- **`combined_gui.py`** — Shared `Tk()` root that launches two `Toplevel` windows: `FoldXAnalyzerGUI` and `PDBAnalyzerApp`.
- **`foldx_analysis.py`** — Core analysis engine (`FoldXAnalyzerGUI` class, ~528 lines).
- **`pdb_analyzer.py`** — Protein structure analysis (`PDBAnalyzerApp` class, ~462 lines) with Biopython and PLIP tabs.

### FoldX Analysis Data Flow

1. User selects a `.fxout` file or folder of files.
2. `find_header_and_read()` reads the tab-separated file (tries UTF-8 → Latin-1 → CP1254), locating the data table at the line starting with `Pdb\t`.
3. `detect_command()` inspects column names to identify which FoldX command produced the file.
4. The appropriate `_process_*()` method handles command-specific data extraction.
5. A matplotlib figure is generated and embedded in a Tkinter tab via `FigureCanvasTkAgg`.
6. Output is saved to `output/{filename_without_ext}/`: a 300 DPI PNG and a CSV.

### Supported FoldX Commands
| Command | Chart Type | Color Logic |
|---|---|---|
| PositionScan | Bar | Green = stabilizing (ΔΔG < 0), Red = destabilizing |
| RepairPDB | Line | Model progression |
| BuildModel | Line | Order sequence |
| AnalyseComplex | Stacked bar | Energy components |
| Stability | Single bar | Energy value |
| Pssm | Heatmap | Numeric matrix |
| RnaScan | Bar | RNA bases A/U/G/C |

### PDB Analyzer
- **Biopython tab**: Identifies CYS residues and neighbors within a configurable distance threshold (default 6.0 Å); checks disulfide bond feasibility between CYS 166 (chain A) and CYS 56 (chain B).
- **PLIP tab**: Detects protein-ligand interactions (H-bonds, hydrophobic, salt bridges, pi-stacking, pi-cation), generates PyMOL scripts, and renders 300 DPI images. Requires PyMOL installed system-wide.

### Key Dependencies
- `pandas`, `numpy`, `matplotlib` — data processing and visualization
- `tkinter` — GUI (built-in)
- `biopython`, `plip` — structure analysis in `pdb_analyzer.py`
- `PyMOL` — external system dependency for PLIP rendering (optional, gracefully handled)

### Sample Data
- `data/` — 11 sample `.fxout` files covering multiple FoldX commands and PDB structures
- `output/` — Pre-generated analysis results showing expected output structure

"""Shared constants for the FoldX Analyzer application."""

# --- FoldX command names ---
CMD_POSITIONSCAN = "PositionScan"
CMD_REPAIRPDB = "RepairPDB"
CMD_BUILDMODEL = "BuildModel"
CMD_ANALYSECOMPLEX = "AnalyseComplex"
CMD_STABILITY = "Stability"
CMD_PSSM = "Pssm"
CMD_RNASCAN = "RnaScan"
CMD_UNKNOWN = "Bilinmiyor"
ALL_COMMANDS = [
    CMD_POSITIONSCAN, CMD_REPAIRPDB, CMD_BUILDMODEL,
    CMD_ANALYSECOMPLEX, CMD_STABILITY, CMD_PSSM, CMD_RNASCAN,
]

# --- Biological sequences ---
AMINO_ACIDS = ['A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I',
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
RNA_BASES = ['A', 'U', 'G', 'C']

# --- Column name candidates (case-insensitive lookup lists) ---
ENERGY_COLS_INTERACTION = ['Interaction Energy', 'interaction energy']
ENERGY_COLS_TOTAL = ['total energy', 'Total Energy', 'Energy', 'energy']
ENERGY_COLS_AVERAGE = ['average', 'Average', 'total energy', 'Total Energy',
                       'Energy', 'energy']
BACKBONE_COLS = ['backbone', 'Backbone']
SIDECHAIN_COLS = ['sidechain', 'Sidechain', 'side chain', 'Side Chain']
INTERACTION_COLS = ['interaction', 'Interaction']

# --- File I/O ---
FILE_ENCODINGS = ['utf-8', 'latin-1', 'cp1254']
FXOUT_EXTENSION = '.fxout'

# --- UI colours ---
BG_COLOR = '#f4f6f9'
STATUS_BG = '#ecf0f1'
STATUS_FG = '#2c3e50'

# --- Graph rendering ---
FIGURE_SIZE = (10, 6)
FIGURE_DPI = 100
SAVE_DPI = 300
BUILDMODEL_MAX_N = 60
ANALYSECOMPLEX_MAX_N = 40
ENERGY_THRESHOLD_GREEN = -0.5   # ΔΔG below this → green bar
ENERGY_THRESHOLD_RED = 0.5      # ΔΔG above this → red bar

# --- PyMOL / PLIP ---
PYMOL_WIDTH = 1200
PYMOL_HEIGHT = 900
PLIP_IMAGE_PREVIEW_SIZE = (800, 600)

# --- PDB analysis defaults ---
DEFAULT_NEIGHBOR_THRESHOLD = 6.0
DEFAULT_DISULFIDE_MIN = 2.0
DEFAULT_DISULFIDE_MAX = 2.2
DEFAULT_TARGET_RESNUM = 166
DEFAULT_TARGET_CHAIN = "A"

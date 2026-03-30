import configparser
import os
import sys

"""Shared constants for the FoldX Analyzer application.
Loads user-configurable settings from config.ini if available.
"""

# Determine the base path of the application (handles PyInstaller case)
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# --- Load configuration ---
config = configparser.ConfigParser()
config_path = os.path.join(BASE_DIR, 'config.ini')

if not os.path.exists(config_path):
    # Create default config.ini if it doesn't exist
    config['Settings'] = {'pymol_executable_path': 'pymol'}
    config['Graph'] = {'figure_dpi': '100', 'save_dpi': '300',
                       'stabilizing_threshold': '-0.5',
                       'destabilizing_threshold': '0.5'}
    config['PDB'] = {'neighbor_threshold': '6.0',
                     'disulfide_min_dist': '2.0',
                     'disulfide_max_dist': '2.2'}
    with open(config_path, 'w') as f:
        config.write(f)
else:
    config.read(config_path)

# --- Configured Settings ---
PYMOL_EXECUTABLE_PATH = config.get('Settings', 'pymol_executable_path', fallback='pymol')

FIGURE_DPI = config.getint('Graph', 'figure_dpi', fallback=100)
SAVE_DPI = config.getint('Graph', 'save_dpi', fallback=300)
ENERGY_THRESHOLD_GREEN = config.getfloat('Graph', 'stabilizing_threshold', fallback=-0.5)
ENERGY_THRESHOLD_RED = config.getfloat('Graph', 'destabilizing_threshold', fallback=0.5)

DEFAULT_NEIGHBOR_THRESHOLD = config.getfloat('PDB', 'neighbor_threshold', fallback=6.0)
DEFAULT_DISULFIDE_MIN = config.getfloat('PDB', 'disulfide_min_dist', fallback=2.0)
DEFAULT_DISULFIDE_MAX = config.getfloat('PDB', 'disulfide_max_dist', fallback=2.2)


# --- FoldX command names (Not configurable by user) ---
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

# --- Graph rendering ---
FIGURE_SIZE = (10, 6)
BUILDMODEL_MAX_N = 60
ANALYSECOMPLEX_MAX_N = 40

# --- PyMOL / PLIP ---
PYMOL_WIDTH = 1200
PYMOL_HEIGHT = 900
PLIP_IMAGE_PREVIEW_SIZE = (800, 600)

# --- PDB analysis defaults ---
DEFAULT_TARGET_RESNUM = 166
DEFAULT_TARGET_CHAIN = "A"

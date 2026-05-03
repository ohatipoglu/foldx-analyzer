import configparser
import logging
import os
import sys
import warnings

"""Shared constants for the FoldX Analyzer application.
Loads user-configurable settings from config.ini if available.
"""

# Determine the base path of the application (handles PyInstaller case)
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def _validate_config(config: configparser.ConfigParser) -> list[str]:
    """Validate config values and return list of warning messages.
    
    Args:
        config: ConfigParser object with loaded settings
        
    Returns:
        List of warning messages for invalid values
    """
    warnings_list = []
    
    try:
        figure_dpi = config.getint('Graph', 'figure_dpi', fallback=100)
        if not (10 <= figure_dpi <= 600):
            warnings_list.append(
                f"config.ini: figure_dpi={figure_dpi} is invalid (10-600), using default 100."
            )
            config.set('Graph', 'figure_dpi', '100')
    except ValueError:
        warnings_list.append("config.ini: figure_dpi is not a valid integer, using default 100.")
        config.set('Graph', 'figure_dpi', '100')
    
    try:
        save_dpi = config.getint('Graph', 'save_dpi', fallback=300)
        if not (72 <= save_dpi <= 1200):
            warnings_list.append(
                f"config.ini: save_dpi={save_dpi} is invalid (72-1200), using default 300."
            )
            config.set('Graph', 'save_dpi', '300')
    except ValueError:
        warnings_list.append("config.ini: save_dpi is not a valid integer, using default 300.")
        config.set('Graph', 'save_dpi', '300')
    
    try:
        neighbor_threshold = config.getfloat('PDB', 'neighbor_threshold', fallback=6.0)
        if not (0.1 <= neighbor_threshold <= 50.0):
            warnings_list.append(
                f"config.ini: neighbor_threshold={neighbor_threshold} is invalid (0.1-50.0), using default 6.0."
            )
            config.set('PDB', 'neighbor_threshold', '6.0')
    except ValueError:
        warnings_list.append("config.ini: neighbor_threshold is not a valid float, using default 6.0.")
        config.set('PDB', 'neighbor_threshold', '6.0')
    
    try:
        disulfide_min = config.getfloat('PDB', 'disulfide_min_dist', fallback=2.0)
        disulfide_max = config.getfloat('PDB', 'disulfide_max_dist', fallback=2.2)
        if not (1.5 <= disulfide_min < disulfide_max <= 3.0):
            warnings_list.append(
                "config.ini: disulfide distances are invalid (min must be >= 1.5, max <= 3.0, min < max), using defaults (2.0/2.2)."
            )
            config.set('PDB', 'disulfide_min_dist', '2.0')
            config.set('PDB', 'disulfide_max_dist', '2.2')
    except ValueError:
        warnings_list.append("config.ini: disulfide distances are not valid floats, using defaults (2.0/2.2).")
        config.set('PDB', 'disulfide_min_dist', '2.0')
        config.set('PDB', 'disulfide_max_dist', '2.2')
    
    return warnings_list


def _create_default_config() -> configparser.ConfigParser:
    """Create a new ConfigParser with default values.
    
    Returns:
        ConfigParser object with default settings
    """
    config = configparser.ConfigParser()
    config['Settings'] = {'pymol_executable_path': 'pymol'}
    config['Graph'] = {
        'figure_dpi': '100',
        'save_dpi': '300',
        'stabilizing_threshold': '-0.5',
        'destabilizing_threshold': '0.5'
    }
    config['PDB'] = {
        'neighbor_threshold': '6.0',
        'disulfide_min_dist': '2.0',
        'disulfide_max_dist': '2.2'
    }
    return config


# --- Load configuration ---
config = configparser.ConfigParser()
config_path = os.path.join(BASE_DIR, 'config.ini')

if not os.path.exists(config_path):
    # Create default config.ini if it doesn't exist
    config = _create_default_config()
    with open(config_path, 'w', encoding='utf-8') as f:
        config.write(f)
    print(f"Created default config.ini at: {config_path}")
else:
    try:
        config.read(config_path, encoding='utf-8')
    except Exception as e:
        warnings.warn(f"Failed to read config.ini: {e}. Using defaults.")
        config = _create_default_config()

# Validate and fix config if needed
validation_warnings = _validate_config(config)
for warning_msg in validation_warnings:
    warnings.warn(warning_msg)

# If validation made changes, save the corrected config
if validation_warnings:
    with open(config_path, 'w', encoding='utf-8') as f:
        config.write(f)

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
CMD_UNKNOWN = "Unknown"
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

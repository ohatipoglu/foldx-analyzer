"""
Pure business-logic layer for FoldX analysis.

No tkinter or UI dependencies — every function here can be called and tested
independently of any running GUI.
"""

import io
import logging
import os
import re

import numpy as np
import pandas as pd

from constants import (
    FILE_ENCODINGS, FXOUT_EXTENSION,
    AMINO_ACIDS, RNA_BASES,
    ENERGY_COLS_INTERACTION, ENERGY_COLS_TOTAL, ENERGY_COLS_AVERAGE,
    BACKBONE_COLS, SIDECHAIN_COLS, INTERACTION_COLS,
    CMD_POSITIONSCAN, CMD_REPAIRPDB, CMD_BUILDMODEL,
    CMD_ANALYSECOMPLEX, CMD_STABILITY, CMD_PSSM, CMD_RNASCAN, CMD_UNKNOWN,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# File reading
# ---------------------------------------------------------------------------

def find_header_and_read(
    file_path: str,
) -> tuple[pd.DataFrame, str, int, list[str]]:
    """Read an .fxout file and locate its tab-separated data table.

    Tries each encoding in FILE_ENCODINGS in order.  Reads the file once
    per attempt (using StringIO so pandas never re-opens the file).

    Returns:
        df              — parsed DataFrame starting at the 'Pdb' header row
        encoding_used   — encoding that succeeded
        header_row_idx  — 0-based line index of the 'Pdb' header row
        header_lines    — lines that precede the data table (metadata)
    """
    last_error: Exception | None = None
    for enc in FILE_ENCODINGS:
        try:
            with open(file_path, 'r', encoding=enc) as f:
                raw = f.read()
            lines = raw.splitlines(keepends=True)
            for i, line in enumerate(lines):
                if line.strip().startswith('Pdb\t'):
                    buf = io.StringIO(raw)
                    # skiprows skips the first i lines; encoding arg is ignored
                    # for StringIO (data already decoded), so omit it.
                    df = pd.read_csv(buf, sep='\t', skiprows=i)
                    header_lines = [ln.rstrip('\n') for ln in lines[:i]]
                    return df, enc, i, header_lines
        except Exception as e:
            logger.debug("find_header_and_read: %s encoding=%s error=%s",
                         file_path, enc, e)
            last_error = e
    raise ValueError(
        f"Cannot read file: {file_path}. Last error: {last_error}"
    )


# ---------------------------------------------------------------------------
# Command detection
# ---------------------------------------------------------------------------

def detect_command(
    df: pd.DataFrame,
    header_lines: list[str],
    file_basename: str = "",
) -> str:
    """Identify which FoldX command produced the given dataframe.

    Args:
        df            — the parsed data table
        header_lines  — lines before the data table (used for keyword checks)
        file_basename — filename without directory (replaces the old
                        self.file_path.get() dependency)
    """
    cols = set(df.columns.str.lower())
    pdb_col = df['Pdb'] if 'Pdb' in df.columns else pd.Series([], dtype=str)
    pdb_sample = str(pdb_col.iloc[0]).lower() if len(pdb_col) > 0 else ""
    header_text = ' '.join(header_lines).lower()

    if 'interaction energy' in cols and any(
            'wt_' in str(p).lower() for p in pdb_col):
        return CMD_POSITIONSCAN

    if 'total energy' in cols and re.search(r'repair_\d+', pdb_sample):
        return CMD_REPAIRPDB

    if 'average' in cols or 'raw' in file_basename.lower():
        return CMD_BUILDMODEL

    if all(x in cols for x in ['backbone', 'sidechain', 'interaction']):
        return CMD_ANALYSECOMPLEX

    if 'total energy' in cols and len(df) == 1:
        return CMD_STABILITY

    if 'pssm' in header_text:
        return CMD_PSSM

    if (any(base in pdb_sample for base in ['a', 'u', 'g', 'c'])
            and 'rna' in header_text):
        return CMD_RNASCAN

    return CMD_UNKNOWN


# ---------------------------------------------------------------------------
# Column / value helpers
# ---------------------------------------------------------------------------

def first_existing_column(
    df: pd.DataFrame,
    candidates: list[str],
) -> str | None:
    """Return the first column from *candidates* found in *df* (case-insensitive)."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


# ---------------------------------------------------------------------------
# Shared scan-row parsing (PositionScan / RnaScan)
# ---------------------------------------------------------------------------

def _parse_scan_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate a scan DataFrame with 'run', 'replicate', and 'is_wt' columns.

    Uses list accumulation rather than df.at[i, col] to avoid relying on
    a contiguous 0-based integer index.
    """
    df = df.copy().reset_index(drop=True)
    is_wt_flags = df['Pdb'].astype(str).str.contains('WT_', na=False).tolist()
    run_list: list[int | None] = []
    rep_list: list[int | None] = []
    prev_run: int | None = None
    prev_rep: int | None = None

    for pdb, is_wt in zip(df['Pdb'].astype(str), is_wt_flags):
        if is_wt:
            run_list.append(prev_run)
            rep_list.append(prev_rep)
        else:
            try:
                parts = pdb.split('_')
                run = int(parts[-2])
                rep = int(parts[-1].split('.')[0])
                prev_run, prev_rep = run, rep
                run_list.append(run)
                rep_list.append(rep)
            except Exception:
                logger.debug("_parse_scan_rows: cannot parse Pdb name '%s'", pdb)
                run_list.append(None)
                rep_list.append(None)

    df['is_wt'] = is_wt_flags
    df['run'] = run_list
    df['replicate'] = rep_list
    return df


def _compute_scan_energies(
    df: pd.DataFrame,
    labels: list[str],
    label_key: str,
) -> list[dict]:
    """Compute ΔΔG for each label by comparing mutant to WT rows.

    Args:
        df        — output of _parse_scan_rows
        labels    — ordered symbol list (AMINO_ACIDS or RNA_BASES)
        label_key — dict key for the symbol ('amino_acid' or 'base')
    """
    energy_col = first_existing_column(df, ENERGY_COLS_INTERACTION)
    if not energy_col:
        logger.warning("_compute_scan_energies: 'Interaction Energy' column not found.")
        return []

    results: list[dict] = []
    for run in df['run'].dropna().unique():
        for rep in df['replicate'].dropna().unique():
            mask = (df['run'] == run) & (df['replicate'] == rep)
            mut = df[mask & ~df['is_wt']]
            wt = df[mask & df['is_wt']]
            if mut.empty or wt.empty:
                continue
            try:
                run_int = int(run)
                if 1 <= run_int <= len(labels):
                    mut_e = pd.to_numeric(mut[energy_col].iloc[0], errors='coerce')
                    wt_e = pd.to_numeric(wt[energy_col].iloc[0], errors='coerce')
                    results.append({
                        label_key: labels[run_int - 1],
                        'energy': float(mut_e - wt_e),
                    })
                else:
                    logger.debug("_compute_scan_energies: invalid run index %d", run_int)
            except Exception as e:
                logger.debug("_compute_scan_energies: run=%s error=%s", run, e)
    return results


# ---------------------------------------------------------------------------
# Per-command processors
# ---------------------------------------------------------------------------

def process_positionscan(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    """Process PositionScan data, aggregating energies by amino acid."""
    df_parsed = _parse_scan_rows(df)
    results_list = _compute_scan_energies(df_parsed, AMINO_ACIDS, 'amino_acid')
    if not results_list:
        return pd.DataFrame(), AMINO_ACIDS

    stats = (pd.DataFrame(results_list)
             .groupby('amino_acid')['energy']
             .agg(['mean', 'std', 'count'])
             .reset_index())
    stats['sem'] = stats['std'] / np.sqrt(stats['count'])
    return stats, AMINO_ACIDS


def process_rnascan(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    """Process RnaScan data, aggregating energies by base."""
    df_parsed = _parse_scan_rows(df)
    results_list = _compute_scan_energies(df_parsed, RNA_BASES, 'base')
    if not results_list:
        return pd.DataFrame(), RNA_BASES

    stats = (pd.DataFrame(results_list).groupby('base')['energy']
             .agg(['mean', 'std', 'count']).reset_index())
    stats['sem'] = stats['std'] / np.sqrt(stats['count'])
    return stats, RNA_BASES


def process_repairpdb(df: pd.DataFrame) -> tuple[pd.DataFrame, None]:
    df = df.copy()
    df['_model'] = df['Pdb'].astype(str).str.extract(r'_(\d+)').astype(float)
    energy_col = first_existing_column(df, ENERGY_COLS_TOTAL)
    if not energy_col:
        logger.warning("process_repairpdb: 'total energy' column not found.")
        return [], None
    results: list[dict] = []
    for _, row in df.iterrows():
        if pd.isna(row['_model']):
            logger.debug("process_repairpdb: no model number for Pdb=%s", row.get('Pdb'))
            continue
        val = pd.to_numeric(row[energy_col], errors='coerce')
        if not pd.isna(val):
            results.append({'model': int(row['_model']), 'energy': float(val)})
    return pd.DataFrame(results), None


def process_buildmodel(df: pd.DataFrame) -> tuple[pd.DataFrame, None]:
    energy_col = first_existing_column(df, ENERGY_COLS_AVERAGE)
    if not energy_col:
        logger.warning("process_buildmodel: no energy column found.")
        return [], None
    pdb_series = df['Pdb'] if 'Pdb' in df.columns else pd.Series(range(len(df)), dtype=str)
    energies = pd.to_numeric(df[energy_col], errors='coerce')
    results = [
        {'pdb': str(pdb), 'energy': float(val)}
        for pdb, val in zip(pdb_series, energies)
        if not pd.isna(val)
    ]
    return pd.DataFrame(results), None


def process_analysecomplex(df: pd.DataFrame) -> tuple[pd.DataFrame, None]:
    back_col = first_existing_column(df, BACKBONE_COLS)
    side_col = first_existing_column(df, SIDECHAIN_COLS)
    inter_col = first_existing_column(df, INTERACTION_COLS)
    if not any([back_col, side_col, inter_col]):
        logger.warning("process_analysecomplex: no backbone/sidechain/interaction columns.")
        return [], None
    results: list[dict] = []
    for _, row in df.iterrows():
        rec: dict = {'pdb': str(row.get('Pdb', ''))}
        for key, col in [('backbone', back_col),
                         ('sidechain', side_col),
                         ('interaction', inter_col)]:
            if col:
                rec[key] = pd.to_numeric(row[col], errors='coerce')
        results.append(rec)
    df_res = pd.DataFrame(results)
    comp_cols = [c for c in ['backbone', 'sidechain', 'interaction']
                 if c in df_res.columns]
    return df_res.dropna(subset=comp_cols, how='all'), None


def process_stability(df: pd.DataFrame) -> tuple[pd.DataFrame, None]:
    energy_col = first_existing_column(df, ENERGY_COLS_TOTAL)
    if not energy_col or len(df) == 0:
        return pd.DataFrame(), None
    val = pd.to_numeric(df[energy_col].iloc[0], errors='coerce')
    if pd.isna(val):
        return pd.DataFrame(), None
    pdb_val = str(df['Pdb'].iloc[0]) if 'Pdb' in df.columns else ''
    return pd.DataFrame([{'pdb': pdb_val, 'energy': float(val)}]), None


def process_pssm(df: pd.DataFrame) -> tuple[pd.DataFrame, None]:
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    numeric_df = df.apply(pd.to_numeric, errors='coerce')
    if numeric_df.notna().sum().sum() == 0:
        return pd.DataFrame(), None
    # Return the original df, as the renderer expects string data for columns/rows
    # and will perform numeric conversion itself.
    return df, None


# ---------------------------------------------------------------------------
# Dispatch table
# ---------------------------------------------------------------------------

PROCESSORS: dict[str, callable] = {
    CMD_POSITIONSCAN:   process_positionscan,
    CMD_REPAIRPDB:      process_repairpdb,
    CMD_BUILDMODEL:     process_buildmodel,
    CMD_ANALYSECOMPLEX: process_analysecomplex,
    CMD_STABILITY:      process_stability,
    CMD_PSSM:           process_pssm,
    CMD_RNASCAN:        process_rnascan,
}


# ---------------------------------------------------------------------------
# File collection
# ---------------------------------------------------------------------------

def collect_files(input_path: str) -> list[str]:
    """Return .fxout file paths from a file or directory path."""
    if os.path.isfile(input_path):
        return [input_path]
    if os.path.isdir(input_path):
        return [
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if f.endswith(FXOUT_EXTENSION)
        ]
    return []

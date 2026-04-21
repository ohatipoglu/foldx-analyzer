import os
import tempfile
import unittest

import pandas as pd

from foldx_core import (
    detect_command,
    find_header_and_read,
    first_existing_column,
    process_analysecomplex,
    process_buildmodel,
    process_positionscan,
    process_repairpdb,
    process_rnascan,
    process_stability,
    CMD_POSITIONSCAN,
    CMD_STABILITY,
    CMD_REPAIRPDB,
    CMD_BUILDMODEL,
    CMD_ANALYSECOMPLEX,
)


class TestDetectCommand(unittest.TestCase):

    def _write_fxout(self, content: str) -> str:
        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.fxout', delete=False, encoding='utf-8'
        )
        tmp.write(content)
        tmp.close()
        return tmp.name

    def test_determine_command_positionscan(self):
        path = self._write_fxout("Pdb\tinteraction energy\nWT_1a22\t-1.5\n")
        try:
            df, _, _, hdr = find_header_and_read(path)
            self.assertEqual(detect_command(df, hdr, "PS_protein.fxout"), CMD_POSITIONSCAN)
        finally:
            os.remove(path)

    def test_determine_command_stability(self):
        path = self._write_fxout("Pdb\ttotal energy\n1a22\t-1.5\n")
        try:
            df, _, _, hdr = find_header_and_read(path)
            self.assertEqual(detect_command(df, hdr, "Stability_of_protein.fxout"), CMD_STABILITY)
        finally:
            os.remove(path)

    def test_determine_command_repairpdb(self):
        path = self._write_fxout("Pdb\ttotal energy\nrepair_0\t-100.0\nrepair_1\t-99.0\n")
        try:
            df, _, _, hdr = find_header_and_read(path)
            self.assertEqual(detect_command(df, hdr, "RepairPDB.fxout"), CMD_REPAIRPDB)
        finally:
            os.remove(path)

    def test_determine_command_buildmodel(self):
        path = self._write_fxout("Pdb\taverage\nmodel_0\t-5.0\nmodel_1\t-4.0\n")
        try:
            df, _, _, hdr = find_header_and_read(path)
            self.assertEqual(detect_command(df, hdr, "raw_BuildModel.fxout"), CMD_BUILDMODEL)
        finally:
            os.remove(path)

    def test_determine_command_analysecomplex(self):
        path = self._write_fxout("Pdb\tbackbone\tsidechain\tinteraction\ncpx_1\t1.0\t2.0\t3.0\n")
        try:
            df, _, _, hdr = find_header_and_read(path)
            self.assertEqual(detect_command(df, hdr, "AC_result.fxout"), CMD_ANALYSECOMPLEX)
        finally:
            os.remove(path)


class TestProcessors(unittest.TestCase):

    def test_process_positionscan_basic(self):
        df = pd.DataFrame({
            'Pdb': ['protein_1_1.pdb', 'WT_1_1.pdb', 'protein_2_1.pdb', 'WT_2_1.pdb'],
            'interaction energy': [-5.0, -4.0, -3.0, -2.0],
        })
        result, labels = process_positionscan(df)
        self.assertFalse(result.empty)
        self.assertIn('amino_acid', result.columns)
        self.assertIn('mean', result.columns)
        self.assertIn('sem', result.columns)
        # run=1 → AMINO_ACIDS[0]='A', ΔΔG = -5 - (-4) = -1
        row_a = result[result['amino_acid'] == 'A']
        self.assertAlmostEqual(float(row_a['mean'].iloc[0]), -1.0)

    def test_process_positionscan_empty(self):
        df = pd.DataFrame({'Pdb': [], 'interaction energy': []})
        result, labels = process_positionscan(df)
        self.assertTrue(result.empty)

    def test_process_repairpdb_basic(self):
        df = pd.DataFrame({
            'Pdb': ['repair_0', 'repair_1', 'repair_2'],
            'total energy': [-100.0, -95.0, -110.0],
        })
        result, _ = process_repairpdb(df)
        self.assertEqual(len(result), 3)
        self.assertIn('model', result.columns)
        self.assertIn('energy', result.columns)
        self.assertListEqual(sorted(result['model'].tolist()), [0, 1, 2])

    def test_process_repairpdb_no_energy_column(self):
        df = pd.DataFrame({'Pdb': ['repair_0'], 'other_col': [1.0]})
        result, _ = process_repairpdb(df)
        # Should return empty list/df, not raise
        self.assertEqual(len(result), 0)

    def test_process_buildmodel_basic(self):
        df = pd.DataFrame({
            'Pdb': ['model_0', 'model_1'],
            'average': [-5.0, -3.0],
        })
        result, _ = process_buildmodel(df)
        self.assertEqual(len(result), 2)
        self.assertIn('energy', result.columns)

    def test_process_buildmodel_skips_nan(self):
        df = pd.DataFrame({
            'Pdb': ['model_0', 'model_1', 'model_2'],
            'average': [-5.0, None, -3.0],
        })
        result, _ = process_buildmodel(df)
        self.assertEqual(len(result), 2)

    def test_process_stability_basic(self):
        df = pd.DataFrame({'Pdb': ['structure_1'], 'total energy': [-10.5]})
        result, _ = process_stability(df)
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(float(result['energy'].iloc[0]), -10.5)

    def test_process_stability_empty(self):
        df = pd.DataFrame({'Pdb': [], 'total energy': []})
        result, _ = process_stability(df)
        self.assertTrue(result.empty)

    def test_process_analysecomplex_basic(self):
        df = pd.DataFrame({
            'Pdb': ['complex_1', 'complex_2'],
            'backbone': [1.0, 2.0],
            'sidechain': [3.0, 4.0],
            'interaction': [5.0, 6.0],
        })
        result, _ = process_analysecomplex(df)
        self.assertEqual(len(result), 2)
        self.assertIn('backbone', result.columns)

    def test_process_rnascan_basic(self):
        df = pd.DataFrame({
            'Pdb': ['rna_1_1.pdb', 'WT_1_1.pdb'],
            'interaction energy': [-3.0, -2.0],
        })
        result, bases = process_rnascan(df)
        self.assertFalse(result.empty)
        self.assertIn('base', result.columns)


class TestHelpers(unittest.TestCase):

    def test_first_existing_column_found(self):
        df = pd.DataFrame({'Total Energy': [1.0], 'Pdb': ['x']})
        col = first_existing_column(df, ['total energy', 'energy'])
        self.assertEqual(col, 'Total Energy')

    def test_first_existing_column_missing(self):
        df = pd.DataFrame({'Pdb': ['x']})
        col = first_existing_column(df, ['total energy'])
        self.assertIsNone(col)

    def test_find_header_and_read_utf8(self):
        content = "Pdb\ttotal energy\nprotein_1\t-5.0\n"
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.fxout', delete=False, encoding='utf-8'
        ) as f:
            f.write(content)
            path = f.name
        try:
            df, encoding_used, _, _ = find_header_and_read(path)
            self.assertEqual(encoding_used, 'utf-8')
            self.assertIn('total energy', df.columns)
        finally:
            os.remove(path)

    def test_find_header_and_read_latin1_fallback(self):
        # 'é' (0xe9) is valid Latin-1 but not valid UTF-8 — forces fallback
        content_bytes = b"Pdb\ttotal energy\nprote\xedn_1\t-5.0\n"
        with tempfile.NamedTemporaryFile(suffix='.fxout', delete=False) as f:
            f.write(content_bytes)
            path = f.name
        try:
            df, encoding_used, _, _ = find_header_and_read(path)
            self.assertIn(encoding_used, ('latin-1', 'cp1254'))
            self.assertIn('total energy', df.columns)
        finally:
            os.remove(path)

    def test_find_header_and_read_invalid_file(self):
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.fxout', delete=False, encoding='utf-8'
        ) as f:
            f.write("no header here\njust some text\n")
            path = f.name
        try:
            with self.assertRaises(ValueError):
                find_header_and_read(path)
        finally:
            os.remove(path)


if __name__ == '__main__':
    unittest.main()

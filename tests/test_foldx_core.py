import os
import tempfile
import unittest

import pandas as pd

from foldx_core import detect_command, CMD_POSITIONSCAN, CMD_STABILITY, find_header_and_read


class TestFoldXAnalyzer(unittest.TestCase):

    def test_determine_command_positionscan(self):
        """Test if PositionScan command is correctly identified from filename/data."""
        # Need to create mock data because detect_command takes a dataframe
        mock_data = "Pdb\tinteraction energy\nWT_1a22\t-1.5\n"
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fxout', delete=False) as tmp:
            tmp.write(mock_data)
            tmp_path = tmp.name

        try:
            df, _, _, header_lines = find_header_and_read(tmp_path)
            command = detect_command(df, header_lines, "PS_protein_mutations.fxout")
            self.assertEqual(command, CMD_POSITIONSCAN)
        finally:
            os.remove(tmp_path)

    def test_determine_command_stability(self):
        """Test if Stability command is correctly identified."""
        mock_data = "Pdb\ttotal energy\n1a22\t-1.5\n"
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fxout', delete=False) as tmp:
            tmp.write(mock_data)
            tmp_path = tmp.name

        try:
            df, _, _, header_lines = find_header_and_read(tmp_path)
            command = detect_command(df, header_lines, "Stability_of_protein.fxout")
            self.assertEqual(command, CMD_STABILITY)
        finally:
            os.remove(tmp_path)


if __name__ == '__main__':
    unittest.main()

"""Unit tests for pdb_analyzer module.

Tests cover PLIP TXT parsing, mutation name extraction,
and basic structural analysis functions.
"""

import os
import tempfile
import unittest
from typing import Any, Dict, List

import numpy as np

from pdb_analyzer import (
    parse_plip_txt_file,
    get_mutation_name,
    create_plip_comparison_word_report,
)

# Bio.PDB tests require biopython installed
try:
    from Bio.PDB import PDBParser
    from pdb_analyzer import (
        load_structure,
        list_all_cys,
        find_neighbors,
        calculate_disulfide_distance,
    )
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


class TestParsePlipTxtFile(unittest.TestCase):
    """Tests for parse_plip_txt_file function."""

    def _create_plip_txt(self, content: str) -> str:
        """Create a temporary PLIP TXT file with given content."""
        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.txt', delete=False, encoding='utf-8'
        )
        tmp.write(content)
        tmp.close()
        return tmp.name

    def test_parse_empty_file(self) -> None:
        """Test parsing an empty PLIP TXT file."""
        path = self._create_plip_txt("")
        try:
            result = parse_plip_txt_file(path)
            self.assertIsInstance(result, dict)
            self.assertEqual(len(result['Hydrophobic']), 0)
            self.assertEqual(len(result['Hydrogen Bond']), 0)
        finally:
            os.remove(path)

    def test_parse_file_not_found(self) -> None:
        """Test parsing a non-existent file."""
        result = parse_plip_txt_file("/nonexistent/path/file.txt")
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result['Hydrophobic']), 0)

    def test_parse_with_hydrophobic_interaction(self) -> None:
        """Test parsing PLIP TXT with hydrophobic interactions."""
        content = """
## Details to all ligand-protein interactions

================================================================================
============================== Ligand: UNK ====================================
================================================================================

**Hydrophobic Interactions**

| RESNR | Restype | Chain | Resnr | Ligtype | Ligchain | Distance |
|-------|---------|-------|-------|---------|----------|----------|
| 50    | LEU     | A     | 100   | UNK     | B        | 3.8      |
| 53    | VAL     | A     | 100   | UNK     | B        | 3.9      |

**Hydrogen Bonds**

| No interactions |
"""
        path = self._create_plip_txt(content)
        try:
            result = parse_plip_txt_file(path, target_residues={50, 53})
            self.assertEqual(len(result['Hydrophobic']), 2)
            self.assertEqual(result['Hydrophobic'][0]['prot_res'], 50)
            self.assertEqual(result['Hydrophobic'][0]['prot_type'], 'LEU')
        finally:
            os.remove(path)

    def test_parse_with_hydrogen_bonds(self) -> None:
        """Test parsing PLIP TXT with hydrogen bonds."""
        content = """
**Hydrogen Bonds**

| RESNR | Restype | Chain | Resnr | Ligtype | Ligchain | Distance | D-H..A |
|-------|---------|-------|-------|---------|----------|----------|--------|
| 63    | SER     | A     | 100   | UNK     | B        | 2.9      | 145    |

**Hydrophobic Interactions**

| No interactions |
"""
        path = self._create_plip_txt(content)
        try:
            result = parse_plip_txt_file(path)
            self.assertEqual(len(result['Hydrogen Bond']), 1)
            self.assertEqual(result['Hydrogen Bond'][0]['prot_res'], 63)
        finally:
            os.remove(path)

    def test_parse_with_target_filter(self) -> None:
        """Test parsing with target residue filtering."""
        content = """
**Hydrophobic Interactions**

| RESNR | Restype | Chain | Resnr | Ligtype | Ligchain | Distance |
|-------|---------|-------|-------|---------|----------|----------|
| 50    | LEU     | A     | 100   | UNK     | B        | 3.8      |
| 90    | PHE     | A     | 100   | UNK     | B        | 3.7      |
"""
        path = self._create_plip_txt(content)
        try:
            result = parse_plip_txt_file(path, target_residues={50})
            self.assertEqual(len(result['Hydrophobic']), 1)
            self.assertEqual(result['Hydrophobic'][0]['prot_res'], 50)
        finally:
            os.remove(path)

    def test_parse_with_chain_filter(self) -> None:
        """Test parsing with chain filtering."""
        content = """
**Hydrophobic Interactions**

| RESNR | Restype | Chain | Resnr | Ligtype | Ligchain | Distance |
|-------|---------|-------|-------|---------|----------|----------|
| 50    | LEU     | A     | 100   | UNK     | B        | 3.8      |
| 50    | LEU     | B     | 100   | UNK     | B        | 3.8      |
"""
        path = self._create_plip_txt(content)
        try:
            result = parse_plip_txt_file(path, target_chain='A')
            self.assertEqual(len(result['Hydrophobic']), 1)
            self.assertEqual(result['Hydrophobic'][0]['prot_chain'], 'A')
        finally:
            os.remove(path)


class TestGetMutationName(unittest.TestCase):
    """Tests for get_mutation_name function."""

    def test_standard_mutation_name(self) -> None:
        """Test extracting mutation name from standard filename."""
        filename = "target_pa_A166C_model.000.00.txt"
        result = get_mutation_name(filename)
        self.assertEqual(result, "A166C")

    def test_mutation_name_without_extension(self) -> None:
        """Test extracting mutation name without .txt extension."""
        filename = "target_pa_G50V_model.000.00"
        result = get_mutation_name(filename)
        self.assertEqual(result, "G50V")

    def test_simple_filename(self) -> None:
        """Test simple filename without mutation pattern."""
        filename = "wildtype.txt"
        result = get_mutation_name(filename)
        self.assertEqual(result, "wildtype")


class TestCreatePlipComparisonWordReport(unittest.TestCase):
    """Tests for create_plip_comparison_word_report function."""

    def test_create_report_basic(self) -> None:
        """Test creating a basic comparison report."""
        all_results: Dict[str, Dict[str, List[Any]]] = {
            'A166C': {'Hydrophobic': [1, 2], 'Hydrogen Bond': [3]},
            'G50V': {'Hydrophobic': [1], 'Hydrogen Bond': [3, 4]},
        }
        with tempfile.NamedTemporaryFile(suffix='.docx', delete=False) as f:
            output_path = f.name
        try:
            result = create_plip_comparison_word_report(
                all_results, output_path, base_residues=[50, 63]
            )
            self.assertTrue(result)
            self.assertTrue(os.path.exists(output_path))
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)

    def test_create_report_empty_results(self) -> None:
        """Test creating report with empty results."""
        all_results: Dict[str, Dict[str, List[Any]]] = {}
        with tempfile.NamedTemporaryFile(suffix='.docx', delete=False) as f:
            output_path = f.name
        try:
            result = create_plip_comparison_word_report(all_results, output_path)
            self.assertTrue(result)
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)


@unittest.skipIf(not BIOPYTHON_AVAILABLE, "Biopython not installed")
class TestLoadStructure(unittest.TestCase):
    """Tests for load_structure function."""

    def _create_minimal_pdb(self) -> str:
        """Create a minimal PDB file for testing."""
        content = """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM      3  C   ALA A   1      12.000  12.000  12.000  1.00  0.00           C
ATOM      4  O   ALA A   1      13.000  13.000  13.000  1.00  0.00           O
TER
END
"""
        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.pdb', delete=False, encoding='utf-8'
        )
        tmp.write(content)
        tmp.close()
        return tmp.name

    def test_load_structure_basic(self) -> None:
        """Test loading a basic PDB structure."""
        path = self._create_minimal_pdb()
        try:
            structure = load_structure(path)
            self.assertIsNotNone(structure)
            self.assertEqual(len(list(structure.get_chains())), 1)
        finally:
            os.remove(path)

    def test_load_structure_invalid_file(self) -> None:
        """Test loading an invalid PDB file."""
        with self.assertRaises(Exception):
            load_structure("/nonexistent/file.pdb")


@unittest.skipIf(not BIOPYTHON_AVAILABLE, "Biopython not installed")
class TestListAllCys(unittest.TestCase):
    """Tests for list_all_cys function."""

    def _create_pdb_with_cys(self) -> str:
        """Create a PDB file with cysteine residues."""
        content = """HEADER    TEST CYS PROTEIN
ATOM      1  N   CYS A  10      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  CYS A  10      11.000  11.000  11.000  1.00  0.00           C
ATOM      3  SG  CYS A  10      12.000  12.000  12.000  1.00  0.00           S
ATOM      4  N   ALA A  20      20.000  20.000  20.000  1.00  0.00           N
ATOM      5  N   CYS B  30      30.000  30.000  30.000  1.00  0.00           N
ATOM      6  SG  CYS B  30      31.000  31.000  31.000  1.00  0.00           S
TER
END
"""
        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.pdb', delete=False, encoding='utf-8'
        )
        tmp.write(content)
        tmp.close()
        return tmp.name

    def test_list_all_cys_basic(self) -> None:
        """Test listing all cysteine residues."""
        path = self._create_pdb_with_cys()
        try:
            structure = load_structure(path)
            cys_list = list_all_cys(structure)
            self.assertEqual(len(cys_list), 2)
            # Check first CYS (chain A, residue 10)
            self.assertEqual(cys_list[0][1], 'A')
            self.assertEqual(cys_list[0][2], 10)
            # Check second CYS (chain B, residue 30)
            self.assertEqual(cys_list[1][1], 'B')
            self.assertEqual(cys_list[1][2], 30)
        finally:
            os.remove(path)

    def test_list_all_cys_no_cys(self) -> None:
        """Test listing cysteines when none exist."""
        content = """HEADER    TEST NO CYS
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
TER
END
"""
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.pdb', delete=False, encoding='utf-8'
        ) as f:
            f.write(content)
            path = f.name
        try:
            structure = load_structure(path)
            cys_list = list_all_cys(structure)
            self.assertEqual(len(cys_list), 0)
        finally:
            os.remove(path)


@unittest.skipIf(not BIOPYTHON_AVAILABLE, "Biopython not installed")
class TestFindNeighbors(unittest.TestCase):
    """Tests for find_neighbors function."""

    def _create_pdb_for_neighbors(self) -> str:
        """Create a PDB file for neighbor testing."""
        content = """HEADER    TEST NEIGHBORS
ATOM      1  N   CYS A 166      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  CYS A 166      11.000  11.000  11.000  1.00  0.00           C
ATOM      3  SG  CYS A 166      12.000  12.000  12.000  1.00  0.00           S
ATOM      4  N   LEU B  50      15.000  15.000  15.000  1.00  0.00           N
ATOM      5  CA  LEU B  50      16.000  16.000  16.000  1.00  0.00           C
ATOM      6  N   VAL B  60      50.000  50.000  50.000  1.00  0.00           N
TER
END
"""
        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.pdb', delete=False, encoding='utf-8'
        )
        tmp.write(content)
        tmp.close()
        return tmp.name

    def test_find_neighbors_basic(self) -> None:
        """Test finding neighbors around target residue."""
        path = self._create_pdb_for_neighbors()
        try:
            structure = load_structure(path)
            # Use large threshold to catch distant residues
            neighbors = find_neighbors(structure, threshold=10.0, target_resnum=166, target_chain='A')
            # Should find at least one neighbor in chain B
            self.assertGreater(len(neighbors), 0)
        finally:
            os.remove(path)

    def test_find_neighbors_no_target(self) -> None:
        """Test finding neighbors when target doesn't exist."""
        path = self._create_pdb_for_neighbors()
        try:
            structure = load_structure(path)
            neighbors = find_neighbors(structure, threshold=5.0, target_resnum=999, target_chain='A')
            self.assertEqual(len(neighbors), 0)
        finally:
            os.remove(path)


@unittest.skipIf(not BIOPYTHON_AVAILABLE, "Biopython not installed")
class TestCalculateDisulfideDistance(unittest.TestCase):
    """Tests for calculate_disulfide_distance function."""

    def _create_pdb_with_disulfide(self) -> str:
        """Create a PDB file with two cysteines for disulfide bridge."""
        # Two CYS residues with SG atoms at known distance
        # Distance between (12, 12, 12) and (22, 22, 22) = sqrt(3*100) = 17.32
        content = """HEADER    TEST DISULFIDE
ATOM      1  N   CYS A 166      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  SG  CYS A 166      12.000  12.000  12.000  1.00  0.00           S
ATOM      3  N   CYS B  56      20.000  20.000  20.000  1.00  0.00           N
ATOM      4  SG  CYS B  56      22.000  22.000  22.000  1.00  0.00           S
TER
END
"""
        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.pdb', delete=False, encoding='utf-8'
        )
        tmp.write(content)
        tmp.close()
        return tmp.name

    def test_calculate_disulfide_distance_basic(self) -> None:
        """Test calculating disulfide bridge distance."""
        path = self._create_pdb_with_disulfide()
        try:
            structure = load_structure(path)
            distance = calculate_disulfide_distance(
                structure, c1_ch='A', c1_r=166, c2_ch='B', c2_r=56
            )
            self.assertIsNotNone(distance)
            if distance is not None:
                expected = np.sqrt(3 * 100)  # ~17.32
                self.assertAlmostEqual(distance, expected, places=2)
        finally:
            os.remove(path)

    def test_calculate_disulfide_distance_missing_cys(self) -> None:
        """Test distance calculation when cysteine is missing."""
        path = self._create_pdb_with_disulfide()
        try:
            structure = load_structure(path)
            distance = calculate_disulfide_distance(
                structure, c1_ch='A', c1_r=999, c2_ch='B', c2_r=56
            )
            self.assertIsNone(distance)
        finally:
            os.remove(path)


if __name__ == '__main__':
    unittest.main()

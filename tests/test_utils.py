"""
tests/test_utils.py
Unit tests for orca_result_analyzer_rust/utils.py (pure Python, no stubs required).
"""

import os
import sys
import importlib.util
import unittest

_SRC = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "orca_result_analyzer_rust", "utils.py"))


def _load_utils():
    spec = importlib.util.spec_from_file_location("orca_utils_mod", _SRC)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["orca_utils_mod"] = mod
    spec.loader.exec_module(mod)
    return mod


_utils = _load_utils()
get_default_export_path = _utils.get_default_export_path


class TestGetDefaultExportPath(unittest.TestCase):

    def test_csv_extension(self):
        result = get_default_export_path("/some/dir/job.out", suffix="_scf_trace", extension=".csv")
        self.assertEqual(result, os.path.join("/some/dir", "job_scf_trace.csv"))

    def test_default_suffix(self):
        result = get_default_export_path("/data/mol.out", extension=".csv")
        self.assertEqual(result, os.path.join("/data", "mol_analyzed.csv"))

    def test_empty_base_path_returns_empty_string(self):
        self.assertEqual(get_default_export_path(""), "")
        self.assertEqual(get_default_export_path(None), "")

    def test_no_directory(self):
        result = get_default_export_path("job.out", suffix="_result", extension=".txt")
        self.assertEqual(result, "job_result.txt")

    def test_no_extension_arg(self):
        result = get_default_export_path("/dir/calc.out", suffix="_data")
        self.assertEqual(result, os.path.join("/dir", "calc_data"))

    def test_preserves_directory(self):
        base = os.path.join("a", "b", "c", "mol.out")
        result = get_default_export_path(base, suffix="_x", extension=".png")
        expected_dir = os.path.join("a", "b", "c")
        self.assertIn(expected_dir, result)
        self.assertTrue(result.endswith("mol_x.png"))


if __name__ == "__main__":
    unittest.main()

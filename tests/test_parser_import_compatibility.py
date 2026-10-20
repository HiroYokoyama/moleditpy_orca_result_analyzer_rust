"""Exercise parser import in a fresh interpreter without the test Qt stubs."""

from pathlib import Path
import subprocess
import sys


def test_structure_parser_import_and_xyz_rows():
    """The geometry parser must import and parse coordinates on Python 3.9."""
    parser_path = (
        Path(__file__).resolve().parents[1]
        / "orca_result_analyzer_rust"
        / "parser_structure.py"
    )
    script = """
import runpy
import sys

parser = runpy.run_path(sys.argv[1])["_StructureParsingMixin"]
assert parser._parse_xyz_row("C 1.0 2.0 3.0") == ("C", [1.0, 2.0, 3.0])
assert parser._parse_xyz_row("invalid row") is None
assert parser._parse_xyz_row("C invalid 2.0 3.0") is None
"""
    result = subprocess.run(
        [sys.executable, "-c", script, str(parser_path)],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr

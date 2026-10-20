"""
tests/test_parser_charges.py
Parser coverage for two charge blocks that the other parser tests do not
reach: the NBO natural-population fallback (used when the summary table is
absent or unparseable) and the frontier-molecular-orbital population table.

Loads OrcaParser directly — no Qt stubs required.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from _parser_loader import load_standalone_parser  # noqa: E402


def _load_parser():
    return load_standalone_parser("orca_parser_charges")


OrcaParser = _load_parser().OrcaParser


def _parse(content):
    p = OrcaParser()
    p.raw_content = content
    p.lines = content.splitlines()
    p.filename = "test.out"
    p.parse_charges()
    return p.data.get("charges", {})


# ---------------------------------------------------------------------------
# NBO natural populations
# ---------------------------------------------------------------------------

NBO_SUMMARY = """
 NATURAL POPULATIONS:  Natural atomic orbital occupancies

  NAO Atom No lang   Type(AO)    Occupancy      Energy
 ----------------------------------------------------
   1    O  1  s      Cor( 1s)     1.99999      -18.9
   2    H  2  s      Val( 1s)     0.60000       -0.5

 Summary of Natural Population Analysis:

                                       Natural Population
              Natural  -----------------------------------------------
   Atom  No    Charge         Core      Valence    Rydberg      Total
 -----------------------------------------------------------------------
      O    1   -0.92345      1.99989     6.91234    0.01122     8.92345
      H    2    0.46172      0.00000     0.53718    0.00110     0.53828
 =======================================================================
"""

# No "Summary of Natural Population Analysis" table -> fallback path.
NBO_FALLBACK = """
 NATURAL POPULATIONS:  Natural atomic orbital occupancies
 ----------------------------------------------------
      O    1   -0.92345
      H    2    0.46172
      H    3    0.46173

 Natural Electron Configuration
"""

# ORCA labels each row with a single "index-symbol" token, e.g. "0-C".
FMO_BLOCK = """
 FRONTIER MOLECULAR ORBITAL POPULATION ANALYSIS

   Atom      HOMO(Mull)   HOMO(Loew)   LUMO(Mull)   LUMO(Loew)
 ---------------------------------------------------------------
   0-O       0.812345     0.798765     0.101234     0.098765
   1-H       0.093827     0.100617     0.449383     0.450617
 ---------------------------------------------------------------
"""

# Some outputs label rows with a bare symbol; indices are then positional.
FMO_BARE_LABELS = """
 FRONTIER MOLECULAR ORBITAL POPULATION ANALYSIS

   Atom      HOMO(Mull)   HOMO(Loew)   LUMO(Mull)   LUMO(Loew)
 ---------------------------------------------------------------
   O         0.812345     0.798765     0.101234     0.098765
   H         0.093827     0.100617     0.449383     0.450617
 ---------------------------------------------------------------
"""


class TestNboCharges(unittest.TestCase):
    def test_the_summary_table_is_read(self):
        nbo = _parse(NBO_SUMMARY).get("NBO")
        self.assertEqual(len(nbo), 2)

    def test_summary_charges_are_read(self):
        nbo = _parse(NBO_SUMMARY).get("NBO")
        self.assertAlmostEqual(nbo[0]["charge"], -0.92345)
        self.assertAlmostEqual(nbo[1]["charge"], 0.46172)

    def test_summary_atom_indices_are_zero_based(self):
        nbo = _parse(NBO_SUMMARY).get("NBO")
        self.assertEqual([a["atom_idx"] for a in nbo], [0, 1])

    def test_summary_symbols_are_read(self):
        nbo = _parse(NBO_SUMMARY).get("NBO")
        self.assertEqual([a["atom_sym"] for a in nbo], ["O", "H"])

    def test_the_fallback_block_is_read_without_a_summary(self):
        nbo = _parse(NBO_FALLBACK).get("NBO")
        self.assertEqual(len(nbo), 3)

    def test_fallback_charges_are_read(self):
        nbo = _parse(NBO_FALLBACK).get("NBO")
        self.assertAlmostEqual(nbo[0]["charge"], -0.92345)

    def test_fallback_atom_indices_are_zero_based(self):
        nbo = _parse(NBO_FALLBACK).get("NBO")
        self.assertEqual([a["atom_idx"] for a in nbo], [0, 1, 2])

    def test_the_fallback_stops_at_the_next_section(self):
        # "Natural Electron Configuration" ends the block
        nbo = _parse(NBO_FALLBACK).get("NBO")
        self.assertEqual(len(nbo), 3)

    def test_unparseable_rows_are_skipped(self):
        text = NBO_FALLBACK.replace("      H    2    0.46172", "      H    x    junk")
        nbo = _parse(text).get("NBO")
        self.assertEqual(len(nbo), 2)

    def test_no_nbo_section_yields_no_nbo_charges(self):
        self.assertNotIn("NBO", _parse("SCF CONVERGED\n"))


class TestFmoCharges(unittest.TestCase):
    def test_the_table_is_read(self):
        fmo = _parse(FMO_BLOCK).get("FMO")
        self.assertEqual(len(fmo), 2)

    def test_both_gauges_are_recorded_for_the_homo(self):
        fmo = _parse(FMO_BLOCK).get("FMO")
        self.assertAlmostEqual(fmo[0]["homo_mulliken"], 0.812345)
        self.assertAlmostEqual(fmo[0]["homo_loewdin"], 0.798765)

    def test_both_gauges_are_recorded_for_the_lumo(self):
        fmo = _parse(FMO_BLOCK).get("FMO")
        self.assertAlmostEqual(fmo[0]["lumo_mulliken"], 0.101234)
        self.assertAlmostEqual(fmo[0]["lumo_loewdin"], 0.098765)

    def test_the_homo_population_doubles_as_the_colouring_value(self):
        fmo = _parse(FMO_BLOCK).get("FMO")
        self.assertAlmostEqual(fmo[0]["charge"], fmo[0]["homo_mulliken"])

    def test_atom_labels_are_split_into_index_and_symbol(self):
        fmo = _parse(FMO_BLOCK).get("FMO")
        self.assertEqual([a["atom_idx"] for a in fmo], [0, 1])
        self.assertEqual([a["atom_sym"] for a in fmo], ["O", "H"])

    def test_a_bare_symbol_label_falls_back_to_positional_indices(self):
        fmo = _parse(FMO_BARE_LABELS).get("FMO")
        self.assertEqual([a["atom_idx"] for a in fmo], [0, 1])
        self.assertEqual([a["atom_sym"] for a in fmo], ["O", "H"])

    def test_unparseable_rows_are_skipped(self):
        text = FMO_BLOCK.replace("0.093827     0.100617", "junk     junk")
        fmo = _parse(text).get("FMO")
        self.assertEqual(len(fmo), 1)

    def test_no_fmo_section_yields_no_fmo_charges(self):
        self.assertNotIn("FMO", _parse("SCF CONVERGED\n"))


if __name__ == "__main__":
    unittest.main()

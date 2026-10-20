"""
tests/test_parser_mixin_gaps.py
Targeted parser tests aimed at the uncovered lines reported in
parser_structure.py (86%) and parser_spectra.py (87%) after the mixin
extraction refactor.

All tests use a standalone OrcaParser loaded without Qt (via _parser_loader),
so they run in CI without any UI dependency.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from _parser_loader import load_standalone_parser  # noqa: E402

_mod = load_standalone_parser("orca_parser_mixin_gaps_test")
OrcaParser = _mod.OrcaParser


def _parse(text):
    p = OrcaParser()
    p.load_from_memory(text, "<test>")
    return p


# ---------------------------------------------------------------------------
# parser_structure.py  —  parse_xyz_content
# ---------------------------------------------------------------------------

_XYZ_HEADER = """\
                       * O   R   C   A *

Program Version 5.0.0
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O     0.0  0.0  0.0
  H     1.0  0.0  0.0

FINAL SINGLE POINT ENERGY      -76.0

                          ****ORCA TERMINATED NORMALLY****
"""


class TestParseXyzContent(unittest.TestCase):
    def _parser(self):
        p = OrcaParser()
        p.lines = []
        p.data = {"scan_steps": []}
        return p

    def test_single_frame_yields_one_step(self):
        p = self._parser()
        xyz = "2\nEnergy -76.0\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        self.assertEqual(len(steps), 1)

    def test_atom_symbols_are_captured(self):
        p = self._parser()
        xyz = "2\nE -76.0\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        self.assertEqual(steps[0]["atoms"], ["O", "H"])

    def test_energy_is_extracted_from_comment(self):
        p = self._parser()
        xyz = "2\nEnergy: -76.123\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        self.assertAlmostEqual(steps[0]["energy"], -76.123, places=3)

    def test_ts_frame_is_excluded(self):
        p = self._parser()
        xyz = "2\nTS E=-76.0\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        self.assertEqual(steps, [])

    def test_ci_frame_is_excluded(self):
        p = self._parser()
        xyz = "2\nCI E=-76.0\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        self.assertEqual(steps, [])

    def test_ci_neb_frame_is_kept(self):
        p = self._parser()
        xyz = "2\nCI-NEB E=-76.0\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        self.assertEqual(len(steps), 1)

    def test_scan_coordinate_in_comment_is_captured(self):
        p = self._parser()
        xyz = "2\nE=-76.0 Dist=1.234\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        self.assertAlmostEqual(steps[0]["dist"], 1.234, places=3)

    def test_malformed_coord_line_is_skipped_not_fatal(self):
        p = self._parser()
        xyz = "2\nE=-76.0\nO NOTANUMBER 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        # H should still be captured even though O was malformed
        self.assertEqual(len(steps), 1)
        self.assertIn("H", steps[0]["atoms"])
        self.assertNotIn("O", steps[0]["atoms"])

    def test_multi_frame_xyz_yields_multiple_steps(self):
        p = self._parser()
        frame = "2\nE=-76.0\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(frame * 3)
        self.assertEqual(len(steps), 3)

    def test_empty_content_returns_empty_list(self):
        p = self._parser()
        self.assertEqual(p.parse_xyz_content(""), [])

    def test_non_numeric_first_line_is_skipped(self):
        p = self._parser()
        xyz = "not-a-count\n2\nE=-76.0\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        self.assertEqual(len(steps), 1)

    def test_coordinate_keyword_in_comment_is_parsed(self):
        p = self._parser()
        xyz = "2\nCoordinate=2.5 E=-76.0\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
        steps = p.parse_xyz_content(xyz)
        self.assertAlmostEqual(steps[0]["scan_coord"], 2.5, places=3)


# ---------------------------------------------------------------------------
# parser_structure.py  —  parse_termination_status
# ---------------------------------------------------------------------------


class TestParseTerminationStatus(unittest.TestCase):
    def _parser_with_lines(self, lines):
        p = OrcaParser()
        p.lines = lines
        p.data = {}
        return p

    def test_normally_terminated_is_recognised(self):
        p = self._parser_with_lines(["****ORCA TERMINATED NORMALLY****"])
        p.parse_termination_status()
        self.assertIn("Terminated normally", p.data["termination_status"])

    def test_error_termination_is_recognised(self):
        p = self._parser_with_lines(["ORCA FINISHED BY ERROR TERMINATION"])
        p.parse_termination_status()
        self.assertEqual(p.data["termination_status"], "ERROR")

    def test_input_error_is_recognised(self):
        p = self._parser_with_lines(["INPUT ERROR: something went wrong"])
        p.parse_termination_status()
        self.assertEqual(p.data["termination_status"], "ERROR")

    def test_error_exclamation_is_recognised(self):
        p = self._parser_with_lines(["ERROR !!! bad"])
        p.parse_termination_status()
        self.assertEqual(p.data["termination_status"], "ERROR")

    def test_file_line_error_pattern_is_recognised(self):
        p = self._parser_with_lines(["[file foo.cpp, line 42]"])
        p.parse_termination_status()
        self.assertEqual(p.data["termination_status"], "ERROR")

    def test_still_running_if_no_terminal_marker(self):
        p = self._parser_with_lines(["SCF cycle 3 ..."])
        p.parse_termination_status()
        self.assertEqual(p.data["termination_status"], "Running")

    def test_empty_lines_does_not_crash(self):
        p = self._parser_with_lines([])
        p.parse_termination_status()  # must not raise
        self.assertEqual(p.data["termination_status"], "Running")


# ---------------------------------------------------------------------------
# parser_structure.py  —  parse_basic  (via full parse)
# ---------------------------------------------------------------------------


class TestParseBasic(unittest.TestCase):
    _BASE = """\
                       * O   R   C   A *

Program Version 5.0.4

Total charge 0
Multiplicity 1

CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O     0.0  0.0  0.0
  H     0.96 0.0  0.0

FINAL SINGLE POINT ENERGY      -76.0

                          ****ORCA TERMINATED NORMALLY****
"""

    def test_version_is_parsed(self):
        p = _parse(self._BASE)
        self.assertEqual(p.data["version"], "5.0.4")

    def test_scf_energy_is_parsed(self):
        p = _parse(self._BASE)
        self.assertAlmostEqual(p.data["scf_energy"], -76.0, places=2)

    def test_charge_is_parsed(self):
        p = _parse(self._BASE)
        self.assertEqual(p.data["charge"], 0)

    def test_multiplicity_is_parsed(self):
        p = _parse(self._BASE)
        self.assertEqual(p.data["mult"], 1)

    def test_atoms_are_parsed(self):
        p = _parse(self._BASE)
        self.assertIn("O", p.data["atoms"])
        self.assertIn("H", p.data["atoms"])

    def test_scan_flag_is_set(self):
        text = self._BASE.replace("FINAL", "RELAXED SURFACE SCAN\nFINAL")
        p = _parse(text)
        self.assertTrue(p.data.get("is_scan", False))

    def test_neb_flag_is_set(self):
        text = self._BASE.replace("FINAL", "NUDGED ELASTIC BAND\nFINAL")
        p = _parse(text)
        self.assertTrue(p.data.get("is_neb", False))

    def test_hurray_convergence_is_recognised(self):
        text = self._BASE + "\nHURRAY\n"
        p = _parse(text)
        self.assertTrue(p.data.get("converged", False))

    def test_scf_converged_sets_flag(self):
        text = self._BASE + "\nSCF CONVERGED\n"
        p = _parse(text)
        self.assertTrue(p.data.get("converged", False))

    def test_neb_trajectory_file_is_captured(self):
        text = (
            self._BASE + "\nCurrent trajectory will be written to ..... orca_neb.trj\n"
        )
        p = _parse(text)
        self.assertIn("neb_trj_file", p.data)


# ---------------------------------------------------------------------------
# parser_spectra.py  —  parse_nmr
# ---------------------------------------------------------------------------

_NMR_SHIELDING_BLOCK = """\
CHEMICAL SHIELDING SUMMARY (PPM)
---------------------------------
  N   Nucleus    Shielding
-----------------------------
  0   C         150.12
  1   H          31.45
  2   H          30.98

"""

_NMR_COUPLING_BLOCK = """\
SUMMARY OF ISOTROPIC COUPLING CONSTANTS (Hz)
-----------------------------------------
    0 C   1 H
---
0 C   0.000   7.500
1 H   7.500   0.000
"""


class TestParseNMR(unittest.TestCase):
    def test_shielding_values_extracted(self):
        text = _NMR_SHIELDING_BLOCK + "\n****ORCA TERMINATED NORMALLY****"
        p = _parse(text)
        shields = {
            e["atom_idx"]: e["shielding"] for e in p.data.get("nmr_shielding", [])
        }
        self.assertIn(0, shields)
        self.assertAlmostEqual(shields[0], 150.12, places=2)

    def test_shielding_atom_symbols_are_captured(self):
        text = _NMR_SHIELDING_BLOCK + "\n****ORCA TERMINATED NORMALLY****"
        p = _parse(text)
        syms = {e["atom_idx"]: e["atom_sym"] for e in p.data.get("nmr_shielding", [])}
        self.assertEqual(syms[1], "H")

    def test_all_shielding_atoms_captured(self):
        text = _NMR_SHIELDING_BLOCK + "\n****ORCA TERMINATED NORMALLY****"
        p = _parse(text)
        self.assertEqual(len(p.data["nmr_shielding"]), 3)

    def test_coupling_constant_is_extracted(self):
        text = _NMR_COUPLING_BLOCK + "\n****ORCA TERMINATED NORMALLY****"
        p = _parse(text)
        couplings = p.data.get("nmr_couplings", [])
        vals = [c["coupling"] for c in couplings]
        self.assertTrue(any(abs(v - 7.5) < 0.1 for v in vals))

    def test_empty_output_yields_empty_shielding_list(self):
        p = _parse("****ORCA TERMINATED NORMALLY****")
        self.assertEqual(p.data.get("nmr_shielding", []), [])


# ---------------------------------------------------------------------------
# parser_spectra.py  —  parse_tddft
# ---------------------------------------------------------------------------

_TDDFT_STATE_BLOCK = """\
STATE  1:  E=   3.456 eV   27880 cm**-1  358.8 nm
    77a ->  78a  :     0.70000 (c= 0.9899)

ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
----------------------------------------------------------
      State   Energy    Wavelength    fosc         T2
            (cm-1)      (nm)
----------------------------------------------------------
  0-1A  ->  1-1A   3.456   27880   358.8   0.1234    0.0
----------------------------------------------------------
"""


class TestParseTDDFT(unittest.TestCase):
    def test_tddft_states_are_extracted(self):
        text = _TDDFT_STATE_BLOCK + "\n****ORCA TERMINATED NORMALLY****"
        p = _parse(text)
        self.assertGreater(len(p.data.get("tddft", [])), 0)

    def test_state_energy_ev_is_captured(self):
        text = _TDDFT_STATE_BLOCK + "\n****ORCA TERMINATED NORMALLY****"
        p = _parse(text)
        states = p.data.get("tddft", [])
        energies = [s["energy_ev"] for s in states]
        self.assertTrue(any(abs(e - 3.456) < 0.01 for e in energies))

    def test_state_wavelength_is_captured(self):
        text = _TDDFT_STATE_BLOCK + "\n****ORCA TERMINATED NORMALLY****"
        p = _parse(text)
        states = p.data.get("tddft", [])
        nms = [s["energy_nm"] for s in states]
        self.assertTrue(any(abs(nm - 358.8) < 0.5 for nm in nms))

    def test_oscillator_strength_is_captured(self):
        text = _TDDFT_STATE_BLOCK + "\n****ORCA TERMINATED NORMALLY****"
        p = _parse(text)
        states = p.data.get("tddft", [])
        oscs = [s.get("osc_len", s.get("osc", 0.0)) for s in states]
        self.assertTrue(any(abs(o - 0.1234) < 0.001 for o in oscs))

    def test_transition_descriptions_are_captured(self):
        text = _TDDFT_STATE_BLOCK + "\n****ORCA TERMINATED NORMALLY****"
        p = _parse(text)
        states = p.data.get("tddft", [])
        trans = [t for s in states for t in s.get("transitions", [])]
        self.assertTrue(any("->" in t for t in trans))

    def test_empty_output_yields_empty_tddft_list(self):
        p = _parse("****ORCA TERMINATED NORMALLY****")
        self.assertEqual(p.data.get("tddft", []), [])


# ---------------------------------------------------------------------------
# parser_spectra.py  —  parse_frequencies  (edge cases)
# ---------------------------------------------------------------------------

_FREQ_BLOCK = """\
VIBRATIONAL FREQUENCIES
-----------------------

   0:         0.00 cm**-1
   1:         0.00 cm**-1
   2:         0.00 cm**-1
   3:      3656.78 cm**-1

"""

_IMAGINARY_FREQ_BLOCK = """\
VIBRATIONAL FREQUENCIES
-----------------------

   0:         0.00 cm**-1
   1:      -143.23 cm**-1
   2:      3656.78 cm**-1

"""


class TestParseFrequencies(unittest.TestCase):
    def test_real_frequencies_are_extracted(self):
        p = _parse(_FREQ_BLOCK + "****ORCA TERMINATED NORMALLY****")
        # frequencies is a list of dicts: {"freq": float, "ir": float, ...}
        freq_entries = p.data.get("frequencies", [])
        self.assertGreater(len(freq_entries), 0)
        freq_vals = [e["freq"] for e in freq_entries]
        self.assertIn(3656.78, freq_vals)

    def test_imaginary_frequency_is_detected_via_negative_value(self):
        p = _parse(_IMAGINARY_FREQ_BLOCK + "****ORCA TERMINATED NORMALLY****")
        freq_entries = p.data.get("frequencies", [])
        freq_vals = [e["freq"] for e in freq_entries]
        # Imaginary mode is stored as a negative value
        self.assertTrue(any(v < -10 for v in freq_vals))

    def test_zero_frequencies_are_tolerated(self):
        p = _parse(_FREQ_BLOCK + "****ORCA TERMINATED NORMALLY****")
        freq_entries = p.data.get("frequencies", [])
        freq_vals = [e["freq"] for e in freq_entries]
        self.assertIn(0.0, freq_vals)


# ---------------------------------------------------------------------------
# parser_structure.py  —  parse_gradient (via full parse)
# ---------------------------------------------------------------------------

_GRAD_BLOCK = """\
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O   0.0  0.0  0.0
  H   1.0  0.0  0.0

CARTESIAN GRADIENT
------------------
   1   O   :    0.0001    0.0002    0.0003
   2   H   :   -0.0001   -0.0002   -0.0003

FINAL SINGLE POINT ENERGY   -76.0

****ORCA TERMINATED NORMALLY****
"""


class TestParseGradient(unittest.TestCase):
    def test_gradient_entries_are_extracted(self):
        p = _parse(_GRAD_BLOCK)
        # gradient data stored under 'gradients', each entry has 'vector'
        grad = p.data.get("gradients", [])
        self.assertGreater(len(grad), 0)

    def test_gradient_has_correct_vector_length(self):
        p = _parse(_GRAD_BLOCK)
        grad = p.data.get("gradients", [])
        for g in grad:
            self.assertEqual(len(g["vector"]), 3)

    def test_first_gradient_values_match(self):
        p = _parse(_GRAD_BLOCK)
        grad = p.data.get("gradients", [])
        self.assertAlmostEqual(grad[0]["vector"][0], 0.0001, places=4)


if __name__ == "__main__":
    unittest.main()

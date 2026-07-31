"""
Regression tests for the v3.11.0 chemistry fixes.

Each of these produced plausible-looking but wrong numbers rather than an
error, so they are pinned against real ORCA output where possible.
"""

import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.gui_harness import load_isolated  # noqa: E402

SAMPLES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_outputs")


def _parse(filename):
    from orca_result_analyzer_rust.parser import OrcaParser

    path = os.path.join(SAMPLES, filename)
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        content = fh.read()
    parser = OrcaParser()
    parser.load_from_memory(content, path)
    return parser


# ---------------------------------------------------------------------------
# Dipole: components are a.u., magnitude is Debye
# ---------------------------------------------------------------------------


class TestDipoleUnits(unittest.TestCase):
    def test_acetone_dipole_matches_experiment_in_debye(self):
        """Acetone's measured gas-phase dipole is 2.88 D. The a.u. vector was
        being displayed under a "Debye" label, which reads as 1.06 D."""
        d = _parse("acetone-opt.out").data["dipoles"]
        self.assertAlmostEqual(d["magnitude_debye"], 2.7871, places=3)
        self.assertGreater(d["magnitude_debye"], 2.5)
        self.assertLess(d["magnitude_debye"], 3.1)

    def test_debye_vector_norm_equals_debye_magnitude(self):
        """The components and the magnitude must be on the same scale."""
        d = _parse("acetone-opt.out").data["dipoles"]
        norm = math.sqrt(sum(c * c for c in d["vector_debye"]))
        self.assertAlmostEqual(norm, d["magnitude_debye"], places=3)

    def test_au_vector_norm_equals_au_magnitude(self):
        d = _parse("acetone-opt.out").data["dipoles"]
        norm = math.sqrt(sum(c * c for c in d["vector_au"]))
        self.assertAlmostEqual(norm, d["magnitude_au"], places=3)

    def test_au_and_debye_differ_by_the_conversion_constant(self):
        """ORCA carries its own a.u.->Debye constant (2.54184) which differs
        from CODATA (2.541746) in the fifth digit. The parser scales the
        vector by ORCA's own ratio so the components stay consistent with the
        magnitude ORCA printed, so only assert it is that conversion."""
        d = _parse("acetone-opt.out").data["dipoles"]
        ratio = d["magnitude_debye"] / d["magnitude_au"]
        self.assertAlmostEqual(ratio, 2.541746473, places=3)

    def test_benzene_dipole_is_zero_by_symmetry(self):
        d = _parse("benzene-opt.out").data["dipoles"]
        self.assertLess(d["magnitude_debye"], 1e-3)


# ---------------------------------------------------------------------------
# IR / Raman intensity columns
# ---------------------------------------------------------------------------


class TestVibrationalIntensities(unittest.TestCase):
    def test_acetone_ir_intensities_are_populated(self):
        freqs = _parse("acetone-opt.out").data["frequencies"]
        self.assertEqual(len(freqs), 30)  # 10 atoms -> 3N
        active = [f for f in freqs if f["ir"] > 0]
        self.assertGreater(len(active), 10)

    def test_acetone_carbonyl_stretch_is_the_strongest_band(self):
        """Acetone's C=O stretch near 1700-1800 cm-1 dominates its IR
        spectrum; picking the wrong column produced a flat or scrambled one."""
        freqs = _parse("acetone-opt.out").data["frequencies"]
        strongest = max(freqs, key=lambda f: f["ir"])
        self.assertGreater(strongest["ir"], 100.0)
        self.assertGreater(strongest["freq"], 1600.0)
        self.assertLess(strongest["freq"], 1900.0)

    def test_the_six_trans_rot_modes_carry_no_intensity(self):
        freqs = _parse("acetone-opt.out").data["frequencies"]
        for f in freqs[:6]:
            self.assertEqual(f["ir"], 0.0)

    def test_intensity_column_is_read_from_the_header(self):
        from orca_result_analyzer_rust.parser import OrcaParser

        parser = OrcaParser()
        parser.lines = [
            "IR SPECTRUM",
            "-----------",
            "",
            " Mode   freq       eps      Int      T**2        TX",
            "       cm**-1   L/(mol*cm) km/mol    a.u.",
        ]
        self.assertEqual(parser._find_intensity_column(0, ("int", "t**2"), 99), 3)

    def test_missing_header_falls_back_to_the_default_column(self):
        from orca_result_analyzer_rust.parser import OrcaParser

        parser = OrcaParser()
        parser.lines = ["IR SPECTRUM", "-----------", "", "", ""]
        self.assertEqual(parser._find_intensity_column(0, ("int",), 3), 3)


# ---------------------------------------------------------------------------
# Imaginary frequency counting
# ---------------------------------------------------------------------------


class TestImaginaryFrequencyCount(unittest.TestCase):
    def test_optimized_minimum_reports_no_imaginary_modes(self):
        for name in ("acetone-opt.out", "benzene-opt.out"):
            data = _parse(name).data
            self.assertEqual(data["thermal"]["imaginary_freq_count"], 0, name)

    def test_small_negative_noise_is_not_counted(self):
        """Unprojected rotations arrive as tiny negative values; counting them
        reported a transition state for a perfectly good minimum."""
        from orca_result_analyzer_rust.parser import IMAGINARY_FREQ_THRESHOLD

        self.assertGreater(IMAGINARY_FREQ_THRESHOLD, 0.0)
        noise = [-0.02, -3.1, 0.0, 120.0]
        counted = [f for f in noise if f < -IMAGINARY_FREQ_THRESHOLD]
        self.assertEqual(counted, [])

    def test_a_genuine_imaginary_mode_is_still_counted(self):
        from orca_result_analyzer_rust.parser import IMAGINARY_FREQ_THRESHOLD

        modes = [-521.4, -0.01, 0.0, 300.0]
        counted = [f for f in modes if f < -IMAGINARY_FREQ_THRESHOLD]
        self.assertEqual(counted, [-521.4])


# ---------------------------------------------------------------------------
# NMR reference standards
# ---------------------------------------------------------------------------


class TestNMRReferenceStandards(unittest.TestCase):
    def setUp(self):
        self.M = load_isolated("nmr_analysis")

    def test_every_reference_of_a_nucleus_gives_the_same_shift(self):
        """delta = delta_ref + (sigma_ref - sigma), so delta_ref + sigma_ref is
        the absolute shielding of the delta=0 point and must be identical
        across a nucleus's reference compounds. 15N's two entries used to
        disagree by 1.7 ppm, so switching reference moved every peak."""
        for nucleus, refs in self.M.DEFAULT_REFERENCE_STANDARDS.items():
            anchors = {
                name: r["delta_ref"] + r["sigma_ref"]
                for name, r in refs.items()
                if name != "No Reference"
            }
            self.assertTrue(anchors, nucleus)
            spread = max(anchors.values()) - min(anchors.values())
            self.assertLess(spread, 1e-6, f"{nucleus} references disagree: {anchors}")

    def test_a_sample_gets_the_same_shift_from_cdcl3_and_tms(self):
        refs = self.M.DEFAULT_REFERENCE_STANDARDS["13C"]
        sigma_sample = 50.0
        tms = refs["TMS"]
        cdcl3 = refs["CDCl3"]
        d_tms = tms["delta_ref"] + (tms["sigma_ref"] - sigma_sample)
        d_cdcl3 = cdcl3["delta_ref"] + (cdcl3["sigma_ref"] - sigma_sample)
        self.assertAlmostEqual(d_tms, d_cdcl3, places=6)

    def test_benzene_carbon_lands_near_its_known_shift(self):
        """TMS-referenced: benzene 13C is 128.4 ppm, so sigma ~ 54 ppm."""
        tms = self.M.DEFAULT_REFERENCE_STANDARDS["13C"]["TMS"]
        delta = tms["delta_ref"] + (tms["sigma_ref"] - 54.0)
        self.assertAlmostEqual(delta, 128.4, delta=1.0)

    def test_the_table_is_defined_once(self):
        """It used to be duplicated verbatim in __init__ and save_settings, so
        editing one silently reclassified a built-in as a custom reference."""
        source_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "orca_result_analyzer_rust",
            "nmr_analysis.py",
        )
        with open(source_path, "r", encoding="utf-8") as fh:
            source = fh.read()
        # Each built-in shielding value is written exactly once, in the module
        # constant; save_settings reads that constant instead of restating it.
        for literal in (
            '"sigma_ref": 31.80',
            '"sigma_ref": 182.40',
            '"sigma_ref": 328.40',
        ):
            self.assertEqual(source.count(literal), 1, literal)
        self.assertIn("default_standards = DEFAULT_REFERENCE_STANDARDS", source)


if __name__ == "__main__":
    unittest.main()

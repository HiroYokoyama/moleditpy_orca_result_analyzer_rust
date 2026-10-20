"""
tests/test_parser_extended.py
Extended parser tests — thermal, orbital energies, IR/Raman spectra,
Mayer charges, scan results table, and recalc_energies.

All tests load OrcaParser directly (no Qt stubs required).
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from _parser_loader import load_standalone_parser  # noqa: E402


def _load_parser():
    return load_standalone_parser("orca_parser_ext")


_parser_mod = _load_parser()
OrcaParser = _parser_mod.OrcaParser


def _parse_method(content, method_name):
    p = OrcaParser()
    p.raw_content = content
    p.lines = content.splitlines()
    p.filename = "test.out"
    getattr(p, method_name)()
    return p


# ---------------------------------------------------------------------------
# TestParseThermal
# ---------------------------------------------------------------------------

# ORCA thermochemistry block.
# parse_thermal looks for "THERMOCHEMISTRY AT" then scans for matching
# keywords; if "Eh" is in the line it extracts the last float before "Eh".
# Temperature line must have "TEMPERATURE", "K", and "..." all present.

_THERMAL_CONTENT = """\
THERMOCHEMISTRY AT 298.15 K
-------------------------------
Temperature   ...   298.15 K
Electronic Energy          ...   -76.12345 Eh
Zero Point Energy          ...    0.02345 Eh
Total Thermal Energy       ...   -76.05678 Eh
Total Enthalpy             ...   -76.04567 Eh
Final Entropy Term         ...   -0.02345 Eh
Final Gibbs Free Energy    ...   -76.04321 Eh
"""

_THERMAL_MINIMAL = """\
THERMOCHEMISTRY AT 298.15 K
Electronic Energy          ...   -100.0 Eh
Final Gibbs Free Energy    ...    -99.5 Eh
"""


class TestParseThermal(unittest.TestCase):
    def test_electronic_energy(self):
        p = _parse_method(_THERMAL_CONTENT, "parse_thermal")
        self.assertAlmostEqual(
            p.data["thermal"]["electronic_energy"], -76.12345, places=4
        )

    def test_zero_point_energy(self):
        p = _parse_method(_THERMAL_CONTENT, "parse_thermal")
        self.assertAlmostEqual(p.data["thermal"]["zpe"], 0.02345, places=4)

    def test_total_enthalpy(self):
        p = _parse_method(_THERMAL_CONTENT, "parse_thermal")
        self.assertAlmostEqual(p.data["thermal"]["enthalpy"], -76.04567, places=4)

    def test_gibbs_free_energy(self):
        p = _parse_method(_THERMAL_CONTENT, "parse_thermal")
        self.assertAlmostEqual(p.data["thermal"]["gibbs"], -76.04321, places=4)

    def test_entropy_term(self):
        p = _parse_method(_THERMAL_CONTENT, "parse_thermal")
        self.assertAlmostEqual(p.data["thermal"]["entropy"], -0.02345, places=4)

    def test_temperature_parsed(self):
        p = _parse_method(_THERMAL_CONTENT, "parse_thermal")
        self.assertAlmostEqual(
            p.data["thermal"].get("temperature", 0.0), 298.15, places=2
        )

    def test_minimal_gibbs(self):
        p = _parse_method(_THERMAL_MINIMAL, "parse_thermal")
        self.assertAlmostEqual(p.data["thermal"]["gibbs"], -99.5, places=1)

    def test_empty_content(self):
        p = _parse_method("", "parse_thermal")
        self.assertEqual(p.data["thermal"], {})


# ---------------------------------------------------------------------------
# TestParseOrbitalEnergies
# ---------------------------------------------------------------------------

# Format (restricted RHF):
#   ORBITAL ENERGIES
#   ----------------
#     NO   OCC          E(Eh)            E(eV)
#      0   2.0000     -20.12345      -547.89012
#      5   0.0000       0.34567         9.40346

_ORBITAL_RESTRICTED = """\
ORBITAL ENERGIES
----------------
   NO   OCC          E(Eh)            E(eV)
    0   2.0000     -20.12345      -547.89012
    1   2.0000      -1.23456       -33.59876
    2   0.0000       0.34567         9.40346
---
"""

# UHF: two sections
_ORBITAL_UHF = """\
SPIN UP ORBITALS
----------------
   NO   OCC          E(Eh)            E(eV)
    0   1.0000     -20.12345      -547.89012
    1   0.0000       0.34567         9.40346
---
SPIN DOWN ORBITALS
------------------
   NO   OCC          E(Eh)            E(eV)
    0   1.0000     -20.11234      -547.58765
    1   0.0000       0.36789        10.01234
---
"""


# UHF as real ORCA actually prints it: an "ORBITAL ENERGIES" header followed
# by "SPIN UP/DOWN ORBITALS" sub-blocks whose column header sits on the very
# next line (NO "---" separator), with a blank line between the two spins.
_ORBITAL_UHF_REAL = """\
ORBITAL ENERGIES
----------------
                 SPIN UP ORBITALS
  NO   OCC          E(Eh)            E(eV)
    0   1.0000     -19.36210      -526.86960
    1   1.0000      -1.34548       -36.61250
    2   0.0000       0.25639         6.97690

                 SPIN DOWN ORBITALS
  NO   OCC          E(Eh)            E(eV)
    0   1.0000     -19.33005      -525.99760
    1   0.0000       0.30000         8.16340
---
"""


class TestParseOrbitalEnergies(unittest.TestCase):
    def test_restricted_count(self):
        p = _parse_method(_ORBITAL_RESTRICTED, "parse_orbital_energies")
        self.assertEqual(len(p.data["orbital_energies"]), 3)

    def test_restricted_energies(self):
        p = _parse_method(_ORBITAL_RESTRICTED, "parse_orbital_energies")
        mos = p.data["orbital_energies"]
        self.assertAlmostEqual(mos[0]["energy_eh"], -20.12345, places=4)
        self.assertAlmostEqual(mos[0]["energy_ev"], -547.89012, places=4)
        self.assertAlmostEqual(mos[2]["energy_eh"], 0.34567, places=4)

    def test_occupation_values(self):
        p = _parse_method(_ORBITAL_RESTRICTED, "parse_orbital_energies")
        mos = p.data["orbital_energies"]
        self.assertAlmostEqual(mos[0]["occupation"], 2.0)
        self.assertAlmostEqual(mos[2]["occupation"], 0.0)

    def test_occupied_vs_virtual_type(self):
        p = _parse_method(_ORBITAL_RESTRICTED, "parse_orbital_energies")
        mos = p.data["orbital_energies"]
        self.assertEqual(mos[0]["type"], "occupied")
        self.assertEqual(mos[2]["type"], "virtual")

    def test_spin_label_restricted(self):
        p = _parse_method(_ORBITAL_RESTRICTED, "parse_orbital_energies")
        for mo in p.data["orbital_energies"]:
            self.assertEqual(mo["spin"], "restricted")

    def test_uhf_both_spins(self):
        p = _parse_method(_ORBITAL_UHF, "parse_orbital_energies")
        spins = {mo["spin"] for mo in p.data["orbital_energies"]}
        self.assertIn("alpha", spins)
        self.assertIn("beta", spins)

    def test_uhf_total_count(self):
        p = _parse_method(_ORBITAL_UHF, "parse_orbital_energies")
        self.assertEqual(len(p.data["orbital_energies"]), 4)

    def test_backward_compat_mos_list(self):
        p = _parse_method(_ORBITAL_RESTRICTED, "parse_orbital_energies")
        self.assertEqual(len(p.data["mos"]), len(p.data["orbital_energies"]))

    def test_empty_content(self):
        p = _parse_method("", "parse_orbital_energies")
        self.assertEqual(p.data["orbital_energies"], [])

    def test_uhf_real_format_both_spins(self):
        """Real ORCA UHF layout (header on next line, no '---') must split spins."""
        p = _parse_method(_ORBITAL_UHF_REAL, "parse_orbital_energies")
        spins = {mo["spin"] for mo in p.data["orbital_energies"]}
        self.assertIn("alpha", spins)
        self.assertIn("beta", spins)
        self.assertNotIn("restricted", spins)

    def test_uhf_real_format_counts(self):
        """Alpha block (3 orbitals) and beta block (2 orbitals) must both be read."""
        p = _parse_method(_ORBITAL_UHF_REAL, "parse_orbital_energies")
        alpha = [o for o in p.data["orbital_energies"] if o["spin"] == "alpha"]
        beta = [o for o in p.data["orbital_energies"] if o["spin"] == "beta"]
        self.assertEqual(len(alpha), 3)
        self.assertEqual(len(beta), 2)

    def test_uhf_real_format_beta_energies(self):
        """Beta manifold must carry its own energies, not be dropped."""
        p = _parse_method(_ORBITAL_UHF_REAL, "parse_orbital_energies")
        beta = [o for o in p.data["orbital_energies"] if o["spin"] == "beta"]
        self.assertAlmostEqual(beta[0]["energy_eh"], -19.33005, places=4)
        self.assertAlmostEqual(beta[1]["energy_eh"], 0.30000, places=4)


# ---------------------------------------------------------------------------
# TestParseFrequenciesIR
# ---------------------------------------------------------------------------

# IR SPECTRUM block starts at ir_start; parser begins reading at ir_start+5.
# Data line format: "  N:   freq   mass   T^2   TX   TY   TZ"
# parser reads parts[0]="N:" parts[3]=intensity

# parse_frequencies IR/Raman: `curr = ir_start + 5` so exactly 5 lines must
# separate "IR SPECTRUM" from the first data row.  The standard ORCA layout is:
#   ir_start+0  IR SPECTRUM
#   ir_start+1  -----------
#   ir_start+2  (blank)
#   ir_start+3  Mode  freq  mass  T^2  TX  TY  TZ
#   ir_start+4  -----------
#   ir_start+5  <first data row>  ← loop starts here

_FREQ_WITH_IR = """\
VIBRATIONAL FREQUENCIES
-----------------------
   0:         0.00 cm**-1
   1:         0.00 cm**-1
   2:         0.00 cm**-1
   3:      1623.45 cm**-1
   4:      3652.84 cm**-1
   5:      3754.93 cm**-1

IR SPECTRUM
-----------

Mode   freq   mass   T^2    TX     TY     TZ
-----------
    3:  1623.45  1.234  25.678  0.00  0.00  5.07
    4:  3652.84  1.123   8.901  0.00  0.00  2.98
    5:  3754.93  1.089   3.456  0.00  0.00  1.86
"""

# As ORCA 5 writes it: the transition dipole derivative in parentheses after
# the intensity columns. 5.07^2 = 25.7049, matching the T**2 column.
_FREQ_WITH_IR_DERIV = """\
VIBRATIONAL FREQUENCIES
-----------------------
   0:         0.00 cm**-1
   1:         0.00 cm**-1
   2:         0.00 cm**-1
   3:      1623.45 cm**-1

IR SPECTRUM
-----------

 Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
    3:  1623.45   0.005000   25.678  25.7049  ( 0.000000  0.000000  5.070000)
"""

_FREQ_WITH_RAMAN = """\
VIBRATIONAL FREQUENCIES
-----------------------
   0:         0.00 cm**-1
   1:         0.00 cm**-1
   2:         0.00 cm**-1
   3:      1623.45 cm**-1
   4:      3652.84 cm**-1

RAMAN SPECTRUM
--------------

Mode   freq   Activity   DepolarP
--------------
    3:  1623.45   45.678   0.750
    4:  3652.84   12.345   0.750
"""


class TestParseFrequenciesIR(unittest.TestCase):
    def test_ir_intensities_assigned(self):
        p = _parse_method(_FREQ_WITH_IR, "parse_frequencies")
        # Mode 3 has IR intensity 25.678
        self.assertAlmostEqual(p.data["frequencies"][3]["ir"], 25.678, places=2)

    def test_ir_multiple_modes(self):
        p = _parse_method(_FREQ_WITH_IR, "parse_frequencies")
        self.assertAlmostEqual(p.data["frequencies"][4]["ir"], 8.901, places=2)
        self.assertAlmostEqual(p.data["frequencies"][5]["ir"], 3.456, places=2)

    def test_modes_without_ir_have_zero(self):
        p = _parse_method(_FREQ_WITH_IR, "parse_frequencies")
        self.assertAlmostEqual(p.data["frequencies"][0]["ir"], 0.0)

    def test_raman_intensities(self):
        p = _parse_method(_FREQ_WITH_RAMAN, "parse_frequencies")
        self.assertAlmostEqual(p.data["frequencies"][3]["raman"], 45.678, places=2)
        self.assertAlmostEqual(p.data["frequencies"][4]["raman"], 12.345, places=2)


class TestParseDipoleDerivative(unittest.TestCase):
    """The (TX TY TZ) group each IR row carries, in atomic units."""

    def test_the_triple_is_read_from_the_parentheses(self):
        line = "  10:    537.95   0.002935   14.83  0.001702  ( 0.040181 -0.009108  0.002223)"
        self.assertEqual(
            OrcaParser._parse_dipole_derivative(line),
            (0.040181, -0.009108, 0.002223),
        )

    def test_a_row_without_the_group_is_not_an_error(self):
        """Older ORCA prints the intensity columns and stops."""
        self.assertIsNone(
            OrcaParser._parse_dipole_derivative("    3:  1623.45  1.234  25.678")
        )

    def test_a_malformed_group_is_rejected_rather_than_guessed(self):
        for bad in (
            "  1:  100.0  ( 0.1 0.2 )",  # only two components
            "  1:  100.0  ( 0.1 0.2 abc)",  # not a number
        ):
            with self.subTest(line=bad):
                self.assertIsNone(OrcaParser._parse_dipole_derivative(bad))

    def test_the_derivative_lands_on_the_mode(self):
        p = _parse_method(_FREQ_WITH_IR_DERIV, "parse_frequencies")
        self.assertEqual(p.data["frequencies"][3]["dipole_deriv"], (0.0, 0.0, 5.07))

    def test_its_square_is_the_t2_column_orca_printed(self):
        """ORCA's T**2 is |dmu/dQ|^2, so the parse can be checked against
        the program's own arithmetic rather than against itself."""
        p = _parse_method(_FREQ_WITH_IR_DERIV, "parse_frequencies")
        mode = p.data["frequencies"][3]
        self.assertAlmostEqual(
            sum(c * c for c in mode["dipole_deriv"]), 25.7049, places=4
        )

    def test_modes_outside_the_ir_table_carry_no_derivative(self):
        p = _parse_method(_FREQ_WITH_IR_DERIV, "parse_frequencies")
        self.assertIsNone(p.data["frequencies"][0].get("dipole_deriv"))


# ---------------------------------------------------------------------------
# TestSpectrumColumnsRealHeaders
# ---------------------------------------------------------------------------

# The headers above are simplified. Real ORCA writes the frequency unit as an
# inline "(cm**-1)" token in the Raman header (and in the ORCA 4 IR header),
# which has no counterpart column in the data rows. Counting it shifted every
# lookup one place right: Raman activity was read from the depolarization
# column, and the ORCA 4 IR intensity from TX.

_REAL_FREQ_HEADER = """\
-----------------------
VIBRATIONAL FREQUENCIES
-----------------------

   0:         0.00 cm**-1
   1:         0.00 cm**-1
   2:         0.00 cm**-1
   3:         0.00 cm**-1
   4:         0.00 cm**-1
   5:         0.00 cm**-1
   6:      1623.45 cm**-1
   7:      3652.84 cm**-1
"""

_REAL_ORCA5_SPECTRA = (
    _REAL_FREQ_HEADER
    + """
-----------
IR SPECTRUM
-----------

 Mode   freq       eps      Int      T**2         TX        TY        TZ
        cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
    6:   1623.45   1.9999   66.66   0.002999   (-0.05  0.01  0.00)
    7:   3652.84   0.5000   11.11   0.000500   ( 0.02  0.00  0.00)

--------------
RAMAN SPECTRUM
--------------

 Mode    freq (cm**-1)   Activity   Depolarization
-------------------------------------------------
    6:      1623.45       45.678       0.750
    7:      3652.84       12.345       0.111
"""
)

_REAL_ORCA4_IR = (
    _REAL_FREQ_HEADER
    + """
-----------
IR SPECTRUM
-----------

 Mode    freq (cm**-1)   T**2         TX         TY         TZ
                          (KM/Mole)   (a.u.)     (a.u.)     (a.u.)
-------------------------------------------------------------------
    6:   1623.45   66.660000  ( -0.050000   0.010000   0.000000)
    7:   3652.84   11.110000  (  0.020000   0.000000   0.000000)
"""
)


class TestSpectrumColumnsRealHeaders(unittest.TestCase):
    def test_orca5_ir_intensity(self):
        p = _parse_method(_REAL_ORCA5_SPECTRA, "parse_frequencies")
        self.assertAlmostEqual(p.data["frequencies"][6]["ir"], 66.66, places=2)
        self.assertAlmostEqual(p.data["frequencies"][7]["ir"], 11.11, places=2)

    def test_orca5_raman_reads_activity_not_depolarization(self):
        p = _parse_method(_REAL_ORCA5_SPECTRA, "parse_frequencies")
        self.assertAlmostEqual(p.data["frequencies"][6]["raman"], 45.678, places=3)
        self.assertAlmostEqual(p.data["frequencies"][7]["raman"], 12.345, places=3)

    def test_orca4_ir_intensity_not_tx(self):
        p = _parse_method(_REAL_ORCA4_IR, "parse_frequencies")
        self.assertAlmostEqual(p.data["frequencies"][6]["ir"], 66.66, places=2)
        self.assertAlmostEqual(p.data["frequencies"][7]["ir"], 11.11, places=2)

    def test_unit_tokens_ignored_in_header_lookup(self):
        p = OrcaParser()
        p.lines = [
            "RAMAN SPECTRUM",
            " Mode    freq (cm**-1)   Activity   Depolarization",
        ]
        self.assertEqual(p._find_intensity_column(0, ("activity",), default=99), 2)


# ---------------------------------------------------------------------------
# TestParseChargesMayer
# ---------------------------------------------------------------------------

# Mayer population analysis format:
#   MAYER POPULATION ANALYSIS
#   ...
#   ATOM    NA        ZA        QA        VA       BVA       FA
#      0 C     6.0000    6.0000   -0.12345   4.01234   3.98765   0.02469

_MAYER_CONTENT = """\
MAYER POPULATION ANALYSIS
---
ATOM    NA        ZA        QA        VA       BVA       FA
   0 C     6.0000    6.0000   -0.12345   4.01234   3.98765   0.02469
   1 H     0.8765    1.0000    0.06172   0.98765   0.97654   0.01111
   2 H     0.8765    1.0000    0.06173   0.98765   0.97654   0.01111
Mayer bond orders
"""

_MAYER_NO_VALENCY = """\
MAYER POPULATION ANALYSIS
---
ATOM    NA    ZA    QA
   0 O     8.0000  8.0000  -0.76543
Mayer bond orders
"""


class TestParseChargesMayer(unittest.TestCase):
    def test_mayer_count(self):
        p = _parse_method(_MAYER_CONTENT, "parse_charges")
        self.assertEqual(len(p.data["charges"]["Mayer"]), 3)

    def test_mayer_charge_values(self):
        p = _parse_method(_MAYER_CONTENT, "parse_charges")
        entry = p.data["charges"]["Mayer"][0]
        self.assertEqual(entry["atom_sym"], "C")
        self.assertAlmostEqual(entry["charge"], -0.12345, places=4)

    def test_mayer_valency_parsed(self):
        p = _parse_method(_MAYER_CONTENT, "parse_charges")
        entry = p.data["charges"]["Mayer"][0]
        self.assertAlmostEqual(entry["valency"], 4.01234, places=4)
        self.assertAlmostEqual(entry["bonded_valency"], 3.98765, places=4)
        self.assertAlmostEqual(entry["free_valency"], 0.02469, places=4)

    def test_mayer_fallback_as_mulliken_when_no_mulliken(self):
        """When Mulliken is absent, Mayer data is copied into Mulliken key."""
        p = _parse_method(_MAYER_NO_VALENCY, "parse_charges")
        # No separate Mulliken block, so Mayer is also stored as Mulliken
        self.assertIn("Mulliken", p.data["charges"])
        self.assertAlmostEqual(
            p.data["charges"]["Mulliken"][0]["charge"], -0.76543, places=4
        )

    def test_empty_content(self):
        p = _parse_method("", "parse_charges")
        self.assertNotIn("Mayer", p.data["charges"])

    def test_a_truncated_row_does_not_take_out_the_others(self):
        """QA is parts[4]; a >= 4 guard let a short row reach an IndexError
        that the blanket handler swallowed, silently losing the row."""
        content = (
            "MAYER POPULATION ANALYSIS\n"
            "---\n"
            "ATOM    NA    ZA    QA\n"
            "   0 C     6.0000  6.0000  -0.12345\n"
            "   1 H     0.8765\n"
            "   2 H     0.8765  1.0000   0.06172\n"
            "Mayer bond orders\n"
        )
        p = _parse_method(content, "parse_charges")
        idxs = [e["atom_idx"] for e in p.data["charges"]["Mayer"]]
        self.assertEqual(idxs, [0, 2])


# ---------------------------------------------------------------------------
# TestParseScanResultsTable
# ---------------------------------------------------------------------------

# parse_scan_results_table searches reversed lines for "Actual Energy",
# then reads coord/energy pairs from lines immediately following it.

_SCAN_RESULTS_CONTENT = """\
The Calculated Surface
Actual Energy
  1.000000    -76.10000
  1.200000    -76.20000
  1.400000    -76.30000

"""

_SCAN_RESULTS_WITH_EXISTING_STEPS = """\
RELAXED SURFACE SCAN STEP   1
  Actual scan coordinate      ...   1.000000
  FINAL SINGLE POINT ENERGY    -76.10000000
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  C      0.00000    0.00000    0.00000
  H      1.00000    0.00000    0.00000

RELAXED SURFACE SCAN STEP   2
  Actual scan coordinate      ...   1.200000
  FINAL SINGLE POINT ENERGY    -76.20000000
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  C      0.00000    0.00000    0.00000
  H      1.20000    0.00000    0.00000

Actual Energy
  1.000000    -76.10000
  1.200000    -76.20000

"""


class TestParseScanResultsTable(unittest.TestCase):
    def test_builds_steps_from_summary(self):
        p = _parse_method(_SCAN_RESULTS_CONTENT, "parse_scan_results_table")
        self.assertEqual(len(p.data["scan_steps"]), 3)

    def test_energies_from_summary(self):
        p = _parse_method(_SCAN_RESULTS_CONTENT, "parse_scan_results_table")
        energies = [s["energy"] for s in p.data["scan_steps"]]
        self.assertAlmostEqual(energies[0], -76.10000, places=4)
        self.assertAlmostEqual(energies[2], -76.30000, places=4)

    def test_coords_from_summary(self):
        p = _parse_method(_SCAN_RESULTS_CONTENT, "parse_scan_results_table")
        coords = [s["scan_coord"] for s in p.data["scan_steps"]]
        self.assertAlmostEqual(coords[0], 1.0, places=4)
        self.assertAlmostEqual(coords[1], 1.2, places=4)
        self.assertAlmostEqual(coords[2], 1.4, places=4)

    def test_step_type_summary(self):
        p = _parse_method(_SCAN_RESULTS_CONTENT, "parse_scan_results_table")
        for step in p.data["scan_steps"]:
            self.assertEqual(step["type"], "scan_step_summary")

    def test_maps_coords_to_existing_steps(self):
        """When scan_steps already exist, scan_coords are mapped into them."""
        p = OrcaParser()
        p.raw_content = _SCAN_RESULTS_WITH_EXISTING_STEPS
        p.lines = _SCAN_RESULTS_WITH_EXISTING_STEPS.splitlines()
        p.parse_trajectory()  # builds existing scan_steps
        p.parse_scan_results_table()  # maps coords

        for step in p.data["scan_steps"]:
            self.assertIsNotNone(step.get("scan_coord"))

    def test_empty_content(self):
        p = _parse_method("", "parse_scan_results_table")
        self.assertEqual(p.data.get("scan_steps", []), [])


# ---------------------------------------------------------------------------
# TestParseTrajectoryNeb
# ---------------------------------------------------------------------------

# NEB PATH SUMMARY: looks for "PATH SUMMARY" after a "----" separator.
# Parser expects: header line with "Image" and "E(Eh)", then data rows.

_NEB_PATH_SUMMARY = """\
----------------------------
         PATH SUMMARY
----------------------------
All forces in Eh/Bohr.
Image Dist.(Ang.)    E(Eh)   dE(kcal/mol)  max(|Fp|)  RMS(Fp)
   1     0.000     -76.10000     0.000     0.001     0.001
   2     0.500     -76.05000    31.374     0.005     0.003
   3     1.000     -75.95000    94.122     0.008     0.005
"""


class TestParseTrajectoryNeb(unittest.TestCase):
    def test_neb_image_count(self):
        p = _parse_method(_NEB_PATH_SUMMARY, "parse_trajectory")
        neb = [s for s in p.data["scan_steps"] if s.get("type") == "neb_image"]
        self.assertEqual(len(neb), 3)

    def test_neb_energies(self):
        p = _parse_method(_NEB_PATH_SUMMARY, "parse_trajectory")
        neb = [s for s in p.data["scan_steps"] if s.get("type") == "neb_image"]
        self.assertAlmostEqual(neb[0]["energy"], -76.10000, places=4)
        self.assertAlmostEqual(neb[1]["energy"], -76.05000, places=4)

    def test_neb_distances(self):
        p = _parse_method(_NEB_PATH_SUMMARY, "parse_trajectory")
        neb = [s for s in p.data["scan_steps"] if s.get("type") == "neb_image"]
        self.assertAlmostEqual(neb[0]["dist"], 0.000, places=3)
        self.assertAlmostEqual(neb[1]["dist"], 0.500, places=3)
        self.assertAlmostEqual(neb[2]["dist"], 1.000, places=3)


# ---------------------------------------------------------------------------
# TestParseOptCycleConvergence
# ---------------------------------------------------------------------------

# Verify convergence info dict is populated for opt cycles.

_OPT_CONVERGENCE = """\
                         **** OPTIMIZATION CYCLE    1 ****
  FINAL SINGLE POINT ENERGY    -76.23456789
  GEOMETRY CONVERGENCE
  ---
  RMS gradient      0.000123  0.000100  YES
  MAX gradient      0.000234  0.000200  YES
  RMS step          0.000345  0.000400  YES
  MAX step          0.000456  0.000600  YES
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.00000    0.00000    0.11779
  H      0.75545    0.00000   -0.47116
  H     -0.75545    0.00000   -0.47116

"""


class TestParseOptCycleConvergence(unittest.TestCase):
    def test_convergence_dict_populated(self):
        p = _parse_method(_OPT_CONVERGENCE, "parse_trajectory")
        step = p.data["scan_steps"][0]
        self.assertIsInstance(step["convergence"], dict)
        self.assertGreater(len(step["convergence"]), 0)

    def test_convergence_yes_flag(self):
        p = _parse_method(_OPT_CONVERGENCE, "parse_trajectory")
        conv = p.data["scan_steps"][0]["convergence"]
        # At least one key should have converged='YES'
        yes_count = sum(1 for v in conv.values() if v.get("converged") == "YES")
        self.assertGreater(yes_count, 0)


# ---------------------------------------------------------------------------
# TestParseGradientsMultipleBlocks
# ---------------------------------------------------------------------------

_TWO_GRADIENT_BLOCKS = """\
CARTESIAN GRADIENT
------------------
   1   O   :    0.00012   0.00000  -0.00034
   2   H   :    0.00003   0.00000   0.00017
-------

CARTESIAN GRADIENT
------------------
   1   O   :    0.00005   0.00000  -0.00012
   2   H   :    0.00001   0.00000   0.00006
-------
"""


class TestParseGradientsMultipleBlocks(unittest.TestCase):
    def test_two_blocks_stored_in_all_gradients(self):
        p = _parse_method(_TWO_GRADIENT_BLOCKS, "parse_gradients")
        self.assertEqual(len(p.data["all_gradients"]), 2)

    def test_default_gradients_is_last_block(self):
        p = _parse_method(_TWO_GRADIENT_BLOCKS, "parse_gradients")
        # The default "gradients" should be from the second block
        self.assertAlmostEqual(p.data["gradients"][0]["vector"][0], 0.00005, places=5)

    def test_first_block_accessible(self):
        p = _parse_method(_TWO_GRADIENT_BLOCKS, "parse_gradients")
        first = p.data["all_gradients"][0]["grads"]
        self.assertAlmostEqual(first[0]["vector"][0], 0.00012, places=5)


if __name__ == "__main__":
    unittest.main()

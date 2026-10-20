"""
tests/test_parser_samples.py
Integration tests for OrcaParser using real ORCA output files.

All tests load actual .out files from tests/sample_outputs/ and assert
concrete values extracted from those files.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from _parser_loader import load_standalone_parser  # noqa: E402

# ---------------------------------------------------------------------------
# Bootstrap: load parser without Qt
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(__file__)
_SAMPLES = os.path.join(_HERE, "sample_outputs")


def _load_parser():
    return load_standalone_parser("orca_parser_standalone")


_parser_mod = _load_parser()
OrcaParser = _parser_mod.OrcaParser


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load(filename):
    path = os.path.join(_SAMPLES, filename)
    p = OrcaParser()
    with open(path, encoding="utf-8", errors="replace") as f:
        content = f.read()
    p.load_from_memory(content, path)
    return p


def _opt_cycles(p):
    return [s for s in p.data["scan_steps"] if s.get("type") == "opt_cycle"]


def _opt_final(p):
    return [s for s in p.data["scan_steps"] if s.get("type") == "opt_final"]


# ---------------------------------------------------------------------------
# Benzene geometry optimization  —  ORCA 6
# ---------------------------------------------------------------------------


class TestBenzeneOptOrca6(unittest.TestCase):
    """benzene-opt.out  (ORCA 6.1.1, 2 cycles, converged)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt.out")

    # --- basic ---

    def test_version(self):
        self.assertEqual(self.p.data["version"], "6.1.1")

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_final_energy(self):
        self.assertAlmostEqual(self.p.data["scf_energy"], -231.930734972084, places=6)

    def test_charge(self):
        self.assertEqual(self.p.data["charge"], 0)

    def test_multiplicity(self):
        self.assertEqual(self.p.data["mult"], 1)

    # --- final geometry (from FINAL ENERGY EVALUATION) ---

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_atom_symbols(self):
        syms = self.p.data["atoms"]
        self.assertEqual(syms.count("C"), 6)
        self.assertEqual(syms.count("H"), 6)

    def test_final_coords_first_atom(self):
        """data['coords'][0] must be the FINAL ENERGY EVALUATION geometry, not cycle 1."""
        x, y, z = self.p.data["coords"][0]
        self.assertAlmostEqual(x, 0.805933, places=4)
        self.assertAlmostEqual(y, -1.143600, places=4)
        self.assertAlmostEqual(z, -0.008261, places=4)

    # --- trajectory ---

    def test_total_scan_steps(self):
        self.assertEqual(len(self.p.data["scan_steps"]), 3)

    def test_opt_cycle_count(self):
        self.assertEqual(len(_opt_cycles(self.p)), 2)

    def test_opt_final_count(self):
        self.assertEqual(len(_opt_final(self.p)), 1)

    def test_opt_final_is_last(self):
        self.assertEqual(self.p.data["scan_steps"][-1]["type"], "opt_final")

    def test_opt_final_coords_match_data_coords(self):
        """opt_final coords must equal data['coords'] (both from FINAL ENERGY EVALUATION)."""
        final_step = _opt_final(self.p)[0]
        self.assertEqual(final_step["atoms"], self.p.data["atoms"])
        for i, (fc, dc) in enumerate(zip(final_step["coords"], self.p.data["coords"])):
            for j in range(3):
                self.assertAlmostEqual(
                    fc[j], dc[j], places=5, msg=f"atom {i} coord {j} mismatch"
                )

    def test_opt_final_energy(self):
        final_step = _opt_final(self.p)[0]
        self.assertAlmostEqual(final_step["energy"], -231.930734972084, places=6)

    def test_cycle1_step_number(self):
        self.assertEqual(_opt_cycles(self.p)[0]["step"], 1)

    def test_cycle2_step_number(self):
        self.assertEqual(_opt_cycles(self.p)[1]["step"], 2)

    def test_cycle1_has_atoms(self):
        self.assertEqual(len(_opt_cycles(self.p)[0]["atoms"]), 12)

    def test_cycle1_coords_differ_from_final(self):
        """Cycle 1 starting geometry must differ from the FINAL ENERGY EVALUATION coords."""
        cycle1_x = _opt_cycles(self.p)[0]["coords"][0][0]
        final_x = self.p.data["coords"][0][0]
        self.assertNotAlmostEqual(cycle1_x, final_x, places=3)

    # --- dipole ---

    def test_dipole_present(self):
        self.assertIsNotNone(self.p.data["dipole"])

    def test_dipole_near_zero(self):
        """Benzene is symmetric — dipole magnitude must be near zero."""
        mag = self.p.data["dipole"]["magnitude"]
        self.assertLess(abs(mag), 0.1)

    # --- charges ---

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)

    def test_mulliken_sum_near_zero(self):
        total = sum(e["charge"] for e in self.p.data["charges"]["Mulliken"])
        self.assertAlmostEqual(total, 0.0, places=2)


# ---------------------------------------------------------------------------
# Benzene geometry optimization  —  ORCA 5
# ---------------------------------------------------------------------------


class TestBenzeneOptOrca5(unittest.TestCase):
    """benzene-opt_5.out  (ORCA 5.0.4, 2 cycles, converged)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt_5.out")

    def test_version(self):
        self.assertEqual(self.p.data["version"], "5.0.4")

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_final_energy(self):
        self.assertAlmostEqual(self.p.data["scf_energy"], -231.930718933873, places=6)

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_final_coords_first_atom(self):
        x, y, z = self.p.data["coords"][0]
        self.assertAlmostEqual(x, 0.805991, places=4)
        self.assertAlmostEqual(y, -1.143635, places=4)
        self.assertAlmostEqual(z, -0.008264, places=4)

    def test_total_scan_steps(self):
        self.assertEqual(len(self.p.data["scan_steps"]), 3)

    def test_opt_cycle_count(self):
        self.assertEqual(len(_opt_cycles(self.p)), 2)

    def test_opt_final_count(self):
        self.assertEqual(len(_opt_final(self.p)), 1)

    def test_opt_final_is_last(self):
        self.assertEqual(self.p.data["scan_steps"][-1]["type"], "opt_final")

    def test_opt_final_coords_match_data_coords(self):
        final_step = _opt_final(self.p)[0]
        for i, (fc, dc) in enumerate(zip(final_step["coords"], self.p.data["coords"])):
            for j in range(3):
                self.assertAlmostEqual(fc[j], dc[j], places=5)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)


# ---------------------------------------------------------------------------
# Acetone geometry optimization  —  ORCA 6
# ---------------------------------------------------------------------------


class TestAcetoneOptOrca6(unittest.TestCase):
    """acetone-opt.out  (ORCA 6.x, 5 cycles, converged)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("acetone-opt.out")

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_atom_count(self):
        # C3H6O = 10 atoms
        self.assertEqual(len(self.p.data["atoms"]), 10)

    def test_atom_symbols(self):
        syms = self.p.data["atoms"]
        self.assertEqual(syms.count("C"), 3)
        self.assertEqual(syms.count("O"), 1)
        self.assertEqual(syms.count("H"), 6)

    def test_total_scan_steps(self):
        # 5 opt_cycles + 1 opt_final
        self.assertEqual(len(self.p.data["scan_steps"]), 6)

    def test_opt_cycle_count(self):
        self.assertEqual(len(_opt_cycles(self.p)), 5)

    def test_opt_final_count(self):
        self.assertEqual(len(_opt_final(self.p)), 1)

    def test_opt_final_is_last(self):
        self.assertEqual(self.p.data["scan_steps"][-1]["type"], "opt_final")

    def test_opt_final_coords_match_data_coords(self):
        final_step = _opt_final(self.p)[0]
        self.assertEqual(final_step["atoms"], self.p.data["atoms"])
        for i, (fc, dc) in enumerate(zip(final_step["coords"], self.p.data["coords"])):
            for j in range(3):
                self.assertAlmostEqual(fc[j], dc[j], places=5)

    def test_final_coords_first_atom(self):
        x, y, z = self.p.data["coords"][0]
        self.assertAlmostEqual(x, -1.274842, places=4)
        self.assertAlmostEqual(y, 0.247667, places=4)
        self.assertAlmostEqual(z, 0.079907, places=4)

    def test_cycle1_coords_differ_from_final(self):
        cycle1_coords = _opt_cycles(self.p)[0]["coords"][0]
        final_coords = self.p.data["coords"][0]
        diffs = [abs(cycle1_coords[j] - final_coords[j]) for j in range(3)]
        self.assertGreater(max(diffs), 0.001)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 10)

    def test_mulliken_sum_near_zero(self):
        total = sum(e["charge"] for e in self.p.data["charges"]["Mulliken"])
        self.assertAlmostEqual(total, 0.0, places=2)

    def test_dipole_present(self):
        self.assertIsNotNone(self.p.data["dipole"])

    def test_dipole_nonzero(self):
        """Acetone is polar — dipole must be significant."""
        mag = self.p.data["dipole"]["magnitude"]
        self.assertGreater(mag, 2.0)

    def test_dipole_magnitude_debye(self):
        """Acetone dipole ~2.787 Debye (parser reads Magnitude (Debye) line)."""
        self.assertAlmostEqual(self.p.data["dipole"]["magnitude"], 2.787, places=2)

    def test_dipole_vector_z_dominant(self):
        """C=O axis is roughly along z — z component should be largest."""
        vec = self.p.data["dipole"]["vector"]
        self.assertGreater(abs(vec[2]), abs(vec[0]))


# ---------------------------------------------------------------------------
# Acetone geometry optimization  —  ORCA 5
# ---------------------------------------------------------------------------


class TestAcetoneOptOrca5(unittest.TestCase):
    """acetone-opt_5.out  (ORCA 5.x, 5 cycles, converged)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("acetone-opt_5.out")

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 10)

    def test_total_scan_steps(self):
        self.assertEqual(len(self.p.data["scan_steps"]), 6)

    def test_opt_cycle_count(self):
        self.assertEqual(len(_opt_cycles(self.p)), 5)

    def test_opt_final_count(self):
        self.assertEqual(len(_opt_final(self.p)), 1)

    def test_opt_final_is_last(self):
        self.assertEqual(self.p.data["scan_steps"][-1]["type"], "opt_final")

    def test_opt_final_coords_match_data_coords(self):
        final_step = _opt_final(self.p)[0]
        for i, (fc, dc) in enumerate(zip(final_step["coords"], self.p.data["coords"])):
            for j in range(3):
                self.assertAlmostEqual(fc[j], dc[j], places=5)

    def test_dipole_present(self):
        self.assertIsNotNone(self.p.data["dipole"])

    def test_dipole_magnitude_debye(self):
        """Acetone dipole ~2.787 Debye."""
        self.assertAlmostEqual(self.p.data["dipole"]["magnitude"], 2.787, places=2)


# ---------------------------------------------------------------------------
# Benzene single-point energy  —  ORCA 6
# ---------------------------------------------------------------------------


class TestBenzeneEneOrca6(unittest.TestCase):
    """benzene-opt-ene.out  (ORCA 6.1.1, SP, no optimization)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-ene.out")

    def test_final_energy(self):
        self.assertAlmostEqual(self.p.data["scf_energy"], -231.930784926115, places=6)

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)

    def test_no_opt_final(self):
        self.assertEqual(len(_opt_final(self.p)), 0)

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_dipole_near_zero(self):
        self.assertIsNotNone(self.p.data["dipole"])
        self.assertLess(abs(self.p.data["dipole"]["magnitude"]), 0.1)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)


# ---------------------------------------------------------------------------
# Benzene single-point energy  —  ORCA 5
# ---------------------------------------------------------------------------


class TestBenzeneEneOrca5(unittest.TestCase):
    """benzene-opt-ene_5.out  (ORCA 5.0.4, SP, no optimization)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-ene_5.out")

    def test_final_energy(self):
        self.assertAlmostEqual(self.p.data["scf_energy"], -231.930718955827, places=6)

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)

    def test_no_opt_final(self):
        self.assertEqual(len(_opt_final(self.p)), 0)

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)


# ---------------------------------------------------------------------------
# Benzene NMR  —  ORCA 6
# ---------------------------------------------------------------------------


class TestBenzeneNmrOrca6(unittest.TestCase):
    """benzene-opt-nmr.out  (ORCA 6.1.1, NMR shielding)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-nmr.out")

    def test_final_energy(self):
        self.assertAlmostEqual(self.p.data["scf_energy"], -231.930784926115, places=6)

    def test_nmr_shielding_present(self):
        self.assertGreater(len(self.p.data["nmr_shielding"]), 0)

    def test_nmr_shielding_count(self):
        # 6 C + 6 H = 12 nuclei
        self.assertEqual(len(self.p.data["nmr_shielding"]), 12)

    def test_nmr_carbon_shielding_range(self):
        """Benzene carbon shielding should be roughly 40-80 ppm."""
        carbon_shields = [
            e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "C"
        ]
        self.assertEqual(len(carbon_shields), 6)
        for s in carbon_shields:
            self.assertGreater(s, 0.0)

    def test_nmr_hydrogen_shielding_range(self):
        """Benzene H shielding should be positive."""
        h_shields = [
            e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "H"
        ]
        self.assertEqual(len(h_shields), 6)
        for s in h_shields:
            self.assertGreater(s, 0.0)

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)


# ---------------------------------------------------------------------------
# Benzene NMR  —  ORCA 5
# ---------------------------------------------------------------------------


class TestBenzeneNmrOrca5(unittest.TestCase):
    """benzene-opt-nmr_5.out  (ORCA 5.0.4, NMR shielding)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-nmr_5.out")

    def test_final_energy(self):
        self.assertAlmostEqual(self.p.data["scf_energy"], -231.930718955827, places=6)

    def test_nmr_shielding_present(self):
        self.assertGreater(len(self.p.data["nmr_shielding"]), 0)

    def test_nmr_shielding_count(self):
        self.assertEqual(len(self.p.data["nmr_shielding"]), 12)

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)


# ---------------------------------------------------------------------------
# Benzene TDDFT excited states  —  ORCA 6
# ---------------------------------------------------------------------------


class TestBenzeneVexOrca6(unittest.TestCase):
    """benzene-opt-vex.out  (ORCA 6.1.1, TD-DFT, 5 states)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-vex.out")

    def test_final_energy(self):
        self.assertAlmostEqual(self.p.data["scf_energy"], -231.728888949460, places=6)

    def test_tddft_state_count(self):
        self.assertEqual(len(self.p.data["tddft"]), 5)

    def test_tddft_first_state_ev(self):
        self.assertAlmostEqual(self.p.data["tddft"][0]["energy_ev"], 5.494, places=2)

    def test_tddft_states_sorted_ascending(self):
        evs = [s["energy_ev"] for s in self.p.data["tddft"]]
        self.assertEqual(evs, sorted(evs))

    def test_tddft_all_states_have_transitions(self):
        for state in self.p.data["tddft"]:
            self.assertGreater(len(state["transitions"]), 0)

    def test_tddft_energy_range(self):
        """All 5 states should be between 5 and 8 eV for benzene."""
        for state in self.p.data["tddft"]:
            self.assertGreater(state["energy_ev"], 5.0)
            self.assertLess(state["energy_ev"], 8.5)

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)


# ---------------------------------------------------------------------------
# Benzene TDDFT excited states  —  ORCA 5
# ---------------------------------------------------------------------------


class TestBenzeneVexOrca5(unittest.TestCase):
    """benzene-opt-vex_5.out  (ORCA 5.0.4, TD-DFT, 5 states)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-vex_5.out")

    def test_final_energy(self):
        self.assertAlmostEqual(self.p.data["scf_energy"], -231.728828351432, places=6)

    def test_tddft_state_count(self):
        self.assertEqual(len(self.p.data["tddft"]), 5)

    def test_tddft_first_state_ev(self):
        self.assertAlmostEqual(self.p.data["tddft"][0]["energy_ev"], 5.494, places=2)

    def test_tddft_states_sorted_ascending(self):
        evs = [s["energy_ev"] for s in self.p.data["tddft"]]
        self.assertEqual(evs, sorted(evs))

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)


# ---------------------------------------------------------------------------
# Cross-version consistency
# ---------------------------------------------------------------------------


class TestVersionConsistency(unittest.TestCase):
    """Verify that ORCA 5 and 6 outputs yield structurally identical results."""

    def test_benzene_opt_same_atom_count(self):
        p6 = _load("benzene-opt.out")
        p5 = _load("benzene-opt_5.out")
        self.assertEqual(len(p6.data["atoms"]), len(p5.data["atoms"]))

    def test_benzene_opt_same_step_count(self):
        p6 = _load("benzene-opt.out")
        p5 = _load("benzene-opt_5.out")
        self.assertEqual(len(p6.data["scan_steps"]), len(p5.data["scan_steps"]))

    def test_benzene_opt_same_step_types(self):
        p6 = _load("benzene-opt.out")
        p5 = _load("benzene-opt_5.out")
        types6 = [s["type"] for s in p6.data["scan_steps"]]
        types5 = [s["type"] for s in p5.data["scan_steps"]]
        self.assertEqual(types6, types5)

    def test_benzene_ene_same_atom_count(self):
        p6 = _load("benzene-opt-ene.out")
        p5 = _load("benzene-opt-ene_5.out")
        self.assertEqual(len(p6.data["atoms"]), len(p5.data["atoms"]))

    def test_benzene_nmr_same_shielding_count(self):
        p6 = _load("benzene-opt-nmr.out")
        p5 = _load("benzene-opt-nmr_5.out")
        self.assertEqual(len(p6.data["nmr_shielding"]), len(p5.data["nmr_shielding"]))

    def test_benzene_vex_same_state_count(self):
        p6 = _load("benzene-opt-vex.out")
        p5 = _load("benzene-opt-vex_5.out")
        self.assertEqual(len(p6.data["tddft"]), len(p5.data["tddft"]))

    def test_acetone_opt_same_step_count(self):
        p6 = _load("acetone-opt.out")
        p5 = _load("acetone-opt_5.out")
        self.assertEqual(len(p6.data["scan_steps"]), len(p5.data["scan_steps"]))

    def test_acetone_opt_same_atom_count(self):
        p6 = _load("acetone-opt.out")
        p5 = _load("acetone-opt_5.out")
        self.assertEqual(len(p6.data["atoms"]), len(p5.data["atoms"]))


# ---------------------------------------------------------------------------
# Basis Set  —  ORCA 6 (benzene-opt-ene.out has "BASIS SET IN INPUT FORMAT")
# ---------------------------------------------------------------------------


class TestBasisSetOrca6(unittest.TestCase):
    """benzene-opt-ene.out — parse_basis_set populates basis_set_shells."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-ene.out")

    def test_basis_set_shells_present(self):
        shells = self.p.data.get("basis_set_shells", [])
        self.assertGreater(len(shells), 0, "basis_set_shells must not be empty")

    def test_basis_set_shells_have_required_keys(self):
        for sh in self.p.data["basis_set_shells"]:
            self.assertIn("atom_idx", sh)
            self.assertIn("l", sh)
            self.assertIn("exps", sh)
            self.assertIn("coeffs", sh)

    def test_basis_set_l_values_valid(self):
        """l=0 (S), 1 (P), 2 (D) etc. — all must be non-negative integers."""
        for sh in self.p.data["basis_set_shells"]:
            self.assertGreaterEqual(sh["l"], 0)

    def test_basis_set_exponents_positive(self):
        for sh in self.p.data["basis_set_shells"]:
            for exp in sh["exps"]:
                self.assertGreater(exp, 0.0)

    def test_basis_set_atom_indices_in_range(self):
        n_atoms = len(self.p.data["atoms"])
        for sh in self.p.data["basis_set_shells"]:
            self.assertGreaterEqual(sh["atom_idx"], 0)
            self.assertLess(sh["atom_idx"], n_atoms)

    def test_basis_set_shell_count_reasonable(self):
        """Benzene (12 atoms) with def2-TZVP should produce >30 shells."""
        self.assertGreater(len(self.p.data["basis_set_shells"]), 30)


class TestBasisSetOrca5(unittest.TestCase):
    """benzene-opt-ene_5.out — basis set parsing on ORCA 5 output."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-ene_5.out")

    def test_basis_set_shells_present(self):
        shells = self.p.data.get("basis_set_shells", [])
        self.assertGreater(len(shells), 0)

    def test_basis_set_shells_have_required_keys(self):
        for sh in self.p.data["basis_set_shells"]:
            self.assertIn("atom_idx", sh)
            self.assertIn("l", sh)
            self.assertIn("exps", sh)
            self.assertIn("coeffs", sh)

    def test_basis_set_l_values_valid(self):
        for sh in self.p.data["basis_set_shells"]:
            self.assertGreaterEqual(sh["l"], 0)

    def test_basis_set_exponents_positive(self):
        for sh in self.p.data["basis_set_shells"]:
            for exp in sh["exps"]:
                self.assertGreater(exp, 0.0)

    def test_basis_set_atom_indices_in_range(self):
        n_atoms = len(self.p.data["atoms"])
        for sh in self.p.data["basis_set_shells"]:
            self.assertGreaterEqual(sh["atom_idx"], 0)
            self.assertLess(sh["atom_idx"], n_atoms)

    def test_basis_set_same_count_as_orca6(self):
        p6 = _load("benzene-opt-ene.out")
        self.assertEqual(
            len(self.p.data["basis_set_shells"]),
            len(p6.data["basis_set_shells"]),
        )


# ---------------------------------------------------------------------------
# Orbital Energies  —  ORCA 6 (benzene-opt-ene.out has ORBITAL ENERGIES)
# ---------------------------------------------------------------------------


class TestOrbitalEnergiesOrca6(unittest.TestCase):
    """benzene-opt-ene.out — parse_orbital_energies."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-ene.out")

    def test_orbital_energies_present(self):
        self.assertGreater(len(self.p.data["orbital_energies"]), 0)

    def test_orbital_energies_have_required_keys(self):
        for orb in self.p.data["orbital_energies"]:
            for key in (
                "index",
                "occupation",
                "energy_eh",
                "energy_ev",
                "spin",
                "type",
            ):
                self.assertIn(key, orb)

    def test_occupied_orbitals_present(self):
        occ = [o for o in self.p.data["orbital_energies"] if o["type"] == "occupied"]
        self.assertGreater(len(occ), 0)

    def test_virtual_orbitals_present(self):
        virt = [o for o in self.p.data["orbital_energies"] if o["type"] == "virtual"]
        self.assertGreater(len(virt), 0)

    def test_homo_energy_negative(self):
        """HOMO energy must be negative (bound state)."""
        occ = [o for o in self.p.data["orbital_energies"] if o["type"] == "occupied"]
        homo = max(occ, key=lambda o: o["energy_eh"])
        self.assertLess(homo["energy_eh"], 0.0)

    def test_lumo_energy_above_homo(self):
        occ = [o for o in self.p.data["orbital_energies"] if o["type"] == "occupied"]
        virt = [o for o in self.p.data["orbital_energies"] if o["type"] == "virtual"]
        homo_e = max(o["energy_eh"] for o in occ)
        lumo_e = min(o["energy_eh"] for o in virt)
        self.assertGreater(lumo_e, homo_e)

    def test_spin_is_restricted(self):
        spins = {o["spin"] for o in self.p.data["orbital_energies"]}
        self.assertIn("restricted", spins)

    def test_mos_backward_compat(self):
        """data['mos'] must be identical to data['orbital_energies']."""
        self.assertEqual(self.p.data["mos"], self.p.data["orbital_energies"])


# ---------------------------------------------------------------------------
# Orbital Energies  —  ORCA 5
# ---------------------------------------------------------------------------


class TestOrbitalEnergiesOrca5(unittest.TestCase):
    """benzene-opt-ene_5.out — parse_orbital_energies on ORCA 5."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-ene_5.out")

    def test_orbital_energies_present(self):
        self.assertGreater(len(self.p.data["orbital_energies"]), 0)

    def test_orbital_energies_have_required_keys(self):
        for orb in self.p.data["orbital_energies"]:
            for key in (
                "index",
                "occupation",
                "energy_eh",
                "energy_ev",
                "spin",
                "type",
            ):
                self.assertIn(key, orb)

    def test_occupied_orbitals_present(self):
        occ = [o for o in self.p.data["orbital_energies"] if o["type"] == "occupied"]
        self.assertGreater(len(occ), 0)

    def test_virtual_orbitals_present(self):
        virt = [o for o in self.p.data["orbital_energies"] if o["type"] == "virtual"]
        self.assertGreater(len(virt), 0)

    def test_homo_energy_negative(self):
        occ = [o for o in self.p.data["orbital_energies"] if o["type"] == "occupied"]
        homo = max(occ, key=lambda o: o["energy_eh"])
        self.assertLess(homo["energy_eh"], 0.0)

    def test_lumo_energy_above_homo(self):
        occ = [o for o in self.p.data["orbital_energies"] if o["type"] == "occupied"]
        virt = [o for o in self.p.data["orbital_energies"] if o["type"] == "virtual"]
        homo_e = max(o["energy_eh"] for o in occ)
        lumo_e = min(o["energy_eh"] for o in virt)
        self.assertGreater(lumo_e, homo_e)

    def test_spin_is_restricted(self):
        spins = {o["spin"] for o in self.p.data["orbital_energies"]}
        self.assertIn("restricted", spins)

    def test_orbital_count_reasonable(self):
        """ORCA 5 ene file may use a larger basis — just assert a sensible count."""
        n = len(self.p.data["orbital_energies"])
        self.assertGreater(n, 20)  # at least occupied + some virtuals

    def test_mos_backward_compat(self):
        self.assertEqual(self.p.data["mos"], self.p.data["orbital_energies"])


# ---------------------------------------------------------------------------
# SCF Trace  —  ORCA 6 (benzene-opt-ene.out: 1 SCF block; benzene-opt.out: 2+)
# ---------------------------------------------------------------------------


class TestScfTraceOrca6(unittest.TestCase):
    """parse_scf_trace — iteration energy traces."""

    @classmethod
    def setUpClass(cls):
        cls.p_ene = _load("benzene-opt-ene.out")
        cls.p_opt = _load("benzene-opt.out")

    def test_scf_trace_present_ene(self):
        self.assertGreater(len(self.p_ene.data["scf_traces"]), 0)

    def test_scf_trace_has_required_keys(self):
        trace = self.p_ene.data["scf_traces"][0]
        self.assertIn("step", trace)
        self.assertIn("iterations", trace)

    def test_scf_iterations_nonempty(self):
        for trace in self.p_ene.data["scf_traces"]:
            self.assertGreater(len(trace["iterations"]), 0)

    def test_scf_each_iteration_has_iter_and_energy(self):
        for trace in self.p_ene.data["scf_traces"]:
            for it in trace["iterations"]:
                self.assertIn("iter", it)
                self.assertIn("energy", it)

    def test_scf_energies_are_negative(self):
        """SCF energies for benzene must be large and negative."""
        for trace in self.p_ene.data["scf_traces"]:
            for it in trace["iterations"]:
                self.assertLess(it["energy"], 0.0)

    def test_scf_opt_has_multiple_cycle_traces(self):
        """An optimization run should produce one SCF trace per cycle."""
        self.assertGreaterEqual(len(self.p_opt.data["scf_traces"]), 2)

    def test_scf_opt_cycle_labels(self):
        labels = [t["step"] for t in self.p_opt.data["scf_traces"]]
        cycle_labels = [
            lbl
            for lbl in labels
            if "Cycle" in lbl
            or "Initial" in lbl
            or "Final" in lbl
            or "Post" in lbl
            or "Property" in lbl
        ]
        self.assertGreater(len(cycle_labels), 0)


# ---------------------------------------------------------------------------
# SCF Trace  —  ORCA 5
# ---------------------------------------------------------------------------


class TestScfTraceOrca5(unittest.TestCase):
    """ORCA 5 SCF trace: benzene-opt-ene_5.out (2 blocks) + benzene-opt_5.out (4 blocks)."""

    @classmethod
    def setUpClass(cls):
        cls.p_ene = _load("benzene-opt-ene_5.out")
        cls.p_opt = _load("benzene-opt_5.out")

    def test_scf_trace_present_ene(self):
        self.assertGreaterEqual(len(self.p_ene.data["scf_traces"]), 1)

    def test_scf_trace_has_required_keys(self):
        trace = self.p_ene.data["scf_traces"][0]
        self.assertIn("step", trace)
        self.assertIn("iterations", trace)

    def test_scf_iterations_nonempty(self):
        for trace in self.p_ene.data["scf_traces"]:
            self.assertGreater(len(trace["iterations"]), 0)

    def test_scf_each_iteration_has_iter_and_energy(self):
        for trace in self.p_ene.data["scf_traces"]:
            for it in trace["iterations"]:
                self.assertIn("iter", it)
                self.assertIn("energy", it)

    def test_scf_energies_are_negative(self):
        for trace in self.p_ene.data["scf_traces"]:
            for it in trace["iterations"]:
                self.assertLess(it["energy"], 0.0)

    def test_scf_opt_has_multiple_cycle_traces(self):
        """benzene-opt_5.out has 4 SCF blocks (2 opt cycles + 1 final eval + aux)."""
        self.assertGreaterEqual(len(self.p_opt.data["scf_traces"]), 2)

    def test_scf_opt_cycle_labels(self):
        labels = [t["step"] for t in self.p_opt.data["scf_traces"]]
        cycle_labels = [
            lbl
            for lbl in labels
            if "Cycle" in lbl
            or "Initial" in lbl
            or "Final" in lbl
            or "Post" in lbl
            or "Property" in lbl
        ]
        self.assertGreater(len(cycle_labels), 0)


# ---------------------------------------------------------------------------
# TDDFT detail  —  ORCA 6 (benzene-opt-vex.out: absorption spectrum tables)
# ---------------------------------------------------------------------------


class TestTddftDetailOrca6(unittest.TestCase):
    """benzene-opt-vex.out — deeper TDDFT field coverage."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-vex.out")

    def test_tddft_state_has_osc_len(self):
        """osc_len (electric dipole length gauge) must be present."""
        for state in self.p.data["tddft"]:
            self.assertIn("osc_len", state)

    def test_tddft_state_has_osc_vel(self):
        """osc_vel (velocity gauge) must be present."""
        for state in self.p.data["tddft"]:
            self.assertIn("osc_vel", state)

    def test_tddft_state_has_energy_nm(self):
        for state in self.p.data["tddft"]:
            self.assertIn("energy_nm", state)
            self.assertGreater(state["energy_nm"], 0.0)

    def test_tddft_state_has_energy_cm(self):
        for state in self.p.data["tddft"]:
            self.assertIn("energy_cm", state)
            self.assertGreater(state["energy_cm"], 0.0)

    def test_tddft_energy_ev_consistent_with_nm(self):
        """eV * nm ≈ 1239.8 (hc in eV·nm)."""
        for state in self.p.data["tddft"]:
            product = state["energy_ev"] * state["energy_nm"]
            self.assertAlmostEqual(product, 1239.84, delta=5.0)

    def test_tddft_osc_non_negative(self):
        for state in self.p.data["tddft"]:
            self.assertGreaterEqual(state["osc"], 0.0)

    def test_tddft_transitions_strings(self):
        for state in self.p.data["tddft"]:
            for t in state.get("transitions", []):
                self.assertIsInstance(t, str)

    def test_tddft_cd_fields_present(self):
        """CD spectrum fields must be present (benzene-opt-vex.out has CD SPECTRUM blocks)."""
        for state in self.p.data["tddft"]:
            self.assertIn("rot_len", state)
            self.assertIn("rot_vel", state)
            self.assertIn("rotatory_strength", state)

    def test_tddft_cd_near_zero_for_benzene(self):
        """Benzene is achiral — all rotatory strengths must be ~0."""
        for state in self.p.data["tddft"]:
            self.assertAlmostEqual(state["rotatory_strength"], 0.0, places=3)


class TestTddftDetailOrca5(unittest.TestCase):
    """benzene-opt-vex_5.out — ORCA 5 TDDFT detail."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-vex_5.out")

    def test_tddft_state_count(self):
        self.assertEqual(len(self.p.data["tddft"]), 5)

    def test_tddft_state_has_energy_ev(self):
        for state in self.p.data["tddft"]:
            self.assertGreater(state["energy_ev"], 0.0)

    def test_tddft_state_has_energy_nm(self):
        for state in self.p.data["tddft"]:
            self.assertIn("energy_nm", state)
            self.assertGreater(state["energy_nm"], 0.0)

    def test_tddft_state_has_energy_cm(self):
        for state in self.p.data["tddft"]:
            self.assertIn("energy_cm", state)
            self.assertGreater(state["energy_cm"], 0.0)

    def test_tddft_state_has_osc_len(self):
        for state in self.p.data["tddft"]:
            self.assertIn("osc_len", state)

    def test_tddft_state_has_osc_vel(self):
        for state in self.p.data["tddft"]:
            self.assertIn("osc_vel", state)

    def test_tddft_osc_non_negative(self):
        for state in self.p.data["tddft"]:
            self.assertGreaterEqual(state["osc"], 0.0)

    def test_tddft_transitions_present(self):
        for state in self.p.data["tddft"]:
            self.assertGreater(len(state.get("transitions", [])), 0)

    def test_tddft_energy_ev_consistent_with_nm(self):
        for state in self.p.data["tddft"]:
            product = state["energy_ev"] * state["energy_nm"]
            self.assertAlmostEqual(product, 1239.84, delta=5.0)

    def test_tddft_first_state_ev_consistent_with_orca6(self):
        """First state eV should match ORCA 6 within 0.05 eV."""
        p6 = _load("benzene-opt-vex.out")
        ev5 = self.p.data["tddft"][0]["energy_ev"]
        ev6 = p6.data["tddft"][0]["energy_ev"]
        self.assertAlmostEqual(ev5, ev6, delta=0.05)

    def test_tddft_cd_fields_present(self):
        """CD spectrum fields must be present (benzene-opt-vex_5.out has CD SPECTRUM block)."""
        for state in self.p.data["tddft"]:
            self.assertIn("rot_len", state)
            self.assertIn("rot_vel", state)
            self.assertIn("rotatory_strength", state)

    def test_tddft_cd_near_zero_for_benzene(self):
        """Benzene is achiral — all rotatory strengths must be ~0."""
        for state in self.p.data["tddft"]:
            self.assertAlmostEqual(state["rotatory_strength"], 0.0, places=3)


# ---------------------------------------------------------------------------
# NMR detail  —  ORCA 6 (benzene-opt-nmr.out has Mayer bond orders)
# ---------------------------------------------------------------------------


class TestNmrDetailOrca6(unittest.TestCase):
    """benzene-opt-nmr.out — NMR shielding fields + Mayer charges."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-nmr.out")

    def test_nmr_shielding_has_required_keys(self):
        for entry in self.p.data["nmr_shielding"]:
            for key in ("atom_idx", "atom_sym", "shielding"):
                self.assertIn(key, entry)

    def test_nmr_shielding_atom_indices_sequential(self):
        idxs = [e["atom_idx"] for e in self.p.data["nmr_shielding"]]
        self.assertEqual(idxs, sorted(idxs))

    def test_nmr_shielding_carbon_value_range(self):
        """Benzene C shielding ~40-80 ppm (absolute)."""
        c_vals = [
            e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "C"
        ]
        for v in c_vals:
            self.assertGreater(v, 10.0)
            self.assertLess(v, 200.0)

    def test_nmr_shielding_hydrogen_value_range(self):
        """Benzene H shielding ~20-35 ppm (absolute)."""
        h_vals = [
            e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "H"
        ]
        for v in h_vals:
            self.assertGreater(v, 10.0)
            self.assertLess(v, 50.0)

    def test_mayer_charges_present(self):
        """Mayer population analysis should be present in nmr output."""
        self.assertIn("Mayer", self.p.data["charges"])

    def test_mayer_charges_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mayer"]), 12)

    def test_mayer_charges_sum_near_zero(self):
        total = sum(e["charge"] for e in self.p.data["charges"]["Mayer"])
        self.assertAlmostEqual(total, 0.0, places=1)


class TestNmrDetailOrca5(unittest.TestCase):
    """benzene-opt-nmr_5.out — ORCA 5 NMR detail."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-nmr_5.out")

    def test_nmr_shielding_count(self):
        self.assertEqual(len(self.p.data["nmr_shielding"]), 12)

    def test_nmr_shielding_has_required_keys(self):
        for entry in self.p.data["nmr_shielding"]:
            for key in ("atom_idx", "atom_sym", "shielding"):
                self.assertIn(key, entry)

    def test_nmr_shielding_atom_indices_sequential(self):
        idxs = [e["atom_idx"] for e in self.p.data["nmr_shielding"]]
        self.assertEqual(idxs, sorted(idxs))

    def test_nmr_carbon_shielding_range(self):
        c_vals = [
            e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "C"
        ]
        self.assertEqual(len(c_vals), 6)
        for v in c_vals:
            self.assertGreater(v, 10.0)
            self.assertLess(v, 200.0)

    def test_nmr_hydrogen_shielding_range(self):
        h_vals = [
            e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "H"
        ]
        self.assertEqual(len(h_vals), 6)
        for v in h_vals:
            self.assertGreater(v, 10.0)
            self.assertLess(v, 50.0)

    def test_nmr_shielding_consistent_with_orca6(self):
        """Mean C shielding diff between ORCA 5 and 6 should be < 5 ppm."""
        p6 = _load("benzene-opt-nmr.out")
        c5 = [
            e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "C"
        ]
        c6 = [e["shielding"] for e in p6.data["nmr_shielding"] if e["atom_sym"] == "C"]
        mean_diff = abs(sum(c5) / len(c5) - sum(c6) / len(c6))
        self.assertLess(mean_diff, 5.0)

    def test_mayer_charges_present(self):
        self.assertIn("Mayer", self.p.data["charges"])

    def test_mayer_charges_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mayer"]), 12)

    def test_mayer_charges_sum_near_zero(self):
        total = sum(e["charge"] for e in self.p.data["charges"]["Mayer"])
        self.assertAlmostEqual(total, 0.0, places=1)


# ---------------------------------------------------------------------------
# Relaxed Surface Scan  —  ORCA 6 (ethane-scan.out)
# ---------------------------------------------------------------------------


class TestRelaxedSurfaceScan(unittest.TestCase):
    """ethane-scan.out — ORCA Relaxed Surface Scan."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("ethane-scan.out")

    def test_scan_step_count(self):
        """Should have 5 scan_steps (the finalized steps for the surface)."""
        scan_steps = [
            s for s in self.p.data.get("scan_steps", []) if s["type"] == "scan_step"
        ]
        self.assertEqual(len(scan_steps), 5)

    def test_total_trajectory_steps(self):
        """Should have 41 total steps in the trajectory (5 scan steps + 31 intermediate cycles + 5 final evals)."""
        self.assertEqual(len(self.p.data.get("scan_steps", [])), 41)

    def test_scan_step_coords_and_energies(self):
        scan_steps = [
            s for s in self.p.data.get("scan_steps", []) if s["type"] == "scan_step"
        ]
        self.assertEqual(scan_steps[0]["scan_coord"], -60.0)
        self.assertAlmostEqual(scan_steps[0]["energy"], -79.791846796153)
        self.assertEqual(scan_steps[-1]["scan_coord"], 60.0)
        self.assertAlmostEqual(scan_steps[-1]["energy"], -79.791847846936)

        for step in scan_steps:
            self.assertEqual(len(step["atoms"]), 8)
            self.assertEqual(len(step["coords"]), 8)
            self.assertIn("scan_coord", step)
            self.assertIn("energy", step)

    def test_opt_final_count(self):
        """Should have exactly 5 opt_final steps (one for each scan step)."""
        opt_finals = [
            s for s in self.p.data.get("scan_steps", []) if s["type"] == "opt_final"
        ]
        self.assertEqual(len(opt_finals), 5)

    def test_opt_final_scan_step_ids(self):
        """Each opt_final must be correctly assigned to its corresponding scan_step_id."""
        opt_finals = [
            s for s in self.p.data.get("scan_steps", []) if s["type"] == "opt_final"
        ]
        ids = [s["scan_step_id"] for s in opt_finals]
        self.assertEqual(ids, [1, 2, 3, 4, 5])

    def test_opt_final_energies(self):
        """Each opt_final must have the correct energy matching the scan step final energy."""
        opt_finals = [
            s for s in self.p.data.get("scan_steps", []) if s["type"] == "opt_final"
        ]
        expected_energies = [
            -79.791846796153,
            -79.789747119085,
            -79.787432724318,
            -79.789747334378,
            -79.791847846936,
        ]
        for opt_f, expected in zip(opt_finals, expected_energies):
            self.assertAlmostEqual(opt_f["energy"], expected, places=8)


# ---------------------------------------------------------------------------
# Open-shell / UHF orbital energies  —  ORCA 6 (o2-opt.out, triplet O2)
# ---------------------------------------------------------------------------


class TestO2OpenShellOrbitals(unittest.TestCase):
    """o2-opt_orb.out (ORCA 6, O2 triplet, UHF) — alpha/beta orbital separation.

    Regression for the open-shell orbital bug: real ORCA prints the
    "NO OCC E(Eh) E(eV)" column header directly under "SPIN UP/DOWN ORBITALS"
    (no "---" separator). The parser used to read only a single "restricted"
    block, mislabeling the alpha orbitals and dropping the entire beta manifold.
    """

    @classmethod
    def setUpClass(cls):
        cls.p = _load("o2-opt_orb.out")

    def test_is_open_shell(self):
        self.assertEqual(self.p.data["mult"], 3)

    def test_both_spins_present(self):
        spins = {o["spin"] for o in self.p.data["orbital_energies"]}
        self.assertIn("alpha", spins)
        self.assertIn("beta", spins)

    def test_no_restricted_label_for_uhf(self):
        spins = {o["spin"] for o in self.p.data["orbital_energies"]}
        self.assertNotIn("restricted", spins)

    def test_alpha_and_beta_counts(self):
        alpha = [o for o in self.p.data["orbital_energies"] if o["spin"] == "alpha"]
        beta = [o for o in self.p.data["orbital_energies"] if o["spin"] == "beta"]
        self.assertEqual(len(alpha), 20)
        self.assertEqual(len(beta), 18)

    def test_beta_manifold_not_dropped(self):
        beta = [o for o in self.p.data["orbital_energies"] if o["spin"] == "beta"]
        self.assertGreater(len(beta), 0, "beta orbitals must not be dropped")

    def test_triplet_occupation_difference(self):
        """Triplet O2 has 2 more alpha electrons than beta (9 vs 7 occupied)."""
        alpha_occ = sum(
            1
            for o in self.p.data["orbital_energies"]
            if o["spin"] == "alpha" and o["type"] == "occupied"
        )
        beta_occ = sum(
            1
            for o in self.p.data["orbital_energies"]
            if o["spin"] == "beta" and o["type"] == "occupied"
        )
        self.assertEqual(alpha_occ, 9)
        self.assertEqual(beta_occ, 7)
        self.assertEqual(alpha_occ - beta_occ, 2)

    def test_alpha_homo_negative(self):
        alpha_occ = [
            o
            for o in self.p.data["orbital_energies"]
            if o["spin"] == "alpha" and o["type"] == "occupied"
        ]
        homo = max(alpha_occ, key=lambda o: o["energy_eh"])
        self.assertLess(homo["energy_eh"], 0.0)

    # --- MO coefficients (o2-opt_orb.out prints the molecular orbital block) ---

    def test_mo_coeffs_present(self):
        self.assertGreater(len(self.p.data["mo_coeffs"]), 0)

    def test_mo_coeffs_both_spins(self):
        spins = {v["spin"] for v in self.p.data["mo_coeffs"].values()}
        self.assertIn("alpha", spins)
        self.assertIn("beta", spins)

    def test_mo_coeffs_spin_keys(self):
        """Keys are namespaced by spin so alpha/beta MOs never collide."""
        keys = list(self.p.data["mo_coeffs"].keys())
        self.assertTrue(any(k.endswith("_alpha") for k in keys))
        self.assertTrue(any(k.endswith("_beta") for k in keys))

    def test_mo_coeffs_have_coefficients(self):
        for mo in self.p.data["mo_coeffs"].values():
            self.assertIn("coeffs", mo)
            self.assertGreater(len(mo["coeffs"]), 0)


# ---------------------------------------------------------------------------
# CHELPG charges  —  ORCA 6 (chelpg-test.inp.out, water)
# ---------------------------------------------------------------------------


class TestChelpgCharges(unittest.TestCase):
    """chelpg-test.inp.out — regression: real ORCA labels the block
    "CHELPG Charges", not "CHELPG ATOMIC CHARGES". The marker used to miss it
    entirely, silently dropping all CHELPG charges."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("chelpg-test.inp.out")

    def test_chelpg_present(self):
        self.assertIn("CHELPG", self.p.data["charges"])

    def test_chelpg_count(self):
        self.assertEqual(len(self.p.data["charges"]["CHELPG"]), 3)

    def test_chelpg_values(self):
        chg = self.p.data["charges"]["CHELPG"]
        self.assertEqual(chg[0]["atom_sym"], "O")
        self.assertAlmostEqual(chg[0]["charge"], -0.709308, places=4)
        self.assertAlmostEqual(chg[1]["charge"], 0.354383, places=4)

    def test_chelpg_sum_near_zero(self):
        total = sum(c["charge"] for c in self.p.data["charges"]["CHELPG"])
        self.assertAlmostEqual(total, 0.0, places=3)


# ---------------------------------------------------------------------------
# MBIS charges  —  ORCA 6 (mbis-test.out, water)
# ---------------------------------------------------------------------------


class TestMbisCharges(unittest.TestCase):
    """mbis-test.out — MBIS ANALYSIS block (charge / population / spin)."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("mbis-test.out")

    def test_mbis_present(self):
        self.assertIn("MBIS", self.p.data["charges"])

    def test_mbis_count(self):
        self.assertEqual(len(self.p.data["charges"]["MBIS"]), 3)

    def test_mbis_charge_and_population(self):
        o = self.p.data["charges"]["MBIS"][0]
        self.assertEqual(o["atom_sym"], "O")
        self.assertAlmostEqual(o["charge"], -0.788429, places=4)
        self.assertAlmostEqual(o["population"], 8.788429, places=4)

    def test_mbis_sum_near_zero(self):
        total = sum(c["charge"] for c in self.p.data["charges"]["MBIS"])
        self.assertAlmostEqual(total, 0.0, places=3)


# ---------------------------------------------------------------------------
# NMR spin-spin couplings  —  ORCA 6 (nmr-coupling-test.out, methane)
# ---------------------------------------------------------------------------


class TestNmrCouplings(unittest.TestCase):
    """nmr-coupling-test.out — SUMMARY OF ISOTROPIC COUPLING CONSTANTS J (Hz).

    Methane (C + 4H = 5 nuclei) -> C(5,2) = 10 unique pairwise couplings.
    """

    @classmethod
    def setUpClass(cls):
        cls.p = _load("nmr-coupling-test.out")

    def test_shielding_count(self):
        self.assertEqual(len(self.p.data["nmr_shielding"]), 5)

    def test_coupling_count(self):
        self.assertEqual(len(self.p.data["nmr_couplings"]), 10)

    def test_couplings_unique_ordered_pairs(self):
        for c in self.p.data["nmr_couplings"]:
            self.assertLess(c["atom_idx1"], c["atom_idx2"])

    def test_ch_coupling_magnitude(self):
        """1J(C,H) couplings (carbon idx 0 to the four H) should all be ~165 Hz."""
        ch = [
            c["coupling"] for c in self.p.data["nmr_couplings"] if c["atom_idx1"] == 0
        ]
        self.assertEqual(len(ch), 4)
        for j in ch:
            self.assertAlmostEqual(j, 165.117, places=1)


# ---------------------------------------------------------------------------
# NBO / Natural Population Analysis  —  ORCA 6 + NBO (nbo-test.out, water)
# ---------------------------------------------------------------------------


class TestNboCharges(unittest.TestCase):
    """nbo-test.out — "Summary of Natural Population Analysis" table.

    Parser reads the detailed table (Charge / Core / Valence / Rydberg /
    Total) and stores 0-based atom indices.
    """

    @classmethod
    def setUpClass(cls):
        cls.p = _load("nbo-test.out")

    def test_nbo_present(self):
        self.assertIn("NBO", self.p.data["charges"])

    def test_nbo_count(self):
        self.assertEqual(len(self.p.data["charges"]["NBO"]), 3)

    def test_nbo_indices_zero_based(self):
        idxs = [c["atom_idx"] for c in self.p.data["charges"]["NBO"]]
        self.assertEqual(idxs, [0, 1, 2])

    def test_nbo_charge_values(self):
        o = self.p.data["charges"]["NBO"][0]
        self.assertEqual(o["atom_sym"], "O")
        self.assertAlmostEqual(o["charge"], -0.88482, places=4)
        self.assertAlmostEqual(
            self.p.data["charges"]["NBO"][1]["charge"], 0.44241, places=4
        )

    def test_nbo_decomposition_fields(self):
        o = self.p.data["charges"]["NBO"][0]
        self.assertAlmostEqual(o["core"], 1.99996, places=4)
        self.assertAlmostEqual(o["valence"], 6.86890, places=4)
        self.assertAlmostEqual(o["rydberg"], 0.01596, places=4)
        self.assertAlmostEqual(o["total"], 8.88482, places=4)

    def test_nbo_sum_near_zero(self):
        total = sum(c["charge"] for c in self.p.data["charges"]["NBO"])
        self.assertAlmostEqual(total, 0.0, places=3)


# ---------------------------------------------------------------------------
# RESP charges  —  ORCA 6 (resp-test.out, water)
# ---------------------------------------------------------------------------


class TestRespCharges(unittest.TestCase):
    """resp-test.out — regression: ORCA 6 supports RESP and labels the block
    "RESP Charges", not "RESP ATOMIC CHARGES". The marker used to miss it
    entirely, silently dropping all RESP charges (same bug class as CHELPG)."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("resp-test.out")

    def test_resp_present(self):
        self.assertIn("RESP", self.p.data["charges"])

    def test_resp_count(self):
        self.assertEqual(len(self.p.data["charges"]["RESP"]), 3)

    def test_resp_values(self):
        chg = self.p.data["charges"]["RESP"]
        self.assertEqual(chg[0]["atom_sym"], "O")
        self.assertAlmostEqual(chg[0]["charge"], -0.709028, places=4)
        self.assertAlmostEqual(chg[1]["charge"], 0.354242, places=4)

    def test_resp_sum_near_zero(self):
        total = sum(c["charge"] for c in self.p.data["charges"]["RESP"])
        self.assertAlmostEqual(total, 0.0, places=3)


# ---------------------------------------------------------------------------
# Spin contamination <S**2>  &  dispersion correction
# ---------------------------------------------------------------------------


class TestSpinContaminationAndDispersion(unittest.TestCase):
    def test_open_shell_s2_present(self):
        """o2-opt_orb.out (triplet) reports <S**2>."""
        p = _load("o2-opt_orb.out")
        s2 = p.data["spin_s2"]
        self.assertIsNotNone(s2)
        self.assertAlmostEqual(s2["actual"], 2.007028, places=4)
        self.assertAlmostEqual(s2["ideal"], 2.0, places=4)
        self.assertAlmostEqual(s2["contamination"], 0.007028, places=4)

    def test_closed_shell_s2_none(self):
        """Restricted benzene has no <S**2> block."""
        p = _load("benzene-opt.out")
        self.assertIsNone(p.data["spin_s2"])

    def test_dispersion_present(self):
        """ethane-scan.out uses a dispersion-corrected functional."""
        p = _load("ethane-scan.out")
        self.assertAlmostEqual(p.data["dispersion"], -0.000882087, places=6)

    def test_dispersion_absent_is_none(self):
        p = _load("benzene-opt.out")
        self.assertIsNone(p.data["dispersion"])


# ---------------------------------------------------------------------------
# Mayer bond orders  —  ORCA 6
# ---------------------------------------------------------------------------


class TestMayerBondOrders(unittest.TestCase):
    def test_o2_single_bond(self):
        """O2 has one O-O bond with double-bond character (~1.8)."""
        p = _load("o2-opt_orb.out")
        bonds = p.data["mayer_bond_orders"]
        self.assertEqual(len(bonds), 1)
        self.assertEqual(bonds[0]["atom_idx1"], 0)
        self.assertEqual(bonds[0]["atom_idx2"], 1)
        self.assertGreater(bonds[0]["order"], 1.5)

    def test_benzene_bond_count(self):
        """Benzene: 6 aromatic C-C + 6 C-H = 12 bonds above threshold."""
        p = _load("benzene-opt.out")
        self.assertEqual(len(p.data["mayer_bond_orders"]), 12)

    def test_benzene_aromatic_cc_order(self):
        """Aromatic C-C bond orders should be ~1.4."""
        p = _load("benzene-opt.out")
        cc = [
            b
            for b in p.data["mayer_bond_orders"]
            if b["atom_sym1"] == "C" and b["atom_sym2"] == "C"
        ]
        self.assertEqual(len(cc), 6)
        for b in cc:
            self.assertAlmostEqual(b["order"], 1.39, delta=0.05)

    def test_indices_ordered(self):
        p = _load("benzene-opt.out")
        for b in p.data["mayer_bond_orders"]:
            self.assertLess(b["atom_idx1"], b["atom_idx2"])


# ---------------------------------------------------------------------------
# NBO orbital list + second-order perturbation  —  ORCA 6 + NBO (nbo-test.out)
# ---------------------------------------------------------------------------


class TestNboAnalysis(unittest.TestCase):
    """nbo-test.out (water) — NATURAL BOND ORBITALS summary + E(2) analysis."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("nbo-test.out")

    def test_nbo_orbital_count(self):
        self.assertEqual(len(self.p.data["nbo_orbitals"]), 24)

    def test_nbo_orbital_types(self):
        types = {o["type"] for o in self.p.data["nbo_orbitals"]}
        self.assertEqual(types, {"CR", "LP", "BD", "BD*", "RY"})

    def test_nbo_bond_orbital_fields(self):
        bd = next(o for o in self.p.data["nbo_orbitals"] if o["type"] == "BD")
        self.assertEqual(bd["index"], 4)
        self.assertEqual(bd["atoms"], "O  1- H  2")
        self.assertAlmostEqual(bd["occupancy"], 1.99864, places=4)
        self.assertAlmostEqual(bd["energy"], -0.60737, places=4)

    def test_core_orbital_high_occupancy(self):
        cr = next(o for o in self.p.data["nbo_orbitals"] if o["type"] == "CR")
        self.assertAlmostEqual(cr["occupancy"], 2.0, delta=0.001)

    def test_nbo_orbital_atom_indices(self):
        """atom_indices (0-based) are extracted for 3D highlighting."""
        bd = next(o for o in self.p.data["nbo_orbitals"] if o["type"] == "BD")
        self.assertEqual(bd["atom_indices"], [0, 1])  # O 1 - H 2 -> [0, 1]
        cr = next(o for o in self.p.data["nbo_orbitals"] if o["type"] == "CR")
        self.assertEqual(cr["atom_indices"], [0])  # O 1 -> [0]

    def test_nbo_bond_hybridization(self):
        """BD O-H: O is ~sp3 hybridized, H is pure s (water)."""
        bd = next(o for o in self.p.data["nbo_orbitals"] if o["type"] == "BD")
        hybrids = bd["hybrids"]
        self.assertEqual(len(hybrids), 2)
        o_comp = next(h for h in hybrids if h["atom_sym"] == "O")
        h_comp = next(h for h in hybrids if h["atom_sym"] == "H")
        self.assertAlmostEqual(o_comp["s_pct"], 22.56, places=1)
        self.assertAlmostEqual(o_comp["p_pct"], 77.29, places=1)
        self.assertTrue(o_comp["label"].startswith("sp"))
        self.assertEqual(h_comp["label"], "s")

    def test_nbo_hybrid_raw_text_captured(self):
        """Each hybrid component keeps the verbatim ORCA text for the detail view."""
        bd = next(o for o in self.p.data["nbo_orbitals"] if o["type"] == "BD")
        o_comp = next(h for h in bd["hybrids"] if h["atom_sym"] == "O")
        self.assertTrue(o_comp["raw"].startswith("s("))
        self.assertIn("22.56%", o_comp["raw"])
        self.assertIn("77.29%", o_comp["raw"])

    def test_nbo_core_and_lone_pair_hybrids(self):
        cr = next(o for o in self.p.data["nbo_orbitals"] if o["type"] == "CR")
        self.assertEqual(cr["hybrids"][0]["label"], "s")  # 1s core
        lps = [o for o in self.p.data["nbo_orbitals"] if o["type"] == "LP"]
        labels = {o["hybrids"][0]["label"] for o in lps if o["hybrids"]}
        self.assertIn("p", labels)  # the pure-p lone pair

    def test_nbo_hybrid_labels_not_explosive(self):
        """Diffuse Rydberg orbitals (tiny %s) must not yield sp^huge labels."""
        import re

        for o in self.p.data["nbo_orbitals"]:
            for h in o.get("hybrids", []):
                for num in re.findall(r"[\d.]+", h["label"]):
                    self.assertLess(
                        float(num), 20.0, f"explosive hybrid label: {h['label']!r}"
                    )

    def test_perturbation_count(self):
        self.assertEqual(len(self.p.data["nbo_perturbation"]), 6)

    def test_perturbation_fields(self):
        first = self.p.data["nbo_perturbation"][0]
        self.assertIn("LP", first["donor"])
        self.assertIn("RY", first["acceptor"])
        self.assertAlmostEqual(first["e2_kcal"], 1.95, places=2)
        self.assertAlmostEqual(first["e_diff"], 1.87, places=2)
        self.assertAlmostEqual(first["fock"], 0.054, places=3)

    def test_non_nbo_file_empty(self):
        p = _load("benzene-opt.out")
        self.assertEqual(p.data["nbo_orbitals"], [])
        self.assertEqual(p.data["nbo_perturbation"], [])


# ---------------------------------------------------------------------------
# Post-HF energy components  —  ORCA 6 (mp2-test.out, ccsd-test.out)
# ---------------------------------------------------------------------------


class TestEnergyComponentsMp2(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = _load("mp2-test.out")

    def _value(self, label):
        for c in self.p.data["energy_components"]:
            if c["label"] == label:
                return c["value"]
        return None

    def test_mp2_correlation(self):
        self.assertAlmostEqual(
            self._value("MP2 correlation energy"), -0.201139520, places=6
        )

    def test_mp2_total(self):
        self.assertAlmostEqual(self._value("MP2 total energy"), -76.162123498, places=6)

    def test_no_ccsd_fields(self):
        labels = {c["label"] for c in self.p.data["energy_components"]}
        self.assertNotIn("E(CCSD(T))", labels)


class TestEnergyComponentsCcsd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = _load("ccsd-test.out")

    def _value(self, label):
        for c in self.p.data["energy_components"]:
            if c["label"] == label:
                return c["value"]
        return None

    def test_reference_energy(self):
        self.assertAlmostEqual(
            self._value("Reference energy E(0)"), -75.960983979, places=6
        )

    def test_ccsd_correlation(self):
        self.assertAlmostEqual(
            self._value("CCSD correlation energy"), -0.210803526, places=6
        )

    def test_triples_correction(self):
        self.assertAlmostEqual(
            self._value("(T) triples correction"), -0.002885334, places=6
        )

    def test_ccsd_t_total(self):
        self.assertAlmostEqual(self._value("E(CCSD(T))"), -76.174672838, places=6)

    def test_t1_diagnostic_is_dimensionless(self):
        comp = next(
            c for c in self.p.data["energy_components"] if c["label"] == "T1 diagnostic"
        )
        self.assertTrue(comp["dimensionless"])
        self.assertAlmostEqual(comp["value"], 0.005944123, places=6)

    def test_ccsd_t_is_final_single_point(self):
        """E(CCSD(T)) should equal the final single point energy."""
        self.assertAlmostEqual(
            self._value("E(CCSD(T))"), self.p.data["scf_energy"], places=6
        )


class TestEnergyComponentsAbsentForDft(unittest.TestCase):
    def test_dft_has_no_components(self):
        p = _load("benzene-opt.out")
        self.assertEqual(p.data["energy_components"], [])


class TestValidOutputsParseQuietly(unittest.TestCase):
    """A clean ORCA file should not make the parser complain.

    Several blocks are scanned line by line with a numeric conversion used to
    tell data rows from the headers and footers around them, so the resulting
    ValueError is ordinary control flow. Reporting it at WARNING put ten lines
    into the log for twenty perfectly good files, which is how a log stops
    being worth reading -- and how a genuine warning goes unnoticed.
    """

    def test_no_warnings_are_logged_for_any_sample(self):
        import glob
        import logging

        records = []

        class _Capture(logging.Handler):
            def emit(self, record):
                if record.levelno >= logging.WARNING:
                    records.append(
                        f"{os.path.basename(record.pathname)}:{record.lineno} "
                        f"{record.getMessage()[:70]}"
                    )

        handler = _Capture()
        root = logging.getLogger()
        root.addHandler(handler)
        previous = root.level
        root.setLevel(logging.WARNING)
        self.addCleanup(root.setLevel, previous)
        self.addCleanup(root.removeHandler, handler)

        samples = sorted(glob.glob(os.path.join(_SAMPLES, "*.out")))
        self.assertGreater(len(samples), 10, "sample corpus went missing")
        for path in samples:
            _load(os.path.basename(path))

        self.assertEqual(
            records, [], "parser warned on valid input:\n" + "\n".join(records)
        )


if __name__ == "__main__":
    unittest.main()

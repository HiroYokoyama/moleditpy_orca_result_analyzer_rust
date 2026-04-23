"""
tests/test_parser_samples.py
Integration tests for OrcaParser using real ORCA output files.

All tests load actual .out files from tests/sample_outputs/ and assert
concrete values extracted from those files.
"""

import os
import sys
import importlib.util
import unittest

# ---------------------------------------------------------------------------
# Bootstrap: load parser without Qt
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(__file__)
_SAMPLES = os.path.join(_HERE, "sample_outputs")
_PARSER_SRC = os.path.normpath(
    os.path.join(_HERE, "..", "orca_result_analyzer_rust", "parser.py")
)


def _load_parser():
    spec = importlib.util.spec_from_file_location("orca_parser_standalone", _PARSER_SRC)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["orca_parser_standalone"] = mod
    spec.loader.exec_module(mod)
    return mod


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
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930734972084, places=6
        )

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
        self.assertAlmostEqual(x,  0.805933, places=4)
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
                self.assertAlmostEqual(fc[j], dc[j], places=5,
                    msg=f"atom {i} coord {j} mismatch")

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
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930718933873, places=6
        )

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_final_coords_first_atom(self):
        x, y, z = self.p.data["coords"][0]
        self.assertAlmostEqual(x,  0.805991, places=4)
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
        self.assertAlmostEqual(y,  0.247667, places=4)
        self.assertAlmostEqual(z,  0.079907, places=4)

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
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930784926115, places=6
        )

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
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930718955827, places=6
        )

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
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930784926115, places=6
        )

    def test_nmr_shielding_present(self):
        self.assertGreater(len(self.p.data["nmr_shielding"]), 0)

    def test_nmr_shielding_count(self):
        # 6 C + 6 H = 12 nuclei
        self.assertEqual(len(self.p.data["nmr_shielding"]), 12)

    def test_nmr_carbon_shielding_range(self):
        """Benzene carbon shielding should be roughly 40-80 ppm."""
        carbon_shields = [
            e["shielding"] for e in self.p.data["nmr_shielding"]
            if e["atom_sym"] == "C"
        ]
        self.assertEqual(len(carbon_shields), 6)
        for s in carbon_shields:
            self.assertGreater(s, 0.0)

    def test_nmr_hydrogen_shielding_range(self):
        """Benzene H shielding should be positive."""
        h_shields = [
            e["shielding"] for e in self.p.data["nmr_shielding"]
            if e["atom_sym"] == "H"
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
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930718955827, places=6
        )

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
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.728888949460, places=6
        )

    def test_tddft_state_count(self):
        self.assertEqual(len(self.p.data["tddft"]), 5)

    def test_tddft_first_state_ev(self):
        self.assertAlmostEqual(
            self.p.data["tddft"][0]["energy_ev"], 5.494, places=2
        )

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
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.728828351432, places=6
        )

    def test_tddft_state_count(self):
        self.assertEqual(len(self.p.data["tddft"]), 5)

    def test_tddft_first_state_ev(self):
        self.assertAlmostEqual(
            self.p.data["tddft"][0]["energy_ev"], 5.494, places=2
        )

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
            for key in ("index", "occupation", "energy_eh", "energy_ev", "spin", "type"):
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
            for key in ("index", "occupation", "energy_eh", "energy_ev", "spin", "type"):
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
        cycle_labels = [l for l in labels if "Cycle" in l or "Initial" in l or "Final" in l or "Post" in l or "Property" in l]
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
        cycle_labels = [l for l in labels if "Cycle" in l or "Initial" in l or "Final" in l or "Post" in l or "Property" in l]
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
        c_vals = [e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "C"]
        for v in c_vals:
            self.assertGreater(v, 10.0)
            self.assertLess(v, 200.0)

    def test_nmr_shielding_hydrogen_value_range(self):
        """Benzene H shielding ~20-35 ppm (absolute)."""
        h_vals = [e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "H"]
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
        c_vals = [e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "C"]
        self.assertEqual(len(c_vals), 6)
        for v in c_vals:
            self.assertGreater(v, 10.0)
            self.assertLess(v, 200.0)

    def test_nmr_hydrogen_shielding_range(self):
        h_vals = [e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "H"]
        self.assertEqual(len(h_vals), 6)
        for v in h_vals:
            self.assertGreater(v, 10.0)
            self.assertLess(v, 50.0)

    def test_nmr_shielding_consistent_with_orca6(self):
        """Mean C shielding diff between ORCA 5 and 6 should be < 5 ppm."""
        p6 = _load("benzene-opt-nmr.out")
        c5 = [e["shielding"] for e in self.p.data["nmr_shielding"] if e["atom_sym"] == "C"]
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
        scan_steps = [s for s in self.p.data.get("scan_steps", []) if s["type"] == "scan_step"]
        self.assertEqual(len(scan_steps), 5)

    def test_total_trajectory_steps(self):
        """Should have 37 total steps in the trajectory (scan steps + intermediate opt cycles)."""
        self.assertEqual(len(self.p.data.get("scan_steps", [])), 37)

    def test_scan_step_coords_and_energies(self):
        scan_steps = [s for s in self.p.data.get("scan_steps", []) if s["type"] == "scan_step"]
        self.assertEqual(scan_steps[0]["scan_coord"], -60.0)
        self.assertAlmostEqual(scan_steps[0]["energy"], -79.791846796153)
        self.assertEqual(scan_steps[-1]["scan_coord"], 60.0)
        self.assertAlmostEqual(scan_steps[-1]["energy"], -79.791847846936)
        
        for step in scan_steps:
            self.assertEqual(len(step["atoms"]), 8)
            self.assertEqual(len(step["coords"]), 8)
            self.assertIn("scan_coord", step)
            self.assertIn("energy", step)


if __name__ == "__main__":
    unittest.main()

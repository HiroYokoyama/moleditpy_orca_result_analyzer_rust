"""
tests/test_cube_and_mo.py
Tests for:
  - CubeVisualizer._parse_cube  (vis.py)  using benzene-opt-ene_MO_21.cube
  - parse_orbital_energies      (parser.py) using benzene-opt-ene.out
  - parse_mo_coeffs             (parser.py) using benzene-opt-ene.out

No Qt / PyVista is needed: _parse_cube is tested via a lightweight
standalone reader, and the parser is loaded directly.
"""

import os
import sys
import importlib.util
import unittest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(__file__)
_SAMPLES = os.path.join(_HERE, "sample_outputs")
_CUBE_FILE = os.path.join(_SAMPLES, "benzene-opt-ene_cubes", "benzene-opt-ene_MO_21.cube")
_ENE_FILE = os.path.join(_SAMPLES, "benzene-opt-ene.out")

_PARSER_SRC = os.path.normpath(
    os.path.join(_HERE, "..", "orca_result_analyzer_rust", "parser.py")
)
_VIS_SRC = os.path.normpath(
    os.path.join(_HERE, "..", "orca_result_analyzer_rust", "vis.py")
)

# ---------------------------------------------------------------------------
# Bootstrap: parser (no Qt needed)
# ---------------------------------------------------------------------------

def _load_parser():
    spec = importlib.util.spec_from_file_location("orca_parser_cube_mo", _PARSER_SRC)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["orca_parser_cube_mo"] = mod
    spec.loader.exec_module(mod)
    return mod

_parser_mod = _load_parser()
OrcaParser = _parser_mod.OrcaParser


def _load_ene():
    p = OrcaParser()
    with open(_ENE_FILE, encoding="utf-8", errors="replace") as f:
        content = f.read()
    p.load_from_memory(content, _ENE_FILE)
    return p


# ---------------------------------------------------------------------------
# Standalone cube reader (mirrors vis._parse_cube, no PyVista required)
# ---------------------------------------------------------------------------

def _parse_cube_standalone(filepath):
    """Minimal cube reader for testing — equivalent to CubeVisualizer._parse_cube."""
    import numpy as np
    with open(filepath, "r") as f:
        lines = f.readlines()

    tokens = lines[2].split()
    n_atoms = abs(int(tokens[0]))
    origin = [float(x) for x in tokens[1:4]]

    nx_raw, *x_vec = lines[3].split()
    ny_raw, *y_vec = lines[4].split()
    nz_raw, *z_vec = lines[5].split()
    nx, ny, nz = int(nx_raw), int(ny_raw), int(nz_raw)

    start_line = 6 + n_atoms
    full_str = " ".join(lines[start_line:])
    data = [float(v) for v in full_str.split()]

    return {
        "n_atoms": n_atoms,
        "origin": origin,
        "dims": (abs(nx), abs(ny), abs(nz)),
        "x_vec": [float(v) for v in x_vec],
        "y_vec": [float(v) for v in y_vec],
        "z_vec": [float(v) for v in z_vec],
        "data": data,
        "is_angstrom": nx < 0,
    }


# ---------------------------------------------------------------------------
# TestParseCube — benzene-opt-ene_MO_21.cube
# ---------------------------------------------------------------------------

class TestParseCube(unittest.TestCase):
    """Real cube file: HOMO-1 (MO 22) of benzene, 40x40x40 grid in Bohr."""

    @classmethod
    def setUpClass(cls):
        cls.meta = _parse_cube_standalone(_CUBE_FILE)

    # --- header ---

    def test_n_atoms(self):
        self.assertEqual(self.meta["n_atoms"], 12)

    def test_dims(self):
        self.assertEqual(self.meta["dims"], (40, 40, 40))

    def test_origin_x(self):
        self.assertAlmostEqual(self.meta["origin"][0], -12.249253, places=4)

    def test_origin_y(self):
        self.assertAlmostEqual(self.meta["origin"][1], -11.834071, places=4)

    def test_origin_z(self):
        self.assertAlmostEqual(self.meta["origin"][2], -7.653851, places=4)

    # --- voxel step vectors ---

    def test_x_step_x(self):
        self.assertAlmostEqual(self.meta["x_vec"][0], 0.628167, places=5)

    def test_x_step_yz_zero(self):
        self.assertAlmostEqual(self.meta["x_vec"][1], 0.0, places=6)
        self.assertAlmostEqual(self.meta["x_vec"][2], 0.0, places=6)

    def test_y_step_y(self):
        self.assertAlmostEqual(self.meta["y_vec"][1], 0.606875, places=5)

    def test_z_step_z(self):
        self.assertAlmostEqual(self.meta["z_vec"][2], 0.392505, places=5)

    # --- units ---

    def test_bohr_units(self):
        """nx > 0 means Bohr (not Angstrom)."""
        self.assertFalse(self.meta["is_angstrom"])

    # --- volumetric data ---

    def test_data_count(self):
        """40 * 40 * 40 = 64000 data points."""
        self.assertEqual(len(self.meta["data"]), 64000)

    def test_data_not_all_zero(self):
        """MO wavefunction must have non-zero values."""
        self.assertTrue(any(abs(v) > 1e-10 for v in self.meta["data"]))

    def test_data_has_positive_and_negative(self):
        """MO wavefunction has both phases."""
        self.assertTrue(any(v > 0 for v in self.meta["data"]))
        self.assertTrue(any(v < 0 for v in self.meta["data"]))


# ---------------------------------------------------------------------------
# TestParseOrbitalEnergies — benzene-opt-ene.out
# ---------------------------------------------------------------------------

class TestParseOrbitalEnergies(unittest.TestCase):
    """ORBITAL ENERGIES section from benzene single-point (ORCA 6.1.1)."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load_ene()
        cls.oe = cls.p.data["orbital_energies"]

    def test_count(self):
        """21 occupied + 11 virtual (ORCA prints first 10 virtuals) = 32."""
        self.assertEqual(len(self.oe), 32)

    def test_mos_alias_same_length(self):
        self.assertEqual(len(self.p.data["mos"]), len(self.oe))

    def test_first_orbital_index(self):
        self.assertEqual(self.oe[0]["index"], 0)

    def test_first_orbital_energy_eh(self):
        self.assertAlmostEqual(self.oe[0]["energy_eh"], -10.188343, places=4)

    def test_first_orbital_energy_ev(self):
        self.assertAlmostEqual(self.oe[0]["energy_ev"], -277.2389, places=2)

    def test_first_orbital_occupation(self):
        self.assertAlmostEqual(self.oe[0]["occupation"], 2.0, places=4)

    def test_first_orbital_type_occupied(self):
        self.assertEqual(self.oe[0]["type"], "occupied")

    def test_first_orbital_spin(self):
        self.assertEqual(self.oe[0]["spin"], "restricted")

    def test_homo_index(self):
        """HOMO = orbital 20 (last doubly occupied for benzene 21-electron pi system)."""
        self.assertEqual(self.oe[20]["index"], 20)

    def test_homo_occupied(self):
        self.assertAlmostEqual(self.oe[20]["occupation"], 2.0, places=4)

    def test_homo_energy_eh(self):
        self.assertAlmostEqual(self.oe[20]["energy_eh"], -0.252649, places=4)

    def test_lumo_index(self):
        self.assertEqual(self.oe[21]["index"], 21)

    def test_lumo_virtual(self):
        self.assertAlmostEqual(self.oe[21]["occupation"], 0.0, places=4)

    def test_lumo_type_virtual(self):
        self.assertEqual(self.oe[21]["type"], "virtual")

    def test_lumo_energy_eh(self):
        self.assertAlmostEqual(self.oe[21]["energy_eh"], -0.005861, places=5)

    def test_lumo_energy_ev(self):
        self.assertAlmostEqual(self.oe[21]["energy_ev"], -0.1595, places=3)

    def test_homo_lumo_gap_positive(self):
        homo_ev = self.oe[20]["energy_ev"]
        lumo_ev = self.oe[21]["energy_ev"]
        self.assertGreater(lumo_ev, homo_ev)

    def test_occupied_count(self):
        occ = [o for o in self.oe if o["type"] == "occupied"]
        self.assertEqual(len(occ), 21)

    def test_virtual_count(self):
        virt = [o for o in self.oe if o["type"] == "virtual"]
        self.assertEqual(len(virt), 11)

    def test_energy_eh_and_energy_same(self):
        """'energy' key must equal energy_eh for backward compatibility."""
        for orb in self.oe:
            self.assertAlmostEqual(orb["energy"], orb["energy_eh"], places=6)


# ---------------------------------------------------------------------------
# TestParseMoCoeffs — benzene-opt-ene.out
# ---------------------------------------------------------------------------

class TestParseMoCoeffs(unittest.TestCase):
    """MOLECULAR ORBITALS (RHF, ROHF) section from benzene single-point."""

    @classmethod
    def setUpClass(cls):
        cls.p = _load_ene()
        cls.mc = cls.p.data["mo_coeffs"]

    def test_total_mo_count(self):
        """Benzene with def2-SVP has 114 basis functions → 114 MOs."""
        self.assertEqual(len(self.mc), 114)

    def test_all_keys_restricted(self):
        """All MOs are restricted (RHF)."""
        for key in self.mc:
            self.assertTrue(key.endswith("_restricted"), f"Unexpected key: {key}")

    def test_mo0_present(self):
        self.assertIn("0_restricted", self.mc)

    def test_mo0_energy(self):
        self.assertAlmostEqual(self.mc["0_restricted"]["energy"], -10.18834, places=4)

    def test_mo0_occupation(self):
        self.assertAlmostEqual(self.mc["0_restricted"]["occ"], 2.0, places=4)

    def test_mo0_spin(self):
        self.assertEqual(self.mc["0_restricted"]["spin"], "restricted")

    def test_mo0_coeff_count(self):
        """114 basis functions → 114 coefficients per MO."""
        self.assertEqual(len(self.mc["0_restricted"]["coeffs"]), 114)

    def test_mo0_first_coeff_atom(self):
        c = self.mc["0_restricted"]["coeffs"][0]
        self.assertEqual(c["atom_idx"], 0)
        self.assertEqual(c["sym"], "C")
        self.assertEqual(c["orb"], "1s")

    def test_mo0_first_coeff_value(self):
        c = self.mc["0_restricted"]["coeffs"][0]
        self.assertAlmostEqual(c["coeff"], 0.406941, places=5)

    def test_homo_mo20_present(self):
        self.assertIn("20_restricted", self.mc)

    def test_homo_mo20_occupied(self):
        self.assertAlmostEqual(self.mc["20_restricted"]["occ"], 2.0, places=4)

    def test_lumo_mo21_present(self):
        self.assertIn("21_restricted", self.mc)

    def test_lumo_mo21_virtual(self):
        self.assertAlmostEqual(self.mc["21_restricted"]["occ"], 0.0, places=4)

    def test_all_mos_have_coeffs(self):
        """Every parsed MO must have at least one coefficient."""
        for key, mo in self.mc.items():
            self.assertGreater(len(mo["coeffs"]), 0, f"MO {key} has no coeffs")

    def test_coeff_fields_present(self):
        """Every coefficient entry must have atom_idx, sym, orb, coeff keys."""
        sample = self.mc["0_restricted"]["coeffs"][0]
        for field in ("atom_idx", "sym", "orb", "coeff"):
            self.assertIn(field, sample)


if __name__ == "__main__":
    unittest.main()

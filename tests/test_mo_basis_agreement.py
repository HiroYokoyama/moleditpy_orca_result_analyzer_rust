"""
tests/test_mo_basis_agreement.py

Regenerates MO 21 of benzene through BasisSetEngine from two ORCA outputs of
the same optimized geometry -- def2-SVP and def2-QZVP -- and checks that the
two orbitals describe the same physical state.

What this covers: the def2-QZVP job carries f and g shells that def2-SVP does
not, so regenerating from it drives the high-angular-momentum branches of
_prepare_definitions/_precompute_shells end to end, from ORCA text through the
parser to grid values, and pins the result against a stored reference. A basis
mismatch, a dropped shell or a coefficient misalignment shows up here.

What this does NOT cover: the *accuracy* of the f/g normalization. No orbital
of this molecule can. Benzene is all light atoms, so def2-QZVP's f and g
functions are polarization functions carrying almost no weight -- 0.27% of MO
21's coefficient weight, and at most 0.49% across all 21 occupied orbitals
(MO 16, used here). Perturbing a g polynomial by 3% moves either orbital's
grid by ~3e-6 relative, below the ~1e-5 text precision of the stored cubes:
there is no headroom to assert on, and every assertion in this module stays
green. A system-level f/g check needs a molecule whose occupied orbitals carry
real high-l weight, not another orbital from this one.

That accuracy is covered directly instead, in test_mo_engine.py, by three
component-level tests that were verified against injected faults: a g0 shape
error trips test_g_z4_nodal_cones, a g+-2 shape error trips
test_components_of_a_shell_are_mutually_orthogonal, and a per-component scale
error trips test_g_components_are_unit_normalized.

Nothing is written to disk -- the grids are evaluated in memory. The stored
cubes are read only for their grid definition and as regression references.
"""

import importlib.util
import os
import sys
import types
import unittest
from unittest.mock import MagicMock

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from _parser_loader import load_standalone_parser  # noqa: E402

_HERE = os.path.dirname(__file__)
_SAMPLES = os.path.join(_HERE, "sample_outputs")
_PKG = os.path.normpath(os.path.join(_HERE, "..", "orca_result_analyzer_rust"))

# The engine works in Bohr; the parser reports shell centers in Angstrom.
BOHR_TO_ANG = 0.529177249

# Both jobs are single points on the same optimized geometry, so the plugin
# writes both cubes on an identical grid -- no interpolation is needed to
# compare them pointwise.
_CASES = {
    "svp": (
        "benzene-opt-ene.out",
        os.path.join("benzene-opt-ene_cubes", "benzene-opt-ene_MO_21.cube"),
    ),
    "qzvp": (
        "benzene-opt-eneQZ.out",
        os.path.join("benzene-opt-eneQZ_cubes", "benzene-opt-eneQZ_MO_21.cube"),
    ),
}

# MO 21 is the LUMO -- the orbital the stored reference cubes hold, so it is
# the one that can be checked against them. MO 16 is added as the occupied
# counterpart: an occupied valence orbital is described almost identically by
# both bases (|cos| 0.9994 against MO 21's 0.979), which affords a far tighter
# threshold, and it carries the largest g weight of any occupied orbital here.
_REFERENCED_MO = 21
_OCCUPIED_MO = 16


def _load_module(name, filename):
    """Load one plugin module against a minimal QtCore stub.

    mo_engine needs QThread to be a real, subclassable class for CalcWorker's
    class body to execute. The stub is torn down afterwards so it cannot leak
    into test modules that install their own, richer PySide6 stand-ins.
    """

    class _QThread:
        def __init__(self, *a, **k):
            pass

    core = types.ModuleType("PySide6.QtCore")
    core.QThread = _QThread
    core.Signal = lambda *a, **k: MagicMock()
    pyside6 = types.ModuleType("PySide6")
    pyside6.QtCore = core

    saved = {k: sys.modules.get(k) for k in ("PySide6", "PySide6.QtCore")}
    sys.modules["PySide6"] = pyside6
    sys.modules["PySide6.QtCore"] = core
    try:
        spec = importlib.util.spec_from_file_location(
            name, os.path.join(_PKG, filename)
        )
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return mod
    finally:
        for key, value in saved.items():
            if value is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = value


_parser_mod = load_standalone_parser("orca_parser_basis_agreement")
_engine_mod = _load_module("orca_mo_engine_basis_agreement", "mo_engine.py")
OrcaParser = _parser_mod.OrcaParser
BasisSetEngine = _engine_mod.BasisSetEngine


def _read_cube(path):
    """Grid definition and volumetric data of a cube file.

    Deliberately independent of vis._parse_cube: this module is checking the
    numbers, so it should not inherit a bug from the reader under test.
    """
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()

    tokens = lines[2].split()
    n_atoms = abs(int(tokens[0]))
    origin = np.array([float(v) for v in tokens[1:4]])

    header = [lines[i].split() for i in (3, 4, 5)]
    dims = tuple(int(h[0]) for h in header)
    vectors = np.array([[float(v) for v in h[1:4]] for h in header])

    values = np.array(
        [float(v) for v in " ".join(lines[6 + n_atoms :]).split()], dtype=float
    )
    return {"origin": origin, "dims": dims, "vectors": vectors, "values": values}


def _grid_points(cube):
    """The cube's sample points as an (N, 3) array in Bohr, in file order.

    Cube files run z fastest, then y, then x -- which is exactly the order
    np.indices produces for the (nx, ny, nz) shape, so the flattened result
    lines up with the file's value sequence.
    """
    idx = np.indices(cube["dims"]).reshape(3, -1).T
    return cube["origin"] + idx @ cube["vectors"]


def _regenerate(out_name, grid, mo_indices):
    """Evaluate the given MOs on `grid`, driving the path the dialog drives."""
    path = os.path.join(_SAMPLES, out_name)
    parser = OrcaParser()
    with open(path, encoding="utf-8", errors="replace") as fh:
        parser.load_from_memory(fh.read(), path)

    shells = []
    for shell in parser.data.get("basis_set_shells", []):
        center_ang = np.array(shell.get("origin", shell.get("center", [0, 0, 0])))
        shells.append(
            {
                "type": shell.get("l", shell.get("type", 0)),
                "center": center_ang / BOHR_TO_ANG,
                "exps": np.array(shell["exps"]),
                "coeffs": np.array(shell["coeffs"]),
            }
        )

    engine = BasisSetEngine(shells)
    coeffs, values = {}, {}
    for mo in mo_indices:
        coeffs[mo] = np.array(
            [
                c["coeff"]
                for c in parser.data["mo_coeffs"][f"{mo}_restricted"]["coeffs"]
            ],
            dtype=float,
        )
        values[mo] = engine.evaluate_mo_on_grid(mo, grid, coeffs[mo])
    return engine, coeffs, values


def _cosine(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


class TestBasisSetAgreement(unittest.TestCase):
    """MOs 16 and 21 regenerated from def2-SVP and def2-QZVP outputs."""

    @classmethod
    def setUpClass(cls):
        cls.cubes = {
            tag: _read_cube(os.path.join(_SAMPLES, cube))
            for tag, (_, cube) in _CASES.items()
        }
        cls.grid = _grid_points(cls.cubes["svp"])
        cls.engines, cls.coeffs, cls.values = {}, {}, {}
        for tag, (out_name, _) in _CASES.items():
            engine, coeffs, values = _regenerate(
                out_name, cls.grid, (_OCCUPIED_MO, _REFERENCED_MO)
            )
            cls.engines[tag], cls.coeffs[tag], cls.values[tag] = engine, coeffs, values

    def test_both_jobs_share_one_grid(self):
        """The pointwise comparisons below are only meaningful if they do."""
        svp, qzvp = self.cubes["svp"], self.cubes["qzvp"]
        self.assertEqual(svp["dims"], qzvp["dims"])
        np.testing.assert_allclose(svp["origin"], qzvp["origin"])
        np.testing.assert_allclose(svp["vectors"], qzvp["vectors"])

    def test_the_two_bases_span_different_angular_momenta(self):
        """Guards the premise: the QZVP path must actually evaluate f and g.

        If a future sample swap quietly dropped the high-l shells, the
        comparison below would still pass while exercising strictly less.
        """
        types_of = {
            tag: {int(s["type"]) for s in engine.shells}
            for tag, engine in self.engines.items()
        }
        self.assertEqual(types_of["svp"], {0, 1, 2})
        self.assertEqual(types_of["qzvp"], {0, 1, 2, 3, 4})

    def test_every_basis_function_gets_a_coefficient(self):
        for tag in _CASES:
            for mo in (_OCCUPIED_MO, _REFERENCED_MO):
                with self.subTest(basis=tag, mo=mo):
                    self.assertEqual(
                        len(self.coeffs[tag][mo]), self.engines[tag].n_basis
                    )

    def test_regeneration_reproduces_the_stored_cube(self):
        """Regression anchor, per basis, before the two are compared."""
        for tag in _CASES:
            with self.subTest(basis=tag):
                self.assertAlmostEqual(
                    _cosine(
                        self.values[tag][_REFERENCED_MO], self.cubes[tag]["values"]
                    ),
                    1.0,
                    places=6,
                )

    def test_svp_and_qzvp_describe_the_same_orbital(self):
        """Same physical orbital, two basis sets: the shapes must match.

        Compared on magnitude because an MO's overall phase is arbitrary --
        ORCA hands these two out with opposite signs. What must hold is that
        the orbitals are parallel, not that they point the same way.

        The occupied orbital is held to a much tighter threshold than the
        LUMO: both bases describe an occupied valence orbital almost
        identically (measured 0.9994), whereas a virtual orbital is shaped
        partly by whatever room the basis leaves it (0.979). Neither residual
        is engine error, so both thresholds sit just below the measured value
        -- tight enough to catch a mis-mapped shell or a shifted coefficient
        vector, loose enough not to chase real basis-set difference.
        """
        for mo, floor in ((_OCCUPIED_MO, 0.999), (_REFERENCED_MO, 0.97)):
            with self.subTest(mo=mo):
                cosine = _cosine(self.values["svp"][mo], self.values["qzvp"][mo])
                self.assertGreater(
                    abs(cosine),
                    floor,
                    f"SVP and QZVP MO {mo} disagree (|cos| = {abs(cosine):.4f}); "
                    "the two jobs should describe the same orbital",
                )

    def test_the_occupied_orbital_is_not_confused_with_its_neighbours(self):
        """MO 16 must match MO 16, and nothing else.

        Occupied orbitals can reorder between basis sets, which would make the
        index-to-index comparison above meaningless while still passing if a
        neighbour happened to look similar. Here the match is unambiguous.
        """
        reference = self.values["svp"][_OCCUPIED_MO]
        engine, _, values = _regenerate(
            _CASES["qzvp"][0], self.grid, range(_OCCUPIED_MO - 2, _OCCUPIED_MO + 3)
        )
        overlaps = {mo: abs(_cosine(reference, v)) for mo, v in values.items()}
        best = max(overlaps, key=overlaps.get)
        self.assertEqual(best, _OCCUPIED_MO)
        others = [v for mo, v in overlaps.items() if mo != _OCCUPIED_MO]
        self.assertLess(max(others), 0.1, f"overlaps: {overlaps}")

    def test_the_orbitals_carry_comparable_weight_on_the_grid(self):
        """A prefactor slip can leave shape intact but change the scale."""
        for mo in (_OCCUPIED_MO, _REFERENCED_MO):
            with self.subTest(mo=mo):
                norms = {
                    tag: float(np.linalg.norm(v[mo])) for tag, v in self.values.items()
                }
                self.assertAlmostEqual(
                    norms["svp"], norms["qzvp"], delta=0.05 * norms["svp"]
                )

    def test_no_cube_is_written_while_regenerating(self):
        """The engine evaluates in memory; only CubeWriter should touch disk."""
        before = set(os.listdir(os.path.join(_SAMPLES, "benzene-opt-eneQZ_cubes")))
        _regenerate(_CASES["qzvp"][0], self.grid[:64], (_REFERENCED_MO,))
        after = set(os.listdir(os.path.join(_SAMPLES, "benzene-opt-eneQZ_cubes")))
        self.assertEqual(before, after)


if __name__ == "__main__":
    unittest.main()

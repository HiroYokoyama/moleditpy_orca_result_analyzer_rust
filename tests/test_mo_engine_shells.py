"""
tests/test_mo_engine_shells.py

Correctness of the Rust GTO engine's angular basis functions.

Engine-local (the Rust crate owns these definitions), but the assertions are
deliberately identical to the pure-Python plugin's
tests/test_mo_engine.py::TestSphericalShellNormalization, so the two engines
are held to the same standard.

Why this exists: a mis-normalized basis function still renders as a plausible
orbital — right symmetry, right sign pattern, lobes roughly in place — so it
cannot be caught by looking at the picture. The f and g tables carried
hand-typed prefactors that were wrong in two different ways for a long time
(f(z^3)'s nodal cone sat at 28.6 degrees instead of 39.2; g2 came out with
about twice the norm of g1).

Requires the compiled extension: run `python build.py` first.
"""

import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from orca_result_analyzer_rust.mo_engine import (
        BasisSetEngine,
        UnsupportedBasisError,
    )
except Exception as exc:  # extension not built
    BasisSetEngine = None
    UnsupportedBasisError = None
    _IMPORT_ERROR = exc
else:
    _IMPORT_ERROR = None


def _shell(l_type, alpha=1.0, coeff=1.0, center=(0.0, 0.0, 0.0)):
    return {
        "type": l_type,
        "center": np.array(center, dtype=float),
        "exps": np.array([alpha], dtype=float),
        "coeffs": np.array([coeff], dtype=float),
    }


@unittest.skipIf(BasisSetEngine is None, f"extension not built: {_IMPORT_ERROR}")
class TestSphericalShellShapes(unittest.TestCase):
    """A pure basis function must be radial * Y_lm.

    Convention-independent: a program's normalization convention fixes an
    overall constant per component, it cannot make the angular part something
    that is not a spherical harmonic.
    """

    # ORCA's pure-shell m ordering: 0, +1, -1, +2, -2, ...
    ORDERS = {
        2: [0, 1, -1, 2, -2],
        3: [0, 1, -1, 2, -2, 3, -3],
        4: [0, 1, -1, 2, -2, 3, -3, 4, -4],
    }

    @staticmethod
    def _real_sph_harm(l, m, theta, phi):
        special = __import__("scipy.special", fromlist=["sph_harm_y"])
        if hasattr(special, "sph_harm_y"):
            y = special.sph_harm_y(l, abs(m), theta, phi)
        else:
            y = special.sph_harm(abs(m), l, phi, theta)
        if m == 0:
            return y.real
        sign = (-1) ** abs(m)
        return math.sqrt(2) * sign * (y.real if m > 0 else y.imag)

    def setUp(self):
        try:
            import scipy.special  # noqa: F401
        except ImportError:
            self.skipTest("scipy not installed")
        rng = np.random.default_rng(11)
        v = rng.normal(size=(2000, 3))
        self.dirs = v / np.linalg.norm(v, axis=1)[:, None]
        self.theta = np.arccos(np.clip(self.dirs[:, 2], -1.0, 1.0))
        self.phi = np.arctan2(self.dirs[:, 1], self.dirs[:, 0])

    def _cos_similarity(self, l_type, component, m):
        eng = BasisSetEngine([_shell(l_type)])
        coeffs = np.zeros(eng.n_basis)
        coeffs[component] = 1.0
        got = np.asarray(eng.evaluate_mo_on_grid(0, self.dirs, coeffs))
        ref = np.array(
            [self._real_sph_harm(l_type, m, t, p) for t, p in zip(self.theta, self.phi)]
        )
        return abs(
            float(np.dot(got, ref)) / (np.linalg.norm(got) * np.linalg.norm(ref))
        )

    def test_every_d_component_matches_its_harmonic(self):
        for component, m in enumerate(self.ORDERS[2]):
            with self.subTest(m=m):
                self.assertAlmostEqual(
                    self._cos_similarity(2, component, m), 1.0, places=6
                )

    def test_every_f_component_matches_its_harmonic(self):
        for component, m in enumerate(self.ORDERS[3]):
            with self.subTest(m=m):
                self.assertAlmostEqual(
                    self._cos_similarity(3, component, m), 1.0, places=6
                )

    def test_every_g_component_matches_its_harmonic(self):
        for component, m in enumerate(self.ORDERS[4]):
            with self.subTest(m=m):
                self.assertAlmostEqual(
                    self._cos_similarity(4, component, m), 1.0, places=6
                )


@unittest.skipIf(BasisSetEngine is None, f"extension not built: {_IMPORT_ERROR}")
class TestNodalAngles(unittest.TestCase):
    """Nodal cones have closed forms and move first when a component is
    distorted, so they are the most sensitive regression check."""

    @staticmethod
    def _sign_changes(l_type, component):
        eng = BasisSetEngine([_shell(l_type)])
        coeffs = np.zeros(eng.n_basis)
        coeffs[component] = 1.0
        thetas = np.linspace(0.5, 89.5, 20000)
        pts = np.stack(
            [
                np.sin(np.radians(thetas)),
                np.zeros_like(thetas),
                np.cos(np.radians(thetas)),
            ],
            axis=1,
        )
        v = np.asarray(eng.evaluate_mo_on_grid(0, pts, coeffs))
        return sorted(thetas[:-1][np.sign(v[:-1]) != np.sign(v[1:])])

    def test_f_z3_nodal_cone(self):
        """f(z^3) changes sign where 5cos^2(theta) = 3."""
        nodes = self._sign_changes(3, 0)
        self.assertEqual(len(nodes), 1)
        self.assertAlmostEqual(
            nodes[0], math.degrees(math.acos(math.sqrt(0.6))), places=1
        )

    def test_g_z4_nodal_cones(self):
        """g(z^4) has two cones, at cos^2(theta) = (3 +/- 2*sqrt(6/5))/7."""
        expected = sorted(
            math.degrees(math.acos(math.sqrt(x)))
            for x in ((3 + 2 * math.sqrt(1.2)) / 7, (3 - 2 * math.sqrt(1.2)) / 7)
        )
        nodes = self._sign_changes(4, 0)
        self.assertEqual(len(nodes), 2)
        for got, want in zip(nodes, expected):
            self.assertAlmostEqual(got, want, places=1)


@unittest.skipIf(BasisSetEngine is None, f"extension not built: {_IMPORT_ERROR}")
class TestUnsupportedShellsAreRejected(unittest.TestCase):
    """An unrenderable shell must block, not be skipped.

    Skipping left the shell's coefficient slots unreserved, so every later
    shell read the wrong slice of the MO coefficient vector.
    """

    def test_s_through_g_are_accepted(self):
        eng = BasisSetEngine([_shell(i) for i in range(5)])
        self.assertEqual(eng.n_basis, 1 + 3 + 5 + 7 + 9)

    def test_an_h_shell_is_rejected(self):
        with self.assertRaises(UnsupportedBasisError):
            BasisSetEngine([_shell(0), _shell(5)])

    def test_the_message_names_the_shell_letter(self):
        with self.assertRaises(UnsupportedBasisError) as ctx:
            BasisSetEngine([_shell(5)])
        self.assertIn("H shell", str(ctx.exception))

    def test_an_out_of_range_type_is_rejected(self):
        with self.assertRaises(UnsupportedBasisError):
            BasisSetEngine([_shell(99)])


if __name__ == "__main__":
    unittest.main()

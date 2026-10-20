"""
tests/test_parser_split.py
Asserts the parser.py -> mixin-module decomposition holds: every previously
public parse_* method still resolves on OrcaParser, each mixin owns exactly
its assigned methods, and no method name is duplicated across mixins.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from _parser_loader import load_standalone_parser  # noqa: E402

_mod = load_standalone_parser("orca_parser_split_check")
OrcaParser = _mod.OrcaParser
ParseCancelled = _mod.ParseCancelled

from _parser_loader import _PKG_DIR  # noqa: E402

sys.path.insert(0, _PKG_DIR)
import parser_structure as _pstruct  # noqa: E402
import parser_electronic as _pelec  # noqa: E402
import parser_properties as _pprops  # noqa: E402
import parser_spectra as _pspec  # noqa: E402

_EXPECTED = {
    "_StructureParsingMixin": (
        _pstruct._StructureParsingMixin,
        [
            "_parse_xyz_row",
            "parse_xyz_content",
            "parse_basic",
            "parse_termination_status",
            "parse_trajectory",
            "parse_scan",
            "parse_gradient",
            "parse_gradients",
            "parse_scan_results_table",
        ],
    ),
    "_ElectronicParsingMixin": (
        _pelec._ElectronicParsingMixin,
        [
            "parse_mo_coeffs",
            "parse_orbital_energies",
            "parse_basis_set",
            "parse_scf_trace",
        ],
    ),
    "_PropertyParsingMixin": (
        _pprops._PropertyParsingMixin,
        [
            "parse_dipole",
            "parse_spin_contamination",
            "parse_dispersion",
            "parse_energy_components",
            "parse_mayer_bond_orders",
            "parse_nbo_orbitals",
            "parse_nbo_perturbation",
            "parse_nbo_hybrids",
            "parse_charges",
        ],
    ),
    "_SpectraParsingMixin": (
        _pspec._SpectraParsingMixin,
        [
            "parse_nmr",
            "parse_tddft",
            "parse_thermal",
            "parse_frequencies",
            "_find_intensity_column",
            "_parse_dipole_derivative",
        ],
    ),
}


class TestEveryMethodStillResolves(unittest.TestCase):
    def test_all_prior_public_names_resolve_on_orca_parser(self):
        for _clsname, (_cls, names) in _EXPECTED.items():
            for name in names:
                with self.subTest(name=name):
                    self.assertTrue(
                        hasattr(OrcaParser, name),
                        f"OrcaParser lost attribute {name!r}",
                    )

    def test_kept_top_level_names_still_present(self):
        self.assertTrue(hasattr(_mod, "ParseCancelled"))
        self.assertTrue(hasattr(_mod, "OrcaParser"))
        self.assertTrue(issubclass(ParseCancelled, Exception))


class TestMixinsOwnExactlyTheirMethods(unittest.TestCase):
    def test_each_mixin_defines_exactly_its_methods(self):
        for clsname, (cls, names) in _EXPECTED.items():
            with self.subTest(cls=clsname):
                own_methods = {
                    k
                    for k, v in vars(cls).items()
                    if not k.startswith("__") and callable(getattr(cls, k))
                }
                self.assertEqual(own_methods, set(names))

    def test_no_method_defined_in_more_than_one_mixin(self):
        seen = {}
        for clsname, (_cls, names) in _EXPECTED.items():
            for name in names:
                self.assertNotIn(
                    name,
                    seen,
                    f"{name!r} defined in both {seen.get(name)!r} and {clsname!r}",
                )
                seen[name] = clsname


class TestMroOrder(unittest.TestCase):
    def test_orca_parser_inherits_in_the_documented_order(self):
        mro_names = [c.__name__ for c in OrcaParser.__mro__]
        expected_order = [
            "OrcaParser",
            "_StructureParsingMixin",
            "_ElectronicParsingMixin",
            "_PropertyParsingMixin",
            "_SpectraParsingMixin",
        ]
        for name in expected_order:
            self.assertIn(name, mro_names)
        indices = [mro_names.index(n) for n in expected_order]
        self.assertEqual(indices, sorted(indices))


if __name__ == "__main__":
    unittest.main()

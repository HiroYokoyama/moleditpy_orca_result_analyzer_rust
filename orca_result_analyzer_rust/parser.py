"""
parser.py — OrcaParser backed by the Rust extension `orca_parser_rs`.

The Rust extension must be compiled before use:
    python build.py

If the extension is absent a RuntimeError is raised with build instructions.
The public API is identical to the original pure-Python parser so all GUI
and analysis modules work without modification.
"""

import logging
import re
import os

from .parser_structure import _StructureParsingMixin
from .parser_electronic import _ElectronicParsingMixin
from .parser_properties import _PropertyParsingMixin
from .parser_spectra import _SpectraParsingMixin

#: 1 atomic unit of electric dipole moment (e*a0) in Debye.
#: Mirrors the Python plugin's module constant — part of the drop-in surface,
#: and the crate keeps its own copy in sync.
AU_TO_DEBYE = 2.541746473

#: Modes below this magnitude are translations/rotations or numerical noise,
#: not genuine imaginary frequencies (cm^-1).
IMAGINARY_FREQ_THRESHOLD = 10.0

try:
    from . import orca_parser_rs as _rs
    _RUST_AVAILABLE = True
except (ImportError, SystemError):
    # Fallback: direct file loading (e.g. tests load parser.py without package context)
    try:
        import importlib.util as _ilu
        _pkg_dir = os.path.dirname(os.path.abspath(__file__))
        _exts = [f for f in os.listdir(_pkg_dir)
                 if "orca_parser_rs" in f and f.endswith((".pyd", ".so"))]
        if _exts:
            _spec = _ilu.spec_from_file_location("orca_parser_rs",
                                                  os.path.join(_pkg_dir, _exts[0]))
            _rs = _ilu.module_from_spec(_spec)
            _spec.loader.exec_module(_rs)
            _RUST_AVAILABLE = True
        else:
            _rs = None
            _RUST_AVAILABLE = False
    except (ImportError, SystemError, OSError, AttributeError, TypeError, ValueError):
        _rs = None
        _RUST_AVAILABLE = False


def _rust_required():
    raise RuntimeError(
        "orca_parser_rs extension not found.\n"
        "Build it first:\n"
        "    pip install maturin\n"
        "    python build.py\n"
        "Then restart the plugin."
    )


class ParseCancelled(Exception):
    """Raised by the loading UI before or after the native parsing pass."""


class OrcaParser(_StructureParsingMixin, _ElectronicParsingMixin, _PropertyParsingMixin, _SpectraParsingMixin):
    """
    Parser for ORCA quantum chemistry output files.
    Parsing is performed by the compiled Rust extension (orca_parser_rs).
    """

    PARSE_STEPS = (
        ("parse_basic", "Reading job information"),
        ("parse_gradients", "Reading gradients"),
        ("parse_trajectory", "Reading trajectory"),
        ("parse_frequencies", "Reading vibrational frequencies"),
        ("parse_thermal", "Reading thermochemistry"),
        ("parse_mo_coeffs", "Reading MO coefficients"),
        ("parse_orbital_energies", "Reading orbital energies"),
        ("parse_charges", "Reading atomic charges"),
        ("parse_mayer_bond_orders", "Reading Mayer bond orders"),
        ("parse_nbo_orbitals", "Reading NBO orbitals"),
        ("parse_nbo_hybrids", "Reading NBO hybrids"),
        ("parse_nbo_perturbation", "Reading NBO perturbation analysis"),
        ("parse_dipole", "Reading dipole moment"),
        ("parse_spin_contamination", "Reading spin contamination"),
        ("parse_dispersion", "Reading dispersion correction"),
        ("parse_energy_components", "Reading energy components"),
        ("parse_tddft", "Reading TD-DFT excitations"),
        ("parse_nmr", "Reading NMR data"),
        ("parse_basis_set", "Reading basis set"),
        ("parse_scf_trace", "Reading SCF convergence"),
        ("parse_scan_results_table", "Reading scan results"),
    )

    def __init__(self):
        self.filename = ""
        self.raw_content = ""
        self.lines = []
        self._reset_data()

    def _reset_data(self):
        self.data = {
            "scf_traces": [],
            "converged": False,
            "scf_energy": None,
            "atoms": [],
            "coords": [],
            "charge": 0,
            "mult": 1,
            "frequencies": [],
            "excitation_energies": [],
            "tddft": [],
            "dipole": None,
            "dipoles": None,
            "nmr_shielding": [],
            "nmr_couplings": [],
            "charges": {},
            "version": None,
            "orbital_energies": [],
            "mos": [],
            "gradients": [],
            "all_gradients": [],
            "thermal": {},
            "scan_steps": [],
            "all_steps": [],
            "termination_status": "Running",
        }

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load_from_memory(self, content, filename="", progress=None):
        self._reset_data()
        self.filename = filename
        self.raw_content = content
        lines = content.splitlines()
        start = 0
        version_line = None
        for i, line in enumerate(lines):
            if re.match(r"^\s*Program Version\s+\d", line):
                start = i
                version_line = line
            elif re.match(r"^\s*\*?\s*JOB NUMBER\s+\d+\s*\*?\s*$", line):
                start = i
        self.lines = lines[start:]
        if version_line and not any("Program Version" in line for line in self.lines):
            self.lines.insert(0, version_line)
        self.parse_all(progress=progress)

    def parse_all(self, progress=None):
        """Native bulk parsing plus the shared NBO/spin compatibility sections.

        The native pass runs once. Cancellation is checked before it and during
        the following Python compatibility steps, not inside the native pass.
        """
        if not _RUST_AVAILABLE:
            _rust_required()
        total = len(self.PARSE_STEPS)
        if progress is not None:
            progress(0, total, self.PARSE_STEPS[0][1])
        content = "\n".join(self.lines) if self.lines else self.raw_content
        self.data = _rs.parse_all(content)
        shared = {"parse_spin_contamination", "parse_nbo_orbitals", "parse_nbo_hybrids", "parse_nbo_perturbation"}
        for done, (method, label) in enumerate(self.PARSE_STEPS):
            if done and progress is not None:
                progress(done, total, label)
            if method in shared:
                getattr(self, method)()
        if progress is not None:
            progress(total, total, "Finished")

    def parse_xyz_content(self, content):
        """Parse multi-frame XYZ content (e.g. NEB trajectory files)."""
        if not _RUST_AVAILABLE:
            _rust_required()
        return _rs.parse_xyz_content(content)

    # ------------------------------------------------------------------
    # Stub aliases kept for API compatibility with analysis modules
    # that call individual parse_* methods directly.
    # Each simply re-runs the full Rust parse pass.
    # ------------------------------------------------------------------

    def _reparse(self):
        self.parse_all()

    parse_basic            = _reparse
    parse_gradients        = _reparse
    parse_gradient         = _reparse
    parse_trajectory       = _reparse
    parse_scan             = _reparse
    parse_frequencies      = _reparse
    parse_thermal          = _reparse
    parse_mo_coeffs        = _reparse
    parse_orbital_energies = _reparse
    parse_charges          = _reparse
    parse_dipole           = _reparse
    parse_tddft            = _reparse
    parse_nmr              = _reparse
    parse_basis_set        = _ElectronicParsingMixin.parse_basis_set
    parse_scf_trace        = _reparse
    parse_scan_results_table = _reparse
    parse_mayer_bond_orders = _reparse
    parse_dispersion = _reparse
    parse_energy_components = _reparse

    parse_termination_status = _reparse

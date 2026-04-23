"""
parser.py — OrcaParser backed by the Rust extension `orca_parser_rs`.

The Rust extension must be compiled before use:
    python build.py

If the extension is absent a RuntimeError is raised with build instructions.
The public API is identical to the original pure-Python parser so all GUI
and analysis modules work without modification.
"""

import logging
import os

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
    except Exception:
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


class OrcaParser:
    """
    Parser for ORCA quantum chemistry output files.
    Parsing is performed by the compiled Rust extension (orca_parser_rs).
    """

    def __init__(self):
        self.filename = ""
        self.raw_content = ""
        self.lines = []
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
        }

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load_from_memory(self, content, filename=""):
        self.filename = filename
        self.raw_content = content
        self.lines = content.splitlines()
        self.parse_all()

    def parse_all(self):
        """Delegate all parsing to the Rust extension."""
        if not _RUST_AVAILABLE:
            _rust_required()
        try:
            self.data = _rs.parse_all(self.raw_content)
        except Exception as e:
            logging.error("[parser.py] Rust parse_all raised: %s", e)
            raise

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
        if self.raw_content:
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
    parse_basis_set        = _reparse
    parse_scf_trace        = _reparse
    parse_scan_results_table = _reparse

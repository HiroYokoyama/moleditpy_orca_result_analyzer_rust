# ORCA Result Analyzer — Test Suite

## Overview

All tests run headlessly without a display server and without requiring PyQt6,
PySide6, or a GPU.  PyQt6 and `matplotlib.backends.backend_qtagg` are stubbed
at module level in tests that need them.  `parser.py` and `utils.py` are
loaded directly via `importlib.util` with no stubs at all, since they have
no Qt or third-party dependencies.

## Running

```bash
# Run all tests
python -m pytest tests/ -v

# Run a single file
python -m pytest tests/test_parser.py -v

# Run a single test
python -m pytest tests/test_parser.py::TestParseBasic::test_scf_energy -v

# Run with coverage
python -m pytest tests/ --cov=orca_result_analyzer --cov-report=term-missing
```

## Test Files

### test_utils.py
**Pure utility functions — no stubs required**

| Test | What is tested |
|---|---|
| `TestGetDefaultExportPath` | Path generation from base path: CSV/PNG extensions, custom suffix, empty input, no-directory case, directory preservation |

---

### test_parser.py
**OrcaParser — pure Python, no stubs required**

`parser.py` imports only `re` and `logging`, so `OrcaParser` is loaded directly
without any stubs.  Each test class calls a single `parse_*` method in isolation
using minimal ORCA output snippets as input.

| Class | Methods tested |
|---|---|
| `TestParseBasic` | SCF energy (last occurrence wins), convergence detection (`SCF CONVERGED`, `OPTIMIZATION CONVERGED`, `HURRAY`), charge, multiplicity, version, Cartesian geometry (Angstrom vs A.U.), `is_scan`/`is_neb` flags, empty content |
| `TestParseScfTrace` | Single block, default "Initial" label, opt-cycle label, two-cycle → two traces, duplicate label → suffixed label (`Initial (2)`), empty content |
| `TestParseDipole` | Vector components, magnitude from explicit line, magnitude calculated as √(x²+y²+z²) when absent, multiple occurrences (last wins), no-dipole output |
| `TestParseCharges` | Mulliken count and values (colon format), Loewdin, Hirshfeld with spin, empty content |
| `TestParseFrequencies` | 6-mode count, zero/real/imaginary values, `cm-1` alternative unit, empty content |
| `TestParseNmr` | Shielding count and values, empty content, J-coupling table (unique pairs only) |
| `TestParseTddft` | State count, eV/nm energies, sorted by energy, transitions captured, empty content |
| `TestParseTrajectory` | Single opt cycle (energy + atoms + type), two cycles, scan step (type, scan_coord, step index), empty content |
| `TestParseGradients` | Colon format (count + values), no-colon format, NORM line excluded, `all_gradients` list, empty content |
| `TestParseXyzContent` | Multi-frame count and energies, atom coordinates, TS filtering, CI filtering, CI-NEB *not* filtered |

---

### test_traj_analysis.py
**TrajectoryResultDialog logic methods — Qt + matplotlib Qt backends stubbed**

`QDialog` and `FigureCanvasQTAgg` are stubbed as inheritable Python classes
(not `MagicMock()` instances) because `TrajectoryResultDialog` and `MplCanvas`
inherit from them.  `rdkit` and `PIL` are **not** stubbed — `matplotlib`
imports PIL internally and installing a fake PIL breaks `matplotlib.__init__`.
Both are guarded with `try/except ImportError` in the source.

Logic methods are tested by calling them as unbound methods on minimal fake
objects, so `TrajectoryResultDialog.__init__` is never invoked.

| Class | What is tested |
|---|---|
| `TestComputeScanPoints` | No scan IDs → returns all steps; groups by `scan_step_id` (last wins); sorted by id; `opt_cycle` preferred over other types; non-opt-cycle fallback; single step per group; empty list |
| `TestUpdateDisplayValues` | kJ/mol relative, kcal/mol relative, eV relative; absolute kJ/mol; unknown unit → factor 1.0 (Eh); minimum energy is zero in relative mode; single step → 0.0 |

**Key design note:** `QDialog = _BaseWidget` (a real stub class) is required —
`QDialog = MagicMock()` (an instance) cannot be used as a base class in Python
and would make `TrajectoryResultDialog` itself a MagicMock, losing all methods.

---

### test_init.py
**Plugin initialization contract — minimal Qt stubs**

`__init__.py` imports `QMessageBox` and `QFileDialog` at module level, so
minimal PyQt6 stubs are installed before loading.  Qt-dependent paths inside
`open_orca_file()` (dialog creation, `QApplication.processEvents`, etc.) are
not exercised here — they require a live Qt application.

| Class | What is tested |
|---|---|
| `TestMetadata` | `PLUGIN_NAME`, `PLUGIN_VERSION`, `PLUGIN_AUTHOR`, `PLUGIN_DESCRIPTION` constants present and non-empty |
| `TestInitialize` | `initialize(ctx)` registers a `.out` file opener with priority 100; registers exactly one drop handler with priority 100 |
| `TestDropHandler` | Non-.out extensions return False; `.log` returns False; non-existent `.out` path returns False; `.out` file without ORCA header content returns False; empty `.out` file returns False; case-insensitive extension |
| `TestInitializeIdempotent` | A fresh context after `initialize()` has exactly 1 opener and 1 handler |

**Not tested (requires Qt):** `open_orca_file()` dialog lifecycle,
`handle_drop()` True-return path (ORCA header detected), menu action callback.

## Mocking Strategy

Each test file installs its own stubs at **module level** to avoid
collection-order dependencies.

| Module | Approach |
|---|---|
| `parser.py` | Loaded directly — no stubs |
| `utils.py` | Loaded directly — no stubs |
| `traj_analysis.py` | Qt and matplotlib Qt backend stubbed; `QDialog`/`FigureCanvasQTAgg` as real inheritable classes |
| `__init__.py` | Minimal PyQt6 (`QMessageBox`, `QFileDialog`) stubbed as `MagicMock()` |

### test_parser_extended.py
**Extended OrcaParser coverage — no stubs required**

| Class | What is tested |
|---|---|
| `TestParseThermal` | electronic energy, ZPE, enthalpy, Gibbs, entropy, temperature; minimal block; empty content |
| `TestParseOrbitalEnergies` | restricted RHF (3 MOs), UHF alpha+beta spins, energy_eh/energy_ev, occupation, type (occupied/virtual), backward-compat `mos` list; empty content |
| `TestParseFrequenciesIR` | IR intensities stored by mode index (T² column); multiple modes; zero for non-IR modes; Raman activity stored by mode index |
| `TestParseChargesMayer` | Mayer QA/VA/BVA/FA values; fallback into Mulliken when no separate Mulliken block; empty content |
| `TestParseScanResultsTable` | builds scan_steps from "Actual Energy" summary; coord/energy mapping; maps coords into existing trajectory steps; empty content |
| `TestParseTrajectoryNeb` | NEB PATH SUMMARY → neb_image steps, energies, distances |
| `TestParseOptCycleConvergence` | convergence dict populated from GEOMETRY CONVERGENCE block; YES flags |
| `TestParseGradientsMultipleBlocks` | two gradient blocks in all_gradients; default = last block; first block accessible |

---

## Coverage Summary

| Module | Tests cover |
|---|---|
| `utils.py` | ~100% |
| `parser.py` | All `parse_*` methods including thermal, orbital energies, IR/Raman, Mayer charges, scan results, NEB path, opt convergence, multi-gradient blocks; ~80-85% of parser logic |
| `traj_analysis.py` | `compute_scan_points`, `update_display_values`, `recalc_energies`; Qt-dependent UI paths excluded |
| `__init__.py` | `initialize()`, `handle_drop()` False paths; `open_orca_file()` excluded |

Modules with Qt UI code (`gui.py`, `scf_analysis.py`, `freq_analysis.py`,
`mo_analysis.py`, `force_analysis.py`, etc.) require a live display to
instantiate and are not covered in this headless suite.

# ORCA Result Analyzer — Test Suite

## Overview

All tests run headlessly without a display server and without requiring PyQt6,
PySide6, or a GPU.  PyQt6 and `matplotlib.backends.backend_qtagg` are stubbed
at module level in tests that need them.  `parser.py` and `utils.py` are
loaded directly via `importlib.util` with no stubs at all, since they have
no Qt or third-party dependencies.

**Current status: 382 tests, all passing. `parser.py` coverage: 79%.
Static analysis (ruff): 9 minor issues (7 × E741 ambiguous names in quantum
chemistry code, 2 × E402 intentional imports-after-metadata in `__init__.py`).
All E722 bare-except violations resolved.**

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

### test_parser_samples.py
**OrcaParser integration tests — real ORCA output files in `sample_outputs/`**

Each test class loads a real `.out` file and asserts concrete parsed values.
Covers both ORCA 5 and ORCA 6 output for every sample type.

#### `benzene-opt.out` / `benzene-opt_5.out` — Geometry optimization (2 cycles)

| Test assertion | ORCA 6 (`benzene-opt.out`) | ORCA 5 (`benzene-opt_5.out`) |
|---|---|---|
| version | `6.1.1` | `5.0.4` |
| converged | `True` | `True` |
| SCF energy (Eh) | `−231.930734972084` | `−231.930718933873` |
| charge | `0` | — |
| multiplicity | `1` | — |
| atom count | `12` | `12` |
| atom symbols | 6 C + 6 H | — |
| first atom coords (Å) | `(0.805933, −1.143600, −0.008261)` | `(0.805991, −1.143635, −0.008264)` |
| total scan steps | `3` | `3` |
| `opt_cycle` count | `2` | `2` |
| `opt_final` count | `1` | `1` |
| `opt_final` is last step | ✓ | ✓ |
| `opt_final` coords = `data["coords"]` | ✓ | ✓ |
| `opt_final` energy (Eh) | `−231.930734972084` | — |
| cycle 1 step number | `1` | — |
| cycle 2 step number | `2` | — |
| cycle 1 atom count | `12` | — |
| cycle 1 coords ≠ final coords | ✓ | — |
| dipole magnitude | `< 0.1` (symmetric) | — |
| Mulliken count | `12` | `12` |
| Mulliken sum | `~0.0` | — |

#### `acetone-opt.out` / `acetone-opt_5.out` — Geometry optimization (5 cycles)

| Test assertion | ORCA 6 (`acetone-opt.out`) | ORCA 5 (`acetone-opt_5.out`) |
|---|---|---|
| converged | `True` | `True` |
| atom count | `10` | `10` |
| atom symbols | 3 C + 1 O + 6 H | — |
| first atom coords (Å) | `(−1.274842, 0.247667, 0.079907)` | — |
| total scan steps | `6` | `6` |
| `opt_cycle` count | `5` | `5` |
| `opt_final` count | `1` | `1` |
| `opt_final` is last step | ✓ | ✓ |
| `opt_final` coords = `data["coords"]` | ✓ | ✓ |
| cycle 1 coords ≠ final coords | ✓ | — |
| Mulliken count | `10` | — |
| Mulliken sum | `~0.0` | — |
| dipole present | ✓ | ✓ |
| dipole magnitude (Debye) | `~2.787` (> 2.0) | `~2.787` |
| dipole z-component dominant | ✓ (C=O along z) | — |

#### `benzene-opt-ene.out` / `benzene-opt-ene_5.out` — Single-point energy

| Test assertion | ORCA 6 (`benzene-opt-ene.out`) | ORCA 5 (`benzene-opt-ene_5.out`) |
|---|---|---|
| SCF energy (Eh) | `−231.930784926115` | `−231.930718955827` |
| converged | `True` | — |
| `opt_cycle` count | `0` | `0` |
| `opt_final` count | `0` | `0` |
| atom count | `12` | `12` |
| dipole magnitude | `< 0.1` (symmetric) | — |
| Mulliken count | `12` | `12` |

#### `benzene-opt-nmr.out` / `benzene-opt-nmr_5.out` — NMR shielding

| Test assertion | ORCA 6 (`benzene-opt-nmr.out`) | ORCA 5 (`benzene-opt-nmr_5.out`) |
|---|---|---|
| SCF energy (Eh) | `−231.930784926115` | `−231.930718955827` |
| NMR shielding count | `12` (6 C + 6 H) | `12` |
| NMR shielding present | ✓ | ✓ |
| C shieldings all positive | ✓ | — |
| H shieldings all positive | ✓ | — |
| `opt_cycle` count | `0` | `0` |

#### `benzene-opt-vex.out` / `benzene-opt-vex_5.out` — TD-DFT excited states

| Test assertion | ORCA 6 (`benzene-opt-vex.out`) | ORCA 5 (`benzene-opt-vex_5.out`) |
|---|---|---|
| SCF energy (Eh) | `−231.728888949460` | `−231.728828351432` |
| TDDFT state count | `5` | `5` |
| first state energy (eV) | `5.494` | `5.494` |
| states sorted ascending | ✓ | ✓ |
| all states have transitions | ✓ | ✓ |
| state energy range | `5.0 – 8.5 eV` | — |
| `osc_len` field present | ✓ (length gauge) | ✓ |
| `osc_vel` field present | ✓ (velocity gauge) | ✓ |
| `energy_nm` > 0 | ✓ | ✓ |
| `energy_cm` > 0 | ✓ | ✓ |
| eV × nm ≈ 1239.84 | ✓ | ✓ |
| `osc` ≥ 0 | ✓ | ✓ |
| CD fields: `rot_len`, `rot_vel`, `rotatory_strength` | ✓ (benzene is achiral → ≈0) | ✓ |
| ORCA 5 vs 6 first-state eV within 0.05 | — | ✓ |
| `opt_cycle` count | `0` | `0` |
| atom count | `12` | `12` |
| Mulliken count | `12` | — |

#### Cross-version consistency (`TestVersionConsistency`)

| Test assertion | Pair compared |
|---|---|
| same atom count | `benzene-opt` ORCA 5 vs 6 |
| same scan step count | `benzene-opt` ORCA 5 vs 6 |
| same step type sequence | `benzene-opt` ORCA 5 vs 6 |
| same atom count | `benzene-opt-ene` ORCA 5 vs 6 |
| same NMR shielding count | `benzene-opt-nmr` ORCA 5 vs 6 |
| same TDDFT state count | `benzene-opt-vex` ORCA 5 vs 6 |
| same scan step count | `acetone-opt` ORCA 5 vs 6 |
| same atom count | `acetone-opt` ORCA 5 vs 6 |

**Key behaviour verified:**
- `data["atoms"]` / `data["coords"]` always reflect the FINAL ENERGY EVALUATION geometry (last `CARTESIAN COORDINATES (ANGSTROEM)` in file), not a cycle snapshot.
- `opt_final` step is appended after all `opt_cycle` steps; its coords equal `data["coords"]`.
- Dipole magnitude is stored in **Debye** (read from the `Magnitude (Debye)` line). The parser scans up to 6 lines forward from `Total Dipole Moment` to find both magnitude lines, preferring Debye. Fallback to a.u. (√(x²+y²+z²)) only when neither magnitude line is found.

#### `benzene-opt-nmr.out` / `benzene-opt-nmr_5.out` — NMR shielding (expanded)

| Test assertion | ORCA 6 | ORCA 5 |
|---|---|---|
| NMR shielding count | `12` (6 C + 6 H) | `12` |
| keys: `atom_idx`, `atom_sym`, `shielding` | ✓ | ✓ |
| atom indices sequential | ✓ | ✓ |
| C shieldings in 10–200 ppm | ✓ | ✓ (6 values) |
| H shieldings in 10–50 ppm | ✓ | ✓ (6 values) |
| `Mayer` charges present | ✓ | ✓ |
| Mayer count = 12 | ✓ | ✓ |
| Mayer sum ≈ 0.0 | ✓ | ✓ |
| ORCA 5 vs 6 mean C shielding within 5 ppm | — | ✓ |

#### `benzene-opt-ene.out` / `benzene-opt-ene_5.out` — Basis set + orbital energies + SCF trace

| Test assertion | ORCA 6 | ORCA 5 |
|---|---|---|
| `basis_set_shells` count > 30 | ✓ | ✓ |
| shells have `atom_idx`, `l`, `exps`, `coeffs` | ✓ | ✓ |
| all `l` values ≥ 0 | ✓ | ✓ |
| all exponents > 0 | ✓ | ✓ |
| all `atom_idx` in range | ✓ | ✓ |
| `orbital_energies` count > 20 | ✓ | ✓ |
| keys: `index`, `occupation`, `energy_eh`, `energy_ev`, `spin`, `type` | ✓ | ✓ |
| occupied and virtual orbitals present | ✓ | ✓ |
| HOMO energy < 0 | ✓ | ✓ |
| LUMO energy > HOMO | ✓ | ✓ |
| spin = `restricted` | ✓ | ✓ |
| `mos` list identical to `orbital_energies` | ✓ | ✓ |
| `scf_traces` count ≥ 1 (ene) | ✓ | ✓ |
| `scf_traces` count ≥ 2 (opt) | ✓ (`benzene-opt`) | ✓ (`benzene-opt_5`) |
| trace keys: `step`, `iterations` | ✓ | ✓ |
| all iteration energies < 0 | ✓ | ✓ |
| opt trace labels contain cycle info | ✓ | ✓ |

---

### test_cube_and_mo.py
**Cube file parsing + MO coefficients + orbital energies — numpy required, no Qt/PyVista**

Uses `benzene-opt-ene_cubes/benzene-opt-ene_MO_21.cube` and `benzene-opt-ene.out`.

| Class | What is tested |
|---|---|
| `TestParseCube` | `CubeVisualizer._parse_cube` logic (via standalone reimplementation): n_atoms=12, dims=(40,40,40), origin (Bohr), x/y/z step vectors, `is_angstrom=False` (nx>0 → Bohr), 64000 data points, non-zero values, both +/− phases present |
| `TestParseOrbitalEnergies` | 32 orbitals (21 occupied + 11 virtual); orbital 0 energy in Eh and eV; HOMO (idx 20) occupied, LUMO (idx 21) virtual; HOMO–LUMO gap positive; `energy` alias equals `energy_eh`; `mos` list same length as `orbital_energies` |
| `TestParseMoCoeffs` | 114 MO entries (all basis functions); all keys end in `_restricted`; MO 0 energy/occ/spin/114 coefficients; first coeff atom/orbital/value; HOMO (MO 20) occupied; LUMO (MO 21) virtual; every MO has ≥1 coeff; all coeff entries have atom_idx/sym/orb/coeff fields |

**Notes:**
- `_parse_cube` is tested without PyVista by replicating its file-reading logic in the test — `_build_grid` (PyVista) is not exercised.
- `benzene-opt-ene.out` uses the merged basis-function label format (`0C  1s  …`) which the parser handles via the `^(\d+)([A-Za-z]+)$` regex branch.
- 114 MO coefficient entries reflect the full def2-SVP basis for C₆H₆ (6×C≈14 + 6×H≈5 functions).

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
| `TestComputeScanPoints` | No scan IDs → returns all steps; groups by `scan_step_id` (last wins); sorted by id; `opt_final` preferred over `opt_cycle`; `opt_cycle` preferred over other types; non-opt-cycle fallback; single step per group; empty list |
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

---

## Sample Output Files

All real ORCA output files are in `tests/sample_outputs/`.

| File | ORCA | Type | Atoms | Key data |
|---|---|---|---|---|
| `benzene-opt.out` | 6.1.1 | Geometry opt (2 cycles) | C₆H₆ (12) | Energy, frequencies, dipole, Mulliken |
| `benzene-opt_5.out` | 5.0.4 | Geometry opt (2 cycles) | C₆H₆ (12) | Same |
| `benzene-opt-ene.out` | 6.1.1 | Single-point energy | C₆H₆ (12) | Energy, orbital energies (32), MO coeffs (114), dipole, Mulliken |
| `benzene-opt-ene_5.out` | 5.0.4 | Single-point energy | C₆H₆ (12) | Energy, dipole, Mulliken |
| `benzene-opt-nmr.out` | 6.1.1 | NMR properties | C₆H₆ (12) | 12 NMR shieldings |
| `benzene-opt-nmr_5.out` | 5.0.4 | NMR properties | C₆H₆ (12) | 12 NMR shieldings |
| `benzene-opt-vex.out` | 6.1.1 | TD-DFT excited states | C₆H₆ (12) | 5 states (5.49–7.72 eV) |
| `benzene-opt-vex_5.out` | 5.0.4 | TD-DFT excited states | C₆H₆ (12) | 5 states |
| `acetone-opt.out` | 6.x | Geometry opt (5 cycles) | C₃H₆O (10) | Energy, dipole (~2.787 Debye), Mulliken |
| `acetone-opt_5.out` | 5.x | Geometry opt (5 cycles) | C₃H₆O (10) | Energy, dipole (~2.787 Debye), Mulliken |
| `benzene-opt-ene_cubes/benzene-opt-ene_MO_21.cube` | 6.1.1 | MO cube (HOMO−1) | — | 40×40×40 grid, Bohr, 64000 points |

---

## Mocking Strategy

Each test file installs its own stubs at **module level** to avoid
collection-order dependencies.

| Module | Approach |
|---|---|
| `parser.py` | Loaded directly — no stubs |
| `utils.py` | Loaded directly — no stubs |
| `vis.py` (_parse_cube only) | Cube reading logic reimplemented in test — no PyVista needed |
| `traj_analysis.py` | Qt and matplotlib Qt backend stubbed; `QDialog`/`FigureCanvasQTAgg` as real inheritable classes |
| `__init__.py` | Minimal PyQt6 (`QMessageBox`, `QFileDialog`) stubbed as `MagicMock()` |

## Coverage Summary

> Measured via `coverage run --source=orca_result_analyzer -m unittest discover`
> across all 382 tests.

| Module | Coverage | Notes |
|---|---|---|
| `utils.py` | **100%** | All branches hit |
| `parser.py` | **79%** | All `parse_*` methods exercised including relaxed surface scans; remaining 21% is NEB image step parsing and NBO/FMO charges |
| `traj_analysis.py` | **~12%** | `compute_scan_points`, `update_display_values`, `recalc_energies`; Qt-dependent UI paths excluded |
| `__init__.py` | **~20%** | `initialize()`, `handle_drop()` False paths; `open_orca_file()` excluded |
| `vis.py` | *partial* | `_parse_cube` file-reading path only; `_build_grid`/`show_iso` require PyVista + display |
| All Qt UI modules | **0%** | `gui.py`, `*_analysis.py`, etc. require a live display to instantiate |

### parser.py — Known coverage gaps (require new sample files)

| Gap | Lines | Why not covered |
|---|---|---|
| `parse_xyz_content` | 36–151 | No NEB `.xyz` multi-frame file in samples |
| Opt cycle alt energy formats | 483–507 | Require non-standard ORCA config |
| NBO / FMO charges | 1110–1282 | Sample files do not include NBO analysis |
| TDDFT short-table format (pattern B) | 1615–1628 | Benzene vex uses arrow format only |
| CD spectrum (rot_vel path) | 1647–1654 | Covered only partially (benzene CD ≈0) |
| `parse_scan_results_table` | 2285–2349 | No scan output in samples |

## Static Analysis (ruff)

| Status | Count | Code | Issue |
|---|---|---|---|
| ✅ Fixed | 0 | E722 | Bare `except` — all resolved |
| ⚠️ Remaining | 7 | E741 | Ambiguous variable name `l` / `I` in quantum chemistry math (angular momentum) |
| ⚠️ Remaining | 2 | E402 | Imports after plugin metadata constants in `__init__.py` (intentional) |


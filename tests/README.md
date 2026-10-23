# ORCA Result Analyzer — Test Suite

~1780 tests across 53 files, covering the ORCA parser, the pure analysis
logic, and every analysis dialog driven headlessly. No Qt binding, display, or
ORCA installation is required to run them.

## Running the tests

```bash
python run_tests.py                        # full suite, exactly as CI runs it
python run_tests.py -k parser              # any pytest args are forwarded
python run_tests.py tests/test_parser.py -v
python run_tests.py tests/test_nmr_dialog.py::TestRecalc
```

**Use `run_tests.py` rather than calling `pytest` directly.** It sets
`PYTEST_DISABLE_PLUGIN_AUTOLOAD=1` so local runs match CI, which installs no Qt
binding. An installed `pytest-qt` would otherwise import a real PySide6 during
collection and defeat the stubs the tests rely on — producing failures and
segfaults that do not reproduce in CI.

Coverage (autoload is off, so the plugin must be named explicitly):

```bash
python run_tests.py -p pytest_cov --cov=orca_result_analyzer --cov-report=term-missing
```

The target is to keep coverage above 80%. CI reports it on every run but does
not fail the build on it.

Treat the number with some care: this is largely a Qt codebase, and a
constructor that runs without raising counts as covered. The parsing, the
chemical-shift and peak arithmetic, the exporters and the settings round-trips
are where the tests actually pin behaviour.

### Python version

CI runs 3.11–3.13. The source uses 3.10+ syntax (`str | None`, `list[str]`), so
on an older interpreter `tests/test_api.py` fails at collection. On 3.9,
`--ignore=tests/test_api.py` lets the rest of the suite run, but treat CI as
authoritative.

## How the suite is organised

Tests are grouped by what they exercise rather than one-to-one with source
files. Per-file test counts are deliberately not listed here — they go stale
within a commit; run `--collect-only` for the current picture.

| Area | Files | What they cover |
|---|---|---|
| **Parser — units** | `test_parser.py`, `test_parser_extended.py`, `test_parser_charges.py`, `test_parser_xyz_regressions.py` | Individual `parse_*` methods against synthetic ORCA text, including the NBO fallback and FMO population blocks |
| **Parser — real output** | `test_parser_samples.py` | The same parsers against real ORCA 5 and 6 files, plus cross-version consistency |
| **Cube / MO engine** | `test_cube_and_mo.py`, `test_mo_engine.py`, `test_mo_overlap.py` | Cube parsing (including MO cubes with a negative atom count), basis-set evaluation |
| **Analysis dialogs** | `test_*_dialog.py`, `test_*_analysis.py`, `test_small_dialogs.py`, `test_dialog_extras.py` | Each `QDialog` driven through `gui_harness` — construction, controls, exports, settings |
| **Main window** | `test_gui_dialog.py`, `test_gui_interaction.py`, `test_gui_launchers.py`, `test_gui_load_file.py`, `test_about_menu.py`, `test_gui_reset_on_load.py` | Button enablement, file loading, 3D picking, drag-and-drop, analysis-window lifecycle |
| **Plugin contract** | `test_init.py`, `test_plugin_integration.py`, `test_api.py` | Entry points, the host `PluginContext` API, and drift against the real main app |
| **Shared helpers** | `test_utils.py` | Pure utility functions |

## The Qt harness

The dialogs need Qt and a live 3D viewer at import time, neither of which
exists in CI. `tests/gui_harness.py` bridges that gap.

`load_isolated("module_name")` loads one source module against subclassable Qt
stand-ins that keep **real state** — spin boxes, combos, buttons, check boxes,
sliders, line edits, trees, tables and lists all remember what you set — then
restores the shared `sys.modules` so other test modules are unaffected. It also
swaps in the genuine numpy and matplotlib for the duration, since other test
modules leave stubs of those behind.

`tests/test_nmr_dialog.py` is a representative example.

Traps worth knowing before you add a test are documented in
[`CLAUDE.md`](../CLAUDE.md) — chiefly that the widget stubs are permissive (a
mistyped method name passes silently), that settings paths are derived from the
module directory and will otherwise write into the package source tree, and
that some Qt loop idioms never terminate against a MagicMock.

### Older stub style

Test modules predating the harness install their own partial PySide6 stubs at
**module level**, to avoid collection-order dependencies. Those still work and
have not been converted. If you touch one, note that any `sys.modules` stub for
a third-party package (`numpy`, `matplotlib`, `PIL`) must be a *fallback used
only when the real package is absent* — overwriting a working package breaks
every later test that does real numerics, because matplotlib resolves numpy and
its Agg backend from `sys.modules` at call time.

## Sample output files

Real ORCA output lives in `tests/sample_outputs/` (19 files). Files ending
`_5` are ORCA 5 counterparts of the same calculation, used for cross-version
consistency checks.

| File | ORCA | Type | Key data |
|---|---|---|---|
| `benzene-opt.out` / `_5` | 6.1.1 / 5.0.4 | Geometry opt (2 cycles), C₆H₆ | Energy, frequencies, dipole, Mulliken |
| `benzene-opt-ene.out` / `_5` | 6.1.1 / 5.0.4 | Single-point energy, C₆H₆ | Orbital energies (32), MO coeffs (114), SCF trace |
| `benzene-opt-nmr.out` / `_5` | 6.1.1 / 5.0.4 | NMR properties, C₆H₆ | 12 shieldings |
| `benzene-opt-vex.out` / `_5` | 6.1.1 / 5.0.4 | TD-DFT, C₆H₆ | 5 states (5.49–7.72 eV) |
| `acetone-opt.out` / `_5` | 6.x / 5.x | Geometry opt (5 cycles), C₃H₆O | Energy, dipole (~2.787 D), Mulliken |
| `ethane-scan.out` | 6.x | Relaxed surface scan | Scan steps and results table |
| `nbo-test.out` | 6.x | NBO analysis | Natural populations, hybrids, E(2) |
| `nmr-coupling-test.out` | 6.x | NMR with spin-spin coupling | Shieldings and J-couplings |
| `chelpg-test.inp.out`, `mbis-test.out`, `resp-test.out` | 6.x | Charge schemes | CHELPG / MBIS / RESP populations |
| `ccsd-test.out` | 6.x | `DLPNO-CCSD(T)` | Correlated energy components |
| `mp2-test.out` | 6.x | `MP2` | Correlated energy components |
| `o2-opt_orb.out` | 6.x | Open-shell `UKS` opt + freq, O₂ | Alpha/beta orbitals, ⟨S²⟩ contamination |
| `benzene-opt-ene_cubes/…_MO_21.cube` | 6.1.1 | MO cube (HOMO−1) | 40×40×40 grid, Bohr |

See `tests/sample_outputs/NOTICE` for the licensing position on these files —
they are ORCA output, excluded from this project's LGPL-3.0-or-later licence.

## Integration test strategy

The plugin is loaded by a host application, so the host's `PluginContext` API
is a contract that can drift. Two tiers guard it:

```
Stub mode (always runs)
  _StubContext mirrors the real PluginContext API.
  No main app required. Catches interface mismatches immediately.

Real-context mode (runs when the main app is available)
  Uses the actual PluginContext from python_molecular_editor.
  Catches drift invisible to the stub.
```

Real-context mode activates when the repos are siblings:

```
<parent>/
    moleditpy_orca_result_analyzer_plugin/   ← this plugin
    python_molecular_editor/                 ← main app
```

### CI jobs

| Job | Python | Main app | Tests |
|---|---|---|---|
| `test` | 3.11, 3.12, 3.13 | Cloned with `\|\| true` | Full suite with coverage; real-context tests run if the clone succeeded |
| `test-integration` | 3.11 | Cloned (hard failure if the clone fails) | `test_plugin_integration.py` + `test_api.py` including the real-context tier |

## Adding tests

- Match the fixture shapes `parser.py` actually produces, or the test passes
  against fiction. The shapes that have caught people out are listed in
  [`CLAUDE.md`](../CLAUDE.md#fixture-shapes).
- Parser behaviour that needs a new ORCA feature belongs in a new sample file
  under `sample_outputs/`, exercised from `test_parser_samples.py`.
- New dialog behaviour goes through `gui_harness.load_isolated`, not a
  hand-rolled PySide6 stub.

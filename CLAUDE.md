# CLAUDE.md

Guidance for Claude Code when working in **moleditpy_orca_result_analyzer_rust**.

## What this repo is

The **Rust-accelerated** build of the ORCA Result Analyzer plugin for MoleditPy.
It is the same plugin as the pure-Python
[`moleditpy_orca_result_analyzer_plugin`](../moleditpy_orca_result_analyzer_plugin),
except the parser and the MO/GTO grid engine are compiled Rust extensions for
speed. The UI/analysis code is kept **in parity** with the Python plugin.

## Source-of-truth split (READ BEFORE EDITING)

The Python plugin `moleditpy_orca_result_analyzer_plugin` is the **canonical
source** for all UI/analysis logic. This repo mirrors it, keeping only the
Rust engine local.

| Layer | Files | Rule |
|---|---|---|
| **Rust engine (local, do NOT sync from Python)** | `orca_result_analyzer_rust/parser.py`, `mo_engine.py`, `crates/orca_parser_rs/`, `crates/orca_mo_rs/` | Edit here; these wrap the compiled `.pyd`/`.so`. Their **public API is a drop-in** for the Python engine: `OrcaParser` (with the same `self.data` contract), `BasisSetEngine`, `CalcWorker`. |
| **Synced from the Python plugin** | every other `*.py`, including parser compatibility mixins (all `*_analysis.py`, `gui.py`, `__init__.py`, `mo_analysis.py`, `spectrum_widget.py`, `vis.py`, `utils.py`, `energy_diag.py`, …) | When the Python plugin changes these, copy them over. `mo_analysis.py` is synced but binds to the **local** `mo_engine` (same class names). |
| **Rust-only** | `logger.py`, `mo_energy_diag.py`, `build.py` | Local; may be unused by synced code but kept. |

**Version:** `orca_result_analyzer_rust/__init__.py` `PLUGIN_VERSION` mirrors the
Python plugin's version (same functionality → same number). The
`pyproject.toml` / `Cargo.toml` versions are the build-artifact versions and are
aligned with the plugin version at 5.0.3. Keep both crate manifests, Cargo.lock
and the project version consistent when releasing.
The supported native wheel range is Python 3.9–3.13.
NBO/spin sections and helper APIs use shared Python mixins; bulk parsing and MO
grid evaluation remain native. Native-pass cancellation takes effect after it returns.

When syncing, do **not** blindly overwrite `parser.py` / `mo_engine.py`, and
after copying dialog modules verify the seam names still resolve
(`utils.get_default_export_path`, `vis.CubeVisualizer`,
`energy_diag.EnergyDiagramDialog`, `mo_engine.BasisSetEngine`/`CalcWorker`).

### Before touching the shell definitions in `crates/orca_mo_rs/`

Read **`../moleditpy_orca_result_analyzer_plugin/docs/MO_CALCULATION.md`**.
The crate's `build_sph_defs` is the same derivation as the Python engine's,
and for the same reason: a mis-normalized basis function still renders as a
plausible orbital — right symmetry, right signs, lobes roughly in place — so
a change there **cannot be verified by looking at the picture**. It carried
hand-typed f and g prefactor tables that were wrong for a long time without
anyone noticing.

Verify numerically instead (`tests/test_mo_engine_shells.py` does all three):
cosine similarity 1.000000 against `scipy.special.sph_harm_y`, unit norm per
component, and the closed-form nodal angles (f(z³) 39.23°, g(z⁴)
30.56°/70.12°).

### Checking parser parity against the Python plugin

The two parsers must produce the same `self.data`. To check, load a sample
with both and diff the dicts key by key — that is how the missing
`mayer_bond_orders` / `termination_status` / `dispersion` /
`energy_components` blocks and the truncated SCF traces were found, none of
which any test was asserting. Known-benign differences: `all_steps` (Rust
populates it for optimisations), and the tie-break order of degenerate TD-DFT
states and scan-step sub-entries.

## Build (required before tests/use)

```bash
pip install maturin          # + Rust toolchain (https://rustup.rs)
python build.py              # dev build, copies .pyd/.so into orca_result_analyzer_rust/
python build.py --release    # optimised
python build.py --wheel      # wheels only
```
`build.py` builds both crates and produces a per-environment plugin zip in
`build/` named `orca_result_analyzer_rust-<PLUGIN_VERSION>-<platform_tag>.zip`.
The compiled `.pyd`/`.so` are **git-ignored — never commit them.**

## Testing

```bash
python run_tests.py          # CI-faithful (disables pytest plugin autoload)
python run_tests.py -k parser
```

Tests stub PyQt6/PySide6 at module level. Locally-installed `pytest-qt` imports
a real binding at collection and defeats those stubs (spurious failures /
segfaults), so **use `run_tests.py`** (it sets `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1`),
not raw `pytest`. CI installs no Qt binding, so it is unaffected.

Test split mirrors the source split: **keep the engine tests** here
(`test_parser*.py`, `test_cube_and_mo.py` — they validate the Rust engine's
output); the dialog/gui/init/utils tests are synced from the Python plugin (they
feed synthetic parser data, so they are engine-independent). When bringing a
Python test over, rewrite the package name `orca_result_analyzer` →
`orca_result_analyzer_rust`.

## CI (`.github/workflows/`)

| Workflow | Trigger | Does |
|---|---|---|
| `tests.yml` | push to `main`, PRs | build extensions + run the suite on an OS×Python matrix |
| `release.yml` | tag push (`v*` / `X.Y.Z`) | build the plugin zip for every OS×Python env and attach all to a GitHub Release |

## Relationship with the main app

Like all MoleditPy plugins, this validates against the `PluginContext` API in
`../python_molecular_editor/moleditpy/src/moleditpy/plugins/plugin_interface.py`.


## Synchronization and release checks

`SYNC_SOURCE.json` pins the canonical Python release commit. Update the shared
modules/tests and the `sync-check` workflow's checkout ref together, then run
`scripts/check_python_sync.py` against that exact revision. Keep native-specific
regressions in `tests/test_native_result_audit.py` or other local test modules;
do not diverge copied shared tests to hide native parsing differences.

Before release, build both native extensions and run the complete suite with
`run_tests.py`, including the canonical audit regressions and native malformed
input checks. Passing Python helper tests alone does not validate the native
bulk parser. The tag release workflow builds the full OS/Python ZIP matrix;
verify all assets and included license/provenance notices after it completes.
Python 3.9 Linux jobs use Ubuntu 24.04 because setup-python has no 3.9 build for
Ubuntu 26.04. Retain the other OS/Python combinations.

## Copied maintenance scripts

See [scripts/COMMON_TOOLS.md](scripts/COMMON_TOOLS.md) for canonical ownership, pinned versions and update instructions. Common script unit tests stay in the canonical repository; retain this repository's own integration checks.

The API checker engine is an exception to Python test-source synchronization: its canonical owner is `moleditpy-plugins/api-checker/`, pinned by `scripts/api-checker.lock.json`. Keep the consumer API integration class in parity, but do not restore shared engine unit tests from the Python checkout. The source-parity checker validates the engine copy by its pinned hash.

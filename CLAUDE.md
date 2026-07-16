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
| **Synced from the Python plugin** | every other `*.py` (all `*_analysis.py`, `gui.py`, `__init__.py`, `mo_analysis.py`, `spectrum_widget.py`, `vis.py`, `utils.py`, `energy_diag.py`, …) | When the Python plugin changes these, copy them over. `mo_analysis.py` is synced but binds to the **local** `mo_engine` (same class names). |
| **Rust-only** | `logger.py`, `mo_energy_diag.py`, `build.py` | Local; may be unused by synced code but kept. |

**Version:** `orca_result_analyzer_rust/__init__.py` `PLUGIN_VERSION` mirrors the
Python plugin's version (same functionality → same number). The
`pyproject.toml` / `Cargo.toml` versions are the build-artifact versions and are
managed separately.

When syncing, do **not** blindly overwrite `parser.py` / `mo_engine.py`, and
after copying dialog modules verify the seam names still resolve
(`utils.get_default_export_path`, `vis.CubeVisualizer`,
`energy_diag.EnergyDiagramDialog`, `mo_engine.BasisSetEngine`/`CalcWorker`).

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

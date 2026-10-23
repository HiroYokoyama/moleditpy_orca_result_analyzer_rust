# MoleditPy ORCA Result Analyzer — Rust Edition

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20726011.svg)](https://doi.org/10.5281/zenodo.20726011)
[![Tests](https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_rust/actions/workflows/tests.yml/badge.svg)](https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_rust/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/HiroYokoyama/moleditpy_orca_result_analyzer_rust/graph/badge.svg)](https://codecov.io/gh/HiroYokoyama/moleditpy_orca_result_analyzer_rust)
[![MoleditPy](https://img.shields.io/badge/MoleditPy->=5.0.0-3577F7)](https://github.com/HiroYokoyama/python_molecular_editor)
[![](https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86)](https://github.com/sponsors/HiroYokoyama)
[![GitHub tag](https://img.shields.io/github/v/tag/HiroYokoyama/moleditpy_orca_result_analyzer_rust?label=version)](https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_rust/tags)
[![GitHub Downloads](https://img.shields.io/github/downloads/HiroYokoyama/moleditpy_orca_result_analyzer_rust/total)](https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_rust/releases)

A comprehensive plugin for MoleditPy to analyze and visualize results from ORCA quantum chemistry calculations.

This edition synchronizes its UI, analysis features and complete test suite with
Python plugin **5.0.3**, commit `d7f79374829da98c2f4c84dc6811504d0f5fc164`.
The bulk output parser and MO grid evaluator remain compiled Rust extensions.
Shared Python mixins provide helper APIs and NBO/spin parsing. See
[PYTHON_SYNC.md](PYTHON_SYNC.md) for the file-by-file classification and
[SYNC_SOURCE.json](SYNC_SOURCE.json) for the synchronization source.

## Build and supported environments

Use **Moleditpy 5.x**, Python **3.9–3.13**, Rust and maturin. Python 3.14/3.15
are not declared supported by this edition's current PyO3 toolchain.
Runtime dependencies use PySide6, RDKit, NumPy, Matplotlib, Pillow and PyVista;
`nmrsim` is optional.

```bash
python -m pip install maturin
python build.py --release
```

The script builds both crates, copies the matching native extensions into
`orca_result_analyzer_rust/`, and writes an environment-specific plugin ZIP in
`build/`. Install the ZIP matching your OS, CPU architecture and Python version.
The release matrix covers Windows x86-64, Linux x86-64 and macOS ARM64.
Compiled files are gitignored and must not be committed.

File loading supports progress and cancellation between parsing sections. The
native bulk pass itself is synchronous; cancellation takes effect when that pass
returns. MO cube generation supports cancellation between grid chunks.



Repo: [https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_rust](https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_rust)

![Main UI](img/main.png)

## Features

### 1. SCF Trace
Real-time convergence visualization for SCF energy cycles.
- **Concatenated View**: View a single continuous plot of all SCF cycles throughout the calculation.
- **Interactive Tools**: Full zoom, pan, and save support via an integrated Matplotlib toolbar.

### 2. MO Analysis
- **Levels**: View orbital energies and occupancy with HOMO/LUMO identification.
- **Visualization**: Generate 3D Cubes (isosurfaces) with **Smooth Shading** and adjustable opacity.
- **Presets**: Save and manage visualization presets (colors, isovalues, styles).
- **Compare MOs**: Show up to four orbitals at once in a separate window, each with
  its own lobe colors, isovalue, opacity and style. Opens on the orbitals selected
  in the table — or on HOMO, LUMO, LUMO+1 and HOMO-1 when nothing is selected —
  and generates any missing cubes on demand. **Sync Iso** puts every slot on
  Orbital 1's contour level so the shapes are comparable. Also available from
  the orbital table's right-click menu.
- **Advanced Support**: Successfully handles S, P, D, F, and G shells (L=4).
- **Feedback**: Integrated warning system identifies missing ORCA output keywords required for cube generation.

![MO Analysis](img/mo.png)

![Compare Molecular Orbitals](img/mo_compare.png)

### 3. Optimization / Scan
Analyze **Geometry Optimizations** and **Relaxed Surface Scans**.
- **Interactive Graph**: Plot Energy vs. Step. Click points to update the 3D structure.
- **Display Modes**: Toggle between **Absolute** and **Relative** energy (kJ/mol, kcal/mol, eV, Eh).
- **Log Scale**: Supports log-scale visualization for relative energy changes.
- **Animation**: Play/Pause trajectory animations with adjustable FPS.
- **Export**: Save plots as images or export the full 3D animation as high-quality **GIFs**.

![Optimization](img/opt.png)

### 4. Forces
Analyze structural forces for the current structure or the entire trajectory.
- **Historical Gradients**: Capture and display force vectors for **every** step of an optimization or scan.
- **Visualization Controls**: **Auto Scale** feature automatically optimizes vector size for visibility, especially for small gradients.
- **Step-by-Step Navigation**: Precise control with `<` and `>` buttons or the trajectory slider.
- **Convergence Tracking**: Multi-line display of **RMS/MAX Gradient** and **RMS/MAX Step**, color-coded (**Green for YES**, **Red for NO**).
- **Modeless Threshold Graph**: View a non-blocking convergence threshold plot across all optimization steps. The graph features color-matched Y-axes for each metric and easy dropdown metric selection. (Shortcut: **Shift+Click** the "Forces" button in the main GUI to open this graph directly, bypassing the force viewer).
- **Data Table**: Full breakdown of gradients, force components, and magnitudes.

![threshold](img/conv_threshold.png)

### 5. Atomic Charges
- **Populations**: Mulliken, Loewdin, Hirshfeld, and **NBO** populations (if available).
- **3D Coloring**: Color atoms in the 3D viewer based on charge value or population type.

![Atomic Charges](img/charge.png)

### 6. Dipole Moment
- **Vector Visualization**: Display the total dipole moment vector magnitude and direction in the 3D viewer.

### 7. Frequencies
Visualize vibrational modes and spectra.
- **IR/Raman**: Stick and broadened spectra plots with interactive peak labels.
- **Visualization**: Animated vibrational modes with vector arrows.

![Frequencies](img/freq.png)

### 8. Thermochemistry
Detailed analysis of thermodynamic properties based on frequency calculations.
- **Broad Summary**: Electronic Energy, ZPE, Enthalpy (H), Gibbs Free Energy (G).
- **Detailed Corrections**: Optional breakdown of **vibrational, rotational, and translational** contributions to energy and entropy.

### 9. TDDFT
Analyze electronic excitations and absorption spectra.
- **Spectra**: **Absorption** and **CD** (Circular Dichroism) spectra with Gaussian broadening.
- **Controls**: Adjustable broadening (Sigma) and peak-stick overlays.

### 10. NMR
Advanced NMR chemical shielding validation and visualization.
- **Stick Spectrum**: Nucleus-specific stick spectra (1H, 13C, etc.) with experimental reference standards (TMS, CDCl3, DMSO-d6, etc.).
- **Multiplet Simulation**: Realistic J-coupling splitting patterns and first-order multiplicity calculations with adjustable Lorentzian line-widths.
- **Custom References**: Add and manage custom reference standards (delta_ref and sigma_ref).
- **Equivalent Atom Merging**: Manually merge equivalent atoms into single peaks. Merges are saved explicitly via the **Save Merges** button (you are prompted to save or discard unsaved merges when closing the dialog).
- **Interactive Sync**: Robust bidirectional synchronization—selecting peaks in the spectrum highlights atoms in 3D (and vice-versa). 3D atom labels can optionally show the chemical shift ("Show shift values" checkbox, off by default); merged atoms show both their original and merged values (δ orig → merged).

![NMR](img/nmr.png)

### 11. Bond Analysis
Comprehensive bonding and orbital analysis:
- **Mayer Bond Orders**: View and inspect parsed Mayer bond-order matrix values in a clean tabular view.
- **NBO Orbitals & Hybridization**: Full breakdown of Natural Bond Orbitals (NBO) list with hybridization composition (%s / %p / %d / %f per atom). Double-clicking an entry displays the raw NBO details.
- **E(2) Perturbation**: Tabulate donor-acceptor second-order perturbation energy \(E^{(2)}\) interaction values.
- **3D Viewer Integration**: Selecting rows in any bond analysis table automatically highlights the involved atom(s) in the 3D viewport.

### 12. Post-HF Energies & Properties
- **Post-HF Energy Components**: Separate panel displaying Reference Energy, Correlation Energy (MP2, CCSD), and Triples Correction (CCSD(T)) components.
- **Properties**: Parse and display physical parameters such as spin contamination expectation values (\(\langle S^2 \rangle\)) and dispersion corrections.

### 13. NICS (optional)
Hand the open output straight to the **ORCA NICS Analyzer** plugin via **Analysis > NICS Analysis...**.
- **Only when installed**: the menu entry appears only while that plugin is loaded, and hides again when it is not.
- **File hand-off**: the file already open in the analyzer is passed over, so no second file dialog is needed.

## Interface & Usability
- **MoleditPy v4 Ready**: Fully integrated with the new `PluginContext` architecture for stable window management and safe main-window interaction.
- **Standalone Launch**: Open the analyzer at any time via **Extensions > ORCA Result Analyzer** without selecting a file first.
- **Drag-and-Drop**: Drag a `.out` file directly onto the analyzer window to open it, or drag a folder to open the "Select from Directory" picker.
- **Logical Workflow**: Buttons grouped by task: Electronic -> Geometry -> Properties -> Spectroscopy.
- **Modeless Design**: All analysis windows are modeless, allowing side-by-side comparison and 3D viewer interaction.
- **Copyable Tables**: Support selecting and copying data (via `Ctrl+C`) from analysis tables to external spreadsheet editors.
- **Interactive Highlighting**: Highlight atoms in the 3D viewer (with VDW scaling) when selecting rows in the Properties or Bond Analysis tables.
- **Keyboard Shortcuts**: `Ctrl+O` (Open File), `Shift+Ctrl+O` (Select from Directory), `Ctrl+R` (Reload), `Ctrl+W` (Close Window).
- **Persistence**: Remembers your presets, NMR references, and explicitly saved merged peaks across sessions.

## Installation
Download the environment-specific ZIP from [Releases](../../releases), or build it locally.

Ensure the `orca_result_analyzer_rust` folder is placed in your MoleditPy plugins directory.

## Requirements
- **MoleditPy**: Version 5.0.0 or higher, below 6.0.0.
- **ORCA Output**: Reads `.out` files. Basis set info (`Print[P_Basis] 2`) is required for MO Cube generation.
- **Dependencies**: `rdkit`, `matplotlib`, `Pillow` (PIL), and `pyvista` (for 3D vectors). `nmrsim` is optional (for J-coupling simulation).

## Development

Run the test suite with `python run_tests.py` (it forwards any pytest
arguments). Use that wrapper rather than calling `pytest` directly — it
disables third-party plugin autoload so local runs match CI.

Tests run on Python 3.9–3.13 with no Qt binding installed; the dialogs are
exercised headlessly through `tests/gui_harness.py`. CI reports coverage on
every run, with a target of **above 80%**.

See [tests/README.md](tests/README.md) for how the suite is organised, the Qt
harness, and the sample ORCA files it parses; [the canonical contribution guide](https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_plugin/blob/main/CONTRIBUTING.md)
for the contribution workflow.

[docs/MO_CALCULATION.md](docs/MO_CALCULATION.md) documents how molecular
orbitals are evaluated — basis-function normalization, the spherical-harmonic
component tables, and how to verify a change to them numerically rather than
by eye. Read it before touching `mo_engine.py`.

## License Note on Sample Output Files

The test fixtures in `tests/sample_outputs/` are ORCA quantum chemistry output files (`.out`) generated by running [ORCA](https://orcaforum.kofo.mpg.de/) on publicly available molecular geometries (benzene, acetone).

ORCA is proprietary software with its own license. The source-code license does
not relicense these generated fixtures. See `tests/sample_outputs/NOTICE` for
provenance.

If you need to regenerate these files, you will need your own ORCA installation.

## Required ORCA Keywords

**For MO Cube Generation:**

*Note: Standard output is usually sufficient for Geometry and Energies, but Basis Set information is strictly required for generating cubes.*

```text
%output
  Print[P_Basis] 2  # Required for Basis Set parsing
  Print[P_Mos] 1    # Ensure MO coefficients are printed
end

```

**For NMR Simulation (J-Coupling):**

```text
! NMR

%eprnmr
  NUCLEI = ALL H {SHIFT, SSALL} # Required for J-coupling (nmrsim)
end
```

## License & Disclaimer

This edition retains the GPLv3 license in [LICENSE](LICENSE). Synced LGPL source
and third-party dependency notices are recorded in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Existing notices and historical
release licenses are retained.

## Rust edition tests and source boundaries

All canonical test files and fixtures are copied into `tests/`; native-specific
parser/MO tests are retained. Python tests remain identical after the package-name substitution. Spherical normalization, orthogonality,
nodal angles, MO coefficient counts, UHF spin selection, final-job isolation,
loading, GUI behavior and cube I/O remain covered.

```bash
python run_tests.py tests/ -p pytest_cov --cov=orca_result_analyzer_rust --cov-report=term-missing
```

Build the extensions first. The runner disables pytest plugin autoload so a local
pytest-qt installation cannot replace the headless Qt test stubs. Coverage reports
Python wrappers/UI only; the native Rust numerical assertions are separate tests.
The build-artifact and plugin versions are aligned at 5.0.3 for this synchronization.

This edition retains its existing GPLv3 license. Synced LGPL source and dependency
notices are documented in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

import os
import json
import logging
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QComboBox,
    QDoubleSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QPushButton,
    QApplication,
    QGroupBox,
    QMessageBox,
    QFileDialog,
    QCheckBox,
    QButtonGroup,
    QAbstractItemView,
)
from PyQt6.QtCore import Qt, QTimer, QObject, QEvent
import pyvista as pv
import numpy as np
from .utils import get_default_export_path

try:
    import nmrsim
    from nmrsim import Multiplet, Spectrum
except ImportError as e:
    logging.warning("NMR: nmrsim not available — multiplet simulation disabled (%s)", e)
    nmrsim = None
    Multiplet = None
    Spectrum = None

# Import RDKit for VDW radii calculation
try:
    from rdkit import Chem

    _pt = Chem.GetPeriodicTable()
    # Base VDW radii (scaled by 0.3 as in moledit core)
    VDW_RADII = {_pt.GetElementSymbol(i): _pt.GetRvdw(i) * 0.3 for i in range(1, 119)}
except ImportError:
    VDW_RADII = {"H": 1.2 * 0.3, "C": 1.7 * 0.3, "N": 1.55 * 0.3, "O": 1.52 * 0.3}

from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT,
)
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from .nmr_custom_ref_dialog import CustomReferenceDialog
from . import PLUGIN_VERSION



class NMRDialog(QDialog):
    """Enhanced NMR Chemical Shielding Dialog with Spectrum"""

    # Class-level constants for Nucleus Mapping and Physics
    ISOTOPE_MAP = {
        "H": "1H",
        "D": "2H",
        "T": "3H",
        "Li": "7Li",
        "Be": "9Be",
        "B": "11B",
        "C": "13C",
        "N": "15N",
        "O": "17O",
        "F": "19F",
        "Na": "23Na",
        "Mg": "25Mg",
        "Al": "27Al",
        "Si": "29Si",
        "P": "31P",
        "S": "33S",
        "Cl": "35Cl",
        "K": "39K",
        "Ca": "43Ca",
        "Sc": "45Sc",
        "Ti": "47Ti",
        "V": "51V",
        "Cr": "53Cr",
        "Mn": "55Mn",
        "Fe": "57Fe",
        "Co": "59Co",
        "Ni": "61Ni",
        "Cu": "63Cu",
        "Zn": "67Zn",
        "Ga": "71Ga",
        "Ge": "73Ge",
        "As": "75As",
        "Se": "77Se",
        "Br": "81Br",
        "Kr": "83Kr",
        "Rb": "87Rb",
        "Sr": "87Sr",
        "Y": "89Y",
        "Zr": "91Zr",
        "Nb": "93Nb",
        "Mo": "95Mo",
        "Tc": "99Tc",
        "Ru": "99Ru",
        "Rh": "103Rh",
        "Pd": "105Pd",
        "Ag": "109Ag",
        "Cd": "113Cd",
        "In": "115In",
        "Sn": "119Sn",
        "Sb": "121Sb",
        "Te": "125Te",
        "I": "127I",
        "Xe": "129Xe",
        "Cs": "133Cs",
        "Ba": "137Ba",
        "La": "139La",
        "W": "183W",
        "Os": "187Os",
        "Pt": "195Pt",
        "Au": "197Au",
        "Hg": "199Hg",
        "Tl": "205Tl",
        "Pb": "207Pb",
    }

    # Gyromagnetic ratios (10^6 rad s^-1 T^-1) approx -> Updated to IUPAC/CODATA standards
    GAMMA = {
        "H": 267.522,
        "1H": 267.522,
        "2H": 41.065,
        "3H": 285.35,
        "Li": 103.96,
        "7Li": 103.96,
        "9Be": 37.59,
        "B": 85.847,
        "11B": 85.847,
        "10B": 28.75,
        "C": 67.262,
        "13C": 67.262,
        "N": -27.116,
        "15N": -27.116,
        "14N": 19.331,
        "O": -36.26,
        "17O": -36.26,
        "F": 251.662,
        "19F": 251.662,
        "Na": 70.76,
        "23Na": 70.76,
        "25Mg": -16.38,
        "27Al": 69.76,
        "Si": -53.267,
        "29Si": -53.267,
        "P": 108.394,
        "31P": 108.394,
        "S": 20.53,
        "33S": 20.53,
        "Cl": 26.21,
        "35Cl": 26.21,
        "K": 12.48,
        "39K": 12.48,
        "Ca": -18.00,
        "43Ca": -18.00,
        "45Sc": 64.99,
        "47Ti": -15.08,
        "51V": 70.33,
        "Cr": -15.12,
        "53Cr": -15.12,
        "55Mn": 66.08,
        "57Fe": 8.66,
        "59Co": 63.17,
        "61Ni": -23.91,
        "63Cu": 70.90,
        "Zn": 16.74,
        "67Zn": 16.74,
        "71Ga": 81.74,
        "73Ge": -9.33,
        "75As": 45.80,
        "77Se": 51.203,
        "81Br": 72.24,
        "Rb": 87.53,
        "87Rb": 87.53,
        "87Sr": -11.57,
        "89Y": -13.11,
        "91Zr": -24.87,
        "93Nb": 65.47,
        "Mo": -17.46,
        "95Mo": -17.46,
        "99Ru": -12.29,
        "103Rh": -8.468,
        "105Pd": -12.24,
        "109Ag": -12.45,
        "Cd": -59.53,
        "113Cd": -59.53,
        "115In": 58.85,
        "119Sn": -99.95,
        "121Sb": 64.18,
        "125Te": -84.71,
        "127I": 53.71,
        "Xe": -73.99,
        "129Xe": -73.99,
        "133Cs": 35.09,
        "137Ba": 29.86,
        "139La": 37.90,
        "183W": 11.13,
        "187Os": 6.08,
        "Pt": 58.385,
        "195Pt": 58.385,
        "197Au": 4.71,
        "199Hg": 47.91,
        "205Tl": 155.19,
        "207Pb": 55.70,
    }

    def __init__(self, parent, data, couplings=None, file_path=None):
        super().__init__(parent)
        self.setWindowTitle(f"Calculated NMR Spectrum (v{PLUGIN_VERSION})")
        self.resize(600, 850)  # More compact width

        # Make dialog modeless (non-blocking)
        self.setWindowModality(Qt.WindowModality.NonModal)

        # Store parent dialog for 3D viewer access
        self.parent_dlg = parent

        self.data = data
        self.couplings = couplings if couplings else []
        self.displayed_data = list(data)

        # Track atom labels in 3D viewer
        self._atom_labels = []

        # Track selected peaks for highlighting
        self.selected_peak_indices = set()
        self.highlight_artists = []
        self.show_all_mode = False  # Track if showing all labels without highlights

        self.last_ref_name = None

        # Reference standards database (delta = reference position, sigma = isotropic shielding)
        # Chemical shift formula: δ_sample = δ_ref + (σ_ref - σ_sample)
        self.reference_standards = {
            "1H": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "TMS": {"delta_ref": 0.0, "sigma_ref": 31.8},
                "CDCl3": {"delta_ref": 7.26, "sigma_ref": 24.5},
                "DMSO-d6": {"delta_ref": 2.50, "sigma_ref": 29.3},
                "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0},
            },
            "13C": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "TMS": {"delta_ref": 0.0, "sigma_ref": 182.4},
                "CDCl3": {"delta_ref": 77.16, "sigma_ref": 105.2},
                "DMSO-d6": {"delta_ref": 39.52, "sigma_ref": 142.9},
                "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0},
            },
            "15N": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "CH3NO2": {"delta_ref": 0.0, "sigma_ref": -135.8},
                "NH3": {"delta_ref": -381.9, "sigma_ref": 244.4},
                "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0},
            },
            "31P": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "H3PO4 (85%)": {"delta_ref": 0.0, "sigma_ref": 328.4},
                "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0},
            },
            "19F": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "CFCl3": {"delta_ref": 0.0, "sigma_ref": 188.5},
                "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0},
            },
        }

        # Current reference values
        self.delta_ref = 0.0
        self.sigma_ref = 0.0

        # Spectrum settings
        self.linewidth = 1.0  # ppm for spectrum
        self.peak_intensity = 1.0

        # Settings file
        self.file_path = file_path
        if file_path:
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            file_dir = os.path.dirname(file_path)
            self.merged_peaks_file = os.path.join(
                file_dir, f"{base_name}-nmr_peak_info.json"
            )
        else:
            # Fallback if no file path provided
            self.merged_peaks_file = os.path.join(
                os.path.dirname(__file__), "nmr_merged_peaks.json"
            )

        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")

        # Track manually merged peaks: [{"indices": [0, 1, 2], "avg_delta": 7.5, ...}]
        self.merged_peaks = []
        self.load_merged_peaks()

        self.load_settings()

        # Initialize current nucleus before UI setup
        self.current_nucleus = "All"

        self.setup_ui()

        # Timer for polling main window selection (sync 3D -> NMR)
        self.sel_timer = QTimer(self)
        self.sel_timer.timeout.connect(self._check_external_selection)
        self.sel_timer.start(200)  # Check every 200ms

        # Custom 3D highlight actors and names
        self._nmr_sphere_actors = []
        self._nmr_label_names = []  # Explicitly track label names for removal



    def _check_external_selection(self):
        """Poll main window for 3D selection changes"""
        if not hasattr(self.parent_dlg, "mw"):
            return

        mw = self.parent_dlg.mw

        # 0. Check if Selection Mode is active (Optimization: Don't hijack selection if user is doing something else)
        # Assuming mw has a 'current_mode' or we check if 'SelectionTool' is active if that structure exists.
        # Fallback: If mw has 'selection_enabled' flag.
        # 0. Check if Selection Mode is active
        # We removed the mw.scene.mode check because it refers to the 2D editor mode.
        # 3D selection should be allowed unless we are in a specific conflicting mode like Measurement.

        if (
            getattr(mw.edit_3d_manager, "measurement_mode", False)
            if hasattr(mw, "edit_3d_manager")
            else False
        ):
            # Measurement mode uses a different selection list
            pass

        indices = set()

        e3d = getattr(mw, "edit_3d_manager", None)

        # Check standard 3D selection
        if e3d and e3d.selected_atoms_3d:
            indices.update(e3d.selected_atoms_3d)

        # Check measurement selection
        if e3d and e3d.selected_atoms_for_measurement:
            for item in e3d.selected_atoms_for_measurement:
                if isinstance(item, int):
                    indices.add(item)

        # 1. State Tracking for Stability
        # Compare current 3D selection with what we LAST knew about/set.
        # If they are identical, NO CHANGE has happened, so we do nothing.
        # This prevents the "Echo" loop where we clear the selection (setting it to empty)
        # and then this poller sees "Empty" vs "My Internal Selection" and clears the graph.

        current_mw_selection_frozen = frozenset(indices)

        # Initialize tracker if missing
        if getattr(self, "_last_synced_mw_selection", None) is None:
            self._last_synced_mw_selection = frozenset()

        # If the 3D selection HAS NOT CHANGED from what we last saw/set, STOP.
        if current_mw_selection_frozen == self._last_synced_mw_selection:
            return

        # Update our tracker to the new state
        self._last_synced_mw_selection = current_mw_selection_frozen

        # Calculate what the NMR selection SHOULD be based on the current 3D selection
        # Note: We use the "Any Member" rule for selecting, but by not expanding the 3D set,
        # unselection of the specific clicked atom correctly clears the indices.
        new_peak_selection = self._calculate_peak_selection_from_atoms(indices)

        # Only update if the peak selection itself has changed
        if new_peak_selection != self.selected_peak_indices:
            # If 3D selection is empty, force clear
            if not indices:
                self.clear_peak_selection()
            else:
                self.selected_peak_indices = new_peak_selection
                self.highlight_selected_peaks()
                # Update visual labels and spheres (Yellow),
                # but tell it NOT to sync back to the main window's selection set (Green).
                self.update_selected_labels(is_external_sync=True)

    def _calculate_peak_selection_from_atoms(self, target_atoms):
        """Helper to determine which peaks should be selected based on atom set"""
        if not getattr(self, "peaks_metadata", None):
            return set()

        new_selection = set()
        target_atoms = {int(i) for i in target_atoms}

        for peak_idx, metadata in enumerate(self.peaks_metadata):
            # metadata: (shift, intensity, is_merged, atom_indices)
            _, _, _, peak_atoms = metadata
            peak_atoms_set = {int(i) for i in peak_atoms}

            # Selection Rule: A peak is selected if ANY of its atoms are in the 3D selection set.
            # This is stable and prevents the "flicker" caused by switching between ANY and SUBSET.
            if not peak_atoms_set.isdisjoint(target_atoms):
                new_selection.add(peak_idx)
        return new_selection

    def select_peaks_by_atom_indices(self, atom_indices):
        """Deprecated/Legacy: Now uses _check_external_selection logic directly"""
        # Kept for potential internal calls, but redirected to robust logic
        new_peaks = self._calculate_peak_selection_from_atoms(atom_indices)
        if new_peaks != self.selected_peak_indices:
            self.selected_peak_indices = new_peaks
            self.highlight_selected_peaks()
            self.update_selected_labels()

    def get_nucleus_key(self, atom_sym):
        """Map atom symbol to nucleus key for reference standards"""
        # Clean the input symbol (remove whitespace, numbers if accidentally passed)
        clean_sym = "".join([c for c in atom_sym if c.isalpha()])
        return self.ISOTOPE_MAP.get(
            clean_sym.upper(), self.ISOTOPE_MAP.get(clean_sym, clean_sym)
        )

    def load_settings(self):
        """Load NMR settings from JSON"""
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r") as f:
                    settings = json.load(f)

                nmr_settings = settings.get("nmr_settings", {})
                self.linewidth = nmr_settings.get("spectrum_linewidth", 1.0)
                self.peak_intensity = nmr_settings.get("peak_intensity", 1.0)

                self.last_ref_name = nmr_settings.get("last_reference", None)

                # Load custom references
                custom_refs = nmr_settings.get("custom_references", {})
                for nucleus, refs in custom_refs.items():
                    if nucleus not in self.reference_standards:
                        self.reference_standards[nucleus] = {}
                    for ref_name, ref_val in refs.items():
                        self.reference_standards[nucleus][ref_name] = ref_val
            except Exception as e:
                print(f"Error loading NMR settings: {e}")

    def merge_selected_peaks(self):
        """Merge selected peaks into a single entry with isotope validation"""
        if len(self.selected_peak_indices) < 2:
            QMessageBox.warning(
                self, "Invalid Selection", "Please select at least 2 peaks to merge."
            )
            return

        selected_indices = []
        atom_symbols = set()

        # 1. 選択されたピークから原子インデックスと元素記号を取得
        if getattr(self, "peaks_metadata", None) is not None and self.peaks_metadata:
            for peak_idx in self.selected_peak_indices:
                if peak_idx < len(self.peaks_metadata):
                    _, _, _, atom_indices = self.peaks_metadata[peak_idx]
                    selected_indices.extend(atom_indices)

                    # 核種のチェック用に元素記号を収集
                    for idx in atom_indices:
                        item = next(
                            (d for d in self.data if d.get("atom_idx", None) == idx),
                            None,
                        )
                        if item:
                            atom_symbols.add(item.get("atom_sym", None))

        # 2. 【物理バリデーション】異なる元素が混ざっていないかチェック
        if len(atom_symbols) > 1:
            QMessageBox.critical(
                self,
                "Physical Inconsistency",
                f"Cannot merge different nuclei: {', '.join(atom_symbols)}. "
                "NMR peaks can only be merged for the same isotope.",
            )
            return

        # 3. 重複排除とソート
        selected_indices = sorted(list(set(selected_indices)))

        # 4. 既存のマージグループとの競合チェック
        new_merged_peaks = []
        conflict_found = False
        for group in self.merged_peaks:
            if any(idx in group["indices"] for idx in selected_indices):
                conflict_found = True
                continue  # 競合する古いグループはスキップ（後でまとめて追加するため）
            new_merged_peaks.append(group)

        if conflict_found:
            reply = QMessageBox.question(
                self,
                "Merge Conflict",
                "Some selected atoms are already part of another group. Replace existing merge?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if reply == QMessageBox.StandardButton.No:
                return

        # 5. 新しいグループの追加と保存
        new_merged_peaks.append({"indices": selected_indices})
        self.merged_peaks = new_merged_peaks
        self.save_merged_peaks()

        # 6. UIのクリーンアップ
        self.clear_peak_selection()
        if hasattr(self.parent_dlg, "mw") and hasattr(self.parent_dlg.mw, "statusBar"):
            self.parent_dlg.mw.statusBar().showMessage(
                f"Merged {len(selected_indices)} atoms into one peak.", 5000
            )

        self.recalc()

    def unmerge_selected_peaks(self):
        """Separate previously merged peaks back into individuals"""
        if not self.selected_peak_indices:
            return

        groups_to_remove = []
        for peak_idx in self.selected_peak_indices:
            if peak_idx < len(self.peaks_metadata):
                _, _, is_merged, atom_indices = self.peaks_metadata[peak_idx]
                if is_merged:
                    # Find which group in self.merged_peaks contains these indices
                    # Since groups are unique per atom, we can just match any index
                    for group in self.merged_peaks:
                        if all(idx in group["indices"] for idx in atom_indices) and len(
                            group["indices"]
                        ) == len(atom_indices):
                            if group not in groups_to_remove:
                                groups_to_remove.append(group)
                            break

        if not groups_to_remove:
            return

        for group in groups_to_remove:
            self.merged_peaks.remove(group)

        self.save_merged_peaks()
        self.clear_peak_selection()
        self.recalc()

    def save_merged_peaks(self):
        """Save merged peaks to JSON file"""
        try:
            with open(self.merged_peaks_file, "w") as f:
                json.dump(self.merged_peaks, f, indent=2)
        except Exception as e:
            print(f"Error saving merged peaks: {e}")

    def load_merged_peaks(self):
        """Load merged peaks from JSON file"""
        if os.path.exists(self.merged_peaks_file):
            try:
                with open(self.merged_peaks_file, "r") as f:
                    self.merged_peaks = json.load(f)
            except Exception as e:
                print(f"Error loading merged peaks: {e}")
                self.merged_peaks = []
        else:
            self.merged_peaks = []

    def save_settings(self):
        """Save NMR settings to JSON"""
        all_settings = {}
        # 既存の設定を読み込む（MO設定などを消さないため）
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r") as f:
                    all_settings = json.load(f)
            except Exception as _e:
                logging.warning("[nmr_analysis.py:538] silenced: %s", _e)

        # Extract custom references only (non-default)
        # Extract custom references only (non-default)
        default_standards = {
            "1H": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "TMS": {"delta_ref": 0.0, "sigma_ref": 31.8},
                "CDCl3": {"delta_ref": 7.26, "sigma_ref": 24.5},
                "DMSO-d6": {"delta_ref": 2.50, "sigma_ref": 29.3},
            },
            "13C": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "TMS": {"delta_ref": 0.0, "sigma_ref": 182.4},
                "CDCl3": {"delta_ref": 77.16, "sigma_ref": 105.2},
                "DMSO-d6": {"delta_ref": 39.52, "sigma_ref": 142.9},
            },
            "15N": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "CH3NO2": {"delta_ref": 0.0, "sigma_ref": -135.8},
                "NH3": {"delta_ref": -381.9, "sigma_ref": 244.4},
            },
            "31P": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "H3PO4 (85%)": {"delta_ref": 0.0, "sigma_ref": 328.4},
            },
            "19F": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "CFCl3": {"delta_ref": 0.0, "sigma_ref": 188.5},
            },
        }

        custom_refs = {}
        for nucleus, refs in self.reference_standards.items():
            for ref_name, ref_val in refs.items():
                if ref_name == "Custom":
                    continue  # Skip custom placeholder
                default_dict = default_standards.get(nucleus, {})
                if ref_name not in default_dict:
                    # Completely custom reference
                    if nucleus not in custom_refs:
                        custom_refs[nucleus] = {}
                    custom_refs[nucleus][ref_name] = ref_val
                elif default_dict.get(ref_name, None) != ref_val:
                    # Modified default reference
                    if nucleus not in custom_refs:
                        custom_refs[nucleus] = {}
                    custom_refs[nucleus][ref_name] = ref_val

        current_nmr_settings = {
            "spectrum_linewidth": self.linewidth,
            "peak_intensity": self.peak_intensity,
            "last_reference": self.last_ref_name,
            "custom_references": custom_refs,  # 既存の変数
        }

        # 全体設定の 'nmr_settings' キーだけを更新
        all_settings["nmr_settings"] = current_nmr_settings

        try:
            with open(self.settings_file, "w") as f:
                json.dump(all_settings, f, indent=2)
        except Exception as e:
            print(f"Error saving NMR settings: {e}")

    def reset_zoom(self, event):
        """Reset plot zoom on double click"""
        if event.dblclick:
            self.toolbar.home()

    def setup_ui(self):
        main_layout = QVBoxLayout(self)

        # 1. Reference & Element Filter Row
        top_row = QHBoxLayout()

        # Nucleus filter with toggle buttons
        nucleus_box = QGroupBox("Nucleus Filter")
        nucleus_layout = QHBoxLayout(nucleus_box)

        # Get available nuclei
        nuclei = ["All"] + sorted(list(set([d["atom_sym"] for d in self.data])))

        # Create button group for exclusive selection
        self.nucleus_button_group = QButtonGroup()
        self.nucleus_buttons = {}

        for nucleus in nuclei:
            btn = QPushButton(nucleus)
            btn.setCheckable(True)
            btn.setAutoDefault(False)
            btn.setStyleSheet("""
                QPushButton {
                    background-color: #f8f8f8;
                    border: 1px solid #ccc;
                    border-radius: 4px;
                    padding: 6px 15px;
                    font-size: 10pt;
                    font-weight: 500;
                }
                QPushButton:checked {
                    background-color: #0066cc;
                    color: white;
                    border-color: #004d99;
                    font-weight: bold;
                }
                QPushButton:hover {
                    border-color: #0066cc;
                    background-color: #eef6ff;
                }
                QPushButton:checked:hover {
                    background-color: #0059b3;
                }
            """)
            btn.toggled.connect(
                lambda checked, n=nucleus: self.on_nucleus_changed(n)
                if checked
                else None
            )
            self.nucleus_button_group.addButton(btn)
            self.nucleus_buttons[nucleus] = btn
            nucleus_layout.addWidget(btn)

        nucleus_layout.addStretch()
        top_row.addWidget(nucleus_box)

        # Reference selection and values (merged into single group)
        ref_box = QGroupBox("Reference Standard")
        ref_layout = QVBoxLayout(ref_box)

        # Reference selection row
        ref_sel_row = QHBoxLayout()
        ref_sel_row.addWidget(QLabel("Standard:"))
        self.combo_ref = QComboBox()
        self.combo_ref.setMinimumWidth(250)
        from PyQt6.QtWidgets import QSizePolicy

        self.combo_ref.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.combo_ref.currentIndexChanged.connect(self.on_ref_change)
        ref_sel_row.addWidget(self.combo_ref)

        btn_add_ref = QPushButton("+ Custom")
        btn_add_ref.setFixedWidth(80)
        btn_add_ref.setAutoDefault(False)
        btn_add_ref.clicked.connect(self.add_custom_reference)
        ref_sel_row.addWidget(btn_add_ref)

        btn_del_ref = QPushButton("Delete")
        btn_del_ref.setFixedWidth(60)
        btn_del_ref.setAutoDefault(False)
        btn_del_ref.setToolTip("Delete selected custom reference")
        btn_del_ref.clicked.connect(self.delete_custom_reference)
        ref_sel_row.addWidget(btn_del_ref)

        # Duplicate combo_ref removed

        ref_layout.addLayout(ref_sel_row)

        # Delta ref (reference peak position)
        delta_row = QHBoxLayout()
        delta_row.addWidget(QLabel("δ_ref (ppm):"))
        self.spin_delta_ref = QDoubleSpinBox()
        self.spin_delta_ref.setRange(-500, 500)
        self.spin_delta_ref.setValue(self.delta_ref)
        self.spin_delta_ref.setDecimals(2)
        self.spin_delta_ref.setToolTip(
            "Reference peak position on chemical shift scale"
        )
        self.spin_delta_ref.valueChanged.connect(self.on_ref_value_change)
        delta_row.addWidget(self.spin_delta_ref)
        ref_layout.addLayout(delta_row)

        # Sigma ref (isotropic shielding)
        sigma_row = QHBoxLayout()
        sigma_row.addWidget(QLabel("σ_ref (ppm):"))
        self.spin_sigma_ref = QDoubleSpinBox()
        self.spin_sigma_ref.setRange(-10000, 20000)
        self.spin_sigma_ref.setValue(self.sigma_ref)
        self.spin_sigma_ref.setDecimals(2)
        self.spin_sigma_ref.setToolTip(
            "Isotropic shielding value of reference compound"
        )
        self.spin_sigma_ref.valueChanged.connect(self.on_ref_value_change)
        sigma_row.addWidget(self.spin_sigma_ref)
        ref_layout.addLayout(sigma_row)

        top_row.addWidget(ref_box)

        main_layout.addLayout(top_row)

        # 2. Spectrum Plot
        spec_group = QGroupBox("NMR Stick Spectrum (δ = δ_ref + σ_ref - σ)")
        spec_layout = QVBoxLayout(spec_group)

        # Spectrum settings row
        spec_settings = QHBoxLayout()

        # Add checkbox to show all labels
        self.chk_show_all_labels = QCheckBox("Show all atom labels")
        self.chk_show_all_labels.stateChanged.connect(self.toggle_all_labels)
        spec_settings.addWidget(self.chk_show_all_labels)

        # Add button to clear selection
        btn_clear_selection = QPushButton("Clear Selection")
        btn_clear_selection.setFixedWidth(120)
        btn_clear_selection.setAutoDefault(False)
        btn_clear_selection.setToolTip("Clear all selected peaks and labels")
        btn_clear_selection.clicked.connect(self.clear_peak_selection)
        spec_settings.addWidget(btn_clear_selection)

        # Add merge selected button
        btn_merge = QPushButton("Merge Selected")
        btn_merge.setFixedWidth(120)
        btn_merge.setAutoDefault(False)
        btn_merge.setToolTip("Merge selected peaks into one entry")
        btn_merge.clicked.connect(self.merge_selected_peaks)
        spec_settings.addWidget(btn_merge)

        btn_unmerge = QPushButton("Unmerge")
        btn_unmerge.setFixedWidth(120)
        btn_unmerge.setAutoDefault(False)
        btn_unmerge.setToolTip("Separate merged peaks back to individuals")
        btn_unmerge.clicked.connect(self.unmerge_selected_peaks)
        spec_settings.addWidget(btn_unmerge)

        # new line for Real Spectrum controls
        real_spec_layout = QHBoxLayout()

        self.chk_real_spectrum = QCheckBox("Simulate Coupling")
        sim_tooltip = "Simulate multiplets using J-couplings"
        if not nmrsim:
            sim_tooltip += " (Requires nmrsim library)"
            self.chk_real_spectrum.setEnabled(False)
        self.chk_real_spectrum.setToolTip(sim_tooltip)
        self.chk_real_spectrum.stateChanged.connect(self.toggle_simulation_controls)
        self.chk_real_spectrum.stateChanged.connect(self.plot_spectrum)
        real_spec_layout.addWidget(self.chk_real_spectrum)

        # Broadening control
        self.lbl_width = QLabel(" Width (Hz):")
        real_spec_layout.addWidget(self.lbl_width)

        self.spin_real_width = QDoubleSpinBox()
        self.spin_real_width.setRange(0.01, 100.0)
        self.spin_real_width.setValue(1.0)
        self.spin_real_width.setSingleStep(0.5)
        self.spin_real_width.setFixedWidth(100)
        if not nmrsim:
            self.spin_real_width.setEnabled(False)
        self.spin_real_width.valueChanged.connect(lambda val: self.plot_spectrum())
        real_spec_layout.addWidget(self.spin_real_width)

        real_spec_layout.addWidget(QLabel(" Spectrometer (MHz):"))
        self.spin_mhz = QDoubleSpinBox()
        self.spin_mhz.setRange(10.0, 2000.0)
        self.spin_mhz.setValue(400.0)
        self.spin_mhz.setSingleStep(10.0)
        self.spin_mhz.setFixedWidth(100)
        if not nmrsim:
            self.spin_mhz.setEnabled(False)
        self.spin_mhz.valueChanged.connect(lambda val: self.plot_spectrum())
        real_spec_layout.addWidget(self.spin_mhz)

        # Initial state check
        self.toggle_simulation_controls()

        real_spec_layout.addStretch()
        spec_layout.addLayout(real_spec_layout)

        # X-Axis Range row
        x_range_layout = QHBoxLayout()
        self.chk_auto_x = QCheckBox("Auto Range")
        self.chk_auto_x.setChecked(False)
        self.chk_auto_x.stateChanged.connect(lambda: self.plot_spectrum())
        x_range_layout.addWidget(self.chk_auto_x)

        x_range_layout.addWidget(QLabel(" X Range:"))
        self.spin_x_max = QDoubleSpinBox()
        self.spin_x_max.setRange(-2000, 2000)
        self.spin_x_max.setValue(12.0)
        self.spin_x_max.valueChanged.connect(lambda: self.plot_spectrum())
        x_range_layout.addWidget(self.spin_x_max)

        x_range_layout.addWidget(QLabel("to"))
        self.spin_x_min = QDoubleSpinBox()
        self.spin_x_min.setRange(-2000, 2000)
        self.spin_x_min.setValue(-1.0)
        self.spin_x_min.valueChanged.connect(lambda: self.plot_spectrum())
        x_range_layout.addWidget(self.spin_x_min)
        x_range_layout.addWidget(QLabel("ppm"))

        def update_x_spins_enabled():
            enabled = not self.chk_auto_x.isChecked()
            self.spin_x_max.setEnabled(enabled)
            self.spin_x_min.setEnabled(enabled)

        self.chk_auto_x.stateChanged.connect(update_x_spins_enabled)
        update_x_spins_enabled()

        x_range_layout.addStretch()
        spec_layout.addLayout(x_range_layout)

        spec_settings.addStretch()

        # Add spec_settings to layout if not already (it was implied in previous context to be added)
        # spec_layout.addLayout(spec_settings)

        btn_export = QPushButton("Export Image")
        btn_export.setAutoDefault(False)
        btn_export.clicked.connect(self.export_spectrum)
        spec_settings.addWidget(btn_export)

        btn_export_csv = QPushButton("Export CSV")
        btn_export_csv.setAutoDefault(False)
        btn_export_csv.clicked.connect(self.export_spectrum_csv)
        spec_settings.addWidget(btn_export_csv)

        spec_layout.addLayout(spec_settings)

        # Matplotlib canvas - adjusted for narrower dialog
        self.figure = Figure(figsize=(5.5, 4))  # Narrower to fit 600px width
        self.canvas = FigureCanvas(self.figure)
        self.canvas.mpl_connect("button_press_event", self.reset_zoom)
        self.canvas.mpl_connect("button_press_event", self.on_peak_click)

        # Add navigation toolbar for zoom/pan
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        spec_layout.addWidget(self.toolbar)
        spec_layout.addWidget(self.canvas)

        main_layout.addWidget(spec_group)

        # 3. Data Table
        table_group = QGroupBox("Chemical Shift Data")
        table_layout = QVBoxLayout(table_group)

        self.table = QTableWidget()
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(
            ["Idx", "Nucleus", "σ (ppm)", "δ (ppm)", "J (Hz)"]
        )
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self.table.verticalHeader().setVisible(False)
        self.table.setMinimumHeight(100)  # Half size for compact layout
        table_layout.addWidget(self.table)

        table_btn_row = QHBoxLayout()
        btn_copy = QPushButton("Copy to Clipboard")
        btn_copy.clicked.connect(self.copy_table)
        table_btn_row.addWidget(btn_copy)

        btn_export_table_csv = QPushButton("Export CSV")
        btn_export_table_csv.clicked.connect(self.export_table_csv)
        table_btn_row.addWidget(btn_export_table_csv)
        table_btn_row.addStretch()
        table_layout.addLayout(table_btn_row)

        main_layout.addWidget(table_group)

        # 4. Bottom Buttons
        bottom_row = QHBoxLayout()
        bottom_row.addStretch()
        btn_close = QPushButton("Close")
        btn_close.setFixedWidth(100)
        btn_close.clicked.connect(self.close)
        bottom_row.addWidget(btn_close)
        main_layout.addLayout(bottom_row)

        # Finally select "All" by default to trigger population
        if "All" in self.nucleus_buttons:
            self.nucleus_buttons["All"].setChecked(True)

    def update_reference_combo(self):
        """Update reference combo box based on current nucleus"""
        if getattr(self, "combo_ref", None) is None:
            return

        current_nucleus = self.current_nucleus
        if current_nucleus == "All":
            # Special case for "All": Only allow "No Reference"
            # We want to visualize raw values (no shifting) across different nuclei
            self.combo_ref.blockSignals(True)
            self.combo_ref.clear()
            self.combo_ref.addItems(["No Reference"])
            self.combo_ref.setCurrentText("No Reference")
            self.combo_ref.blockSignals(False)

            # Set values to 0,0
            self.delta_ref = 0.0
            self.sigma_ref = 0.0

            # Update spinboxes
            self.spin_delta_ref.blockSignals(True)
            self.spin_sigma_ref.blockSignals(True)
            self.spin_delta_ref.setValue(0.0)
            self.spin_sigma_ref.setValue(0.0)
            self.spin_delta_ref.blockSignals(False)
            self.spin_sigma_ref.blockSignals(False)

            # Disable inputs
            self.spin_delta_ref.setEnabled(False)
            self.spin_sigma_ref.setEnabled(False)
            return

        # Map atom symbol to nucleus key (e.g., "H" -> "1H")
        current_nucleus = self.get_nucleus_key(current_nucleus)

        # Save current selection to preserve it if possible
        current_ref = (
            self.combo_ref.currentText() if self.combo_ref.count() > 0 else None
        )

        # Block signals to prevent triggering on_ref_change during population
        self.combo_ref.blockSignals(True)

        # Clear and repopulate
        self.combo_ref.clear()
        refs = self.reference_standards.get(current_nucleus, {})

        # Create list of items
        items = list(refs.keys())

        if "No Reference" not in items:
            items.insert(0, "No Reference")

        # Always ensure "Custom" is in the list
        if "Custom" not in items:
            items.append("Custom")

        self.combo_ref.addItems(items)

        target_ref = None

        # 直前の選択がある、かつ "No Reference" ではない（Allモード由来ではない）場合
        if current_ref and current_ref in items and current_ref != "No Reference":
            target_ref = current_ref
        # 記憶していたリファレンスがある場合
        elif self.last_ref_name and self.last_ref_name in items:
            target_ref = self.last_ref_name
        # デフォルト(TMS)
        elif "TMS" in items:
            target_ref = "TMS"
        # それ以外
        elif items:
            target_ref = items[0]

        if target_ref:
            self.combo_ref.setCurrentText(target_ref)
            ref_data = refs.get(target_ref, {"delta_ref": 0.0, "sigma_ref": 0.0})
        else:
            self.combo_ref.setCurrentText("Custom")  # または No Reference
            ref_data = {"delta_ref": 0.0, "sigma_ref": 0.0}

        # Update internal values
        self.delta_ref = ref_data.get("delta_ref", 0.0)
        self.sigma_ref = ref_data.get("sigma_ref", 0.0)

        # Update spinboxes (these updates won't trigger recalc since we blocked combo signals)
        self.spin_delta_ref.blockSignals(True)
        self.spin_sigma_ref.blockSignals(True)
        self.spin_delta_ref.setValue(self.delta_ref)
        self.spin_sigma_ref.setValue(self.sigma_ref)
        # Only allow editing if "Custom" is selected
        is_custom = self.combo_ref.currentText() == "Custom"
        self.spin_delta_ref.setEnabled(is_custom)
        self.spin_sigma_ref.setEnabled(is_custom)
        self.spin_delta_ref.blockSignals(False)
        self.spin_sigma_ref.blockSignals(False)

        # Auto Range if "No Reference" is active
        if self.combo_ref.currentText() == "No Reference":
            if getattr(self, "chk_auto_x", None) is not None:
                self.chk_auto_x.setChecked(True)

        # Re-enable combo signals
        self.combo_ref.blockSignals(False)

    def on_ref_change(self):
        """Handle reference standard change"""
        current_nucleus = self.current_nucleus
        if current_nucleus != "All":
            self.last_ref_name = self.combo_ref.currentText()
        if current_nucleus == "All":
            # Force 0,0 for All view
            self.delta_ref = 0.0
            self.sigma_ref = 0.0
            self.recalc()
            return

        # Map atom symbol to nucleus key
        nucleus_key = self.get_nucleus_key(current_nucleus)

        ref_name = self.combo_ref.currentText()
        self.last_ref_name = ref_name

        # Auto Range if "No Reference" is selected
        if ref_name == "No Reference":
            if getattr(self, "chk_auto_x", None) is not None:
                self.chk_auto_x.setChecked(True)

        refs = self.reference_standards.get(nucleus_key, {})
        ref_data = refs.get(ref_name, {"delta_ref": 0.0, "sigma_ref": 0.0})

        self.delta_ref = ref_data["delta_ref"]
        self.sigma_ref = ref_data["sigma_ref"]

        # Block signals only during setValue to prevent triggering on_ref_value_change
        self.spin_delta_ref.blockSignals(True)
        self.spin_sigma_ref.blockSignals(True)
        self.spin_delta_ref.setValue(self.delta_ref)
        self.spin_sigma_ref.setValue(self.sigma_ref)
        self.spin_delta_ref.blockSignals(False)
        self.spin_sigma_ref.blockSignals(False)

        # Only allow editing if "Custom" is selected
        is_custom = ref_name == "Custom"
        self.spin_delta_ref.setEnabled(is_custom)
        self.spin_sigma_ref.setEnabled(is_custom)

        self.recalc()

    def on_ref_value_change(self):
        """Handle manual reference value changes"""
        self.delta_ref = self.spin_delta_ref.value()
        self.sigma_ref = self.spin_sigma_ref.value()

        # Update the reference standard dict
        current_nucleus = self.current_nucleus
        ref_name = self.combo_ref.currentText()
        if current_nucleus != "All":
            # Map atom symbol to nucleus key
            nucleus_key = self.get_nucleus_key(current_nucleus)
            if nucleus_key not in self.reference_standards:
                self.reference_standards[nucleus_key] = {}
            self.reference_standards[nucleus_key][ref_name] = {
                "delta_ref": self.delta_ref,
                "sigma_ref": self.sigma_ref,
            }
        self.recalc()

    def add_custom_reference(self):
        """Add a custom reference standard using custom dialog"""
        # Get available nuclei from data
        available_nuclei = sorted(list(set([d["atom_sym"] for d in self.data])))
        if not available_nuclei:
            available_nuclei = ["1H", "13C", "15N", "31P", "19F"]

        # Show custom dialog
        dialog = CustomReferenceDialog(self, available_nuclei)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return

        ref_name, nucleus_data = dialog.get_reference_data()

        # Add to standards for each nucleus
        for raw_nucleus, values in nucleus_data.items():
            # Ensure we use the standardized key (e.g. "H" -> "1H")
            nucleus = self.get_nucleus_key(raw_nucleus)

            if nucleus not in self.reference_standards:
                self.reference_standards[nucleus] = {}
            self.reference_standards[nucleus][ref_name] = values

        self.save_settings()

        # If we are in "All" mode, switch to the specific nucleus of the added reference
        # This prevents the reference from being hidden by the "All" restriction ("No Reference" only)
        if self.current_nucleus == "All" and nucleus_data:
            # Pick the first nucleus added (e.g. "1H" or "H")
            first_raw = list(nucleus_data.keys())[0]
            # Get the key used for buttons (e.g. "H" or "1H") - keys in nucleus_buttons usually match atom symbols
            # We should try to find the matching button.
            target_btn_key = None

            # Try standardized key first
            std_key = self.get_nucleus_key(first_raw)
            if std_key in self.nucleus_buttons:
                target_btn_key = std_key
            elif first_raw in self.nucleus_buttons:
                target_btn_key = first_raw

            if target_btn_key:
                self.nucleus_buttons[target_btn_key].setChecked(True)
                # Manually trigger the mode change handler since setChecked might not if via code (depends on signal)
                # But usually the group handles it. Let's explicitly call the update to be safe.
                self.current_nucleus = target_btn_key

        # Force UI update
        self.update_reference_combo()

        # Select the new reference
        # Now that we switched inputs (if we were in All), this should work
        self.combo_ref.setCurrentText(ref_name)

        if self.parent_dlg and hasattr(self.parent_dlg, "mw"):
            self.parent_dlg.mw.statusBar().showMessage(
                f"Added reference '{ref_name}' for {len(nucleus_data)} nucleus/nuclei.",
                5000,
            )
        else:
            print(f"Added reference '{ref_name}'")

    def delete_custom_reference(self):
        """Delete currently selected custom reference"""
        if not getattr(self, "current_nucleus", None):
            return

        current_nucleus = (
            self.get_nucleus_key(self.current_nucleus)
            if self.current_nucleus != "All"
            else "1H"
        )
        current_ref = self.combo_ref.currentText()

        if not current_ref:
            return

        # Check if it is a built-in standard
        default_standards = {
            "1H": ["TMS", "CDCl3", "DMSO-d6"],
            "13C": ["TMS", "CDCl3", "DMSO-d6"],
            "15N": ["CH3NO2", "NH3"],
            "31P": ["H3PO4 (85%)"],
            "19F": ["CFCl3"],
        }

        # Also check hardcoded defaults in save_settings to be safe
        is_default = False
        if (
            current_nucleus in default_standards
            and current_ref in default_standards[current_nucleus]
        ):
            is_default = True

        # Prevent deletion of "Custom" placeholder AND "No Reference"
        if current_ref in ["Custom", "No Reference"]:
            is_default = True

        if is_default:
            QMessageBox.warning(
                self,
                "Cannot Delete",
                f"'{current_ref}' is a built-in standard and cannot be deleted.",
            )
            return

        # Confirm deletion
        reply = QMessageBox.question(
            self,
            "Confirm Deletion",
            f"Are you sure you want to delete custom reference '{current_ref}'?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )

        if reply == QMessageBox.StandardButton.Yes:
            # Remove from storage
            if current_nucleus in self.reference_standards:
                if current_ref in self.reference_standards[current_nucleus]:
                    del self.reference_standards[current_nucleus][current_ref]
                    self.save_settings()
                    self.update_reference_combo()
                    if self.parent_dlg and hasattr(self.parent_dlg, "mw"):
                        self.parent_dlg.mw.statusBar().showMessage(
                            f"Reference '{current_ref}' removed.", 5000
                        )
                    else:
                        print(f"Reference '{current_ref}' removed.")
                else:
                    # Should not happen if UI is consistent
                    QMessageBox.warning(
                        self, "Error", "Reference not found in storage."
                    )

    def on_spectrum_settings_change(self):
        """Handle spectrum visualization parameter changes"""
        self.linewidth = self.spin_linewidth.value()
        self.peak_intensity = self.spin_intensity.value()
        self.save_settings()
        self.plot_spectrum()

    def on_nucleus_changed(self, nucleus):
        """Handle nucleus button toggle"""
        self.current_nucleus = nucleus
        self.clear_peak_selection()

        if getattr(self, "chk_auto_x", None) is not None:
            self.chk_auto_x.blockSignals(True)
            self.chk_auto_x.setChecked(False)  # ここで強制解除
            self.chk_auto_x.blockSignals(False)

            # Auto RangeをOFFにしたので、入力欄（スピンボックス）は有効化しておく
            if getattr(self, "spin_x_max", None) is not None:
                self.spin_x_max.setEnabled(True)
            if getattr(self, "spin_x_min", None) is not None:
                self.spin_x_min.setEnabled(True)

        # Update default X range for this nucleus
        self.update_x_range_defaults(nucleus)

        self.apply_filter()

        # User Request: Auto-enable coupling when switching to any nucleus (if coupling exists)
        # We do this AFTER apply_filter (which calls recalc) so we know if coupons exist for THIS nucleus.
        if (
            nucleus != "All"
            and getattr(self, "chk_real_spectrum", None) is not None
            and self.chk_real_spectrum.isEnabled()
        ):
            self.chk_real_spectrum.blockSignals(True)
            self.chk_real_spectrum.setChecked(True)
            self.chk_real_spectrum.blockSignals(False)
            # Since we blocked signals, manually trigger one plot update if we just checked it
            self.plot_spectrum()

    def update_x_range_defaults(self, nucleus):
        """Set appropriate default X range based on nucleus type"""
        if getattr(self, "spin_x_max", None) is None or not hasattr(self, "spin_x_min"):
            return

        # Block signals to prevent intermediate plotting (fixes "slope graph" artifact)
        self.spin_x_max.blockSignals(True)
        self.spin_x_min.blockSignals(True)

        # Ranges based on typical chemical shifts (ppm)
        defaults = {
            "H": (12, -1),
            "1H": (12, -1),
            "C": (220, -10),
            "13C": (220, -10),
            "N": (500, -400),
            "15N": (500, -400),
            "F": (0, -200),
            "19F": (0, -200),
            "P": (250, -150),
            "31P": (250, -150),
            "O": (1000, -500),
            "17O": (1000, -500),
            "Si": (100, -200),
            "29Si": (100, -200),
            "Pt": (1000, -5000),
        }

        # Get nucleus key to handle isotopes
        key = self.get_nucleus_key(nucleus)
        # Try key, then generic element (e.g. "Pt" if "195Pt" not in defaults), then fallback
        import re

        element_only = re.sub(r"[^A-Z]", "", key)

        range_vals = defaults.get(key, defaults.get(element_only, (500, -500)))

        self.spin_x_max.setValue(range_vals[0])
        self.spin_x_min.setValue(range_vals[1])

        self.spin_x_max.blockSignals(False)
        self.spin_x_min.blockSignals(False)

    def apply_filter(self):
        """Filter data by nucleus and update UI"""
        if getattr(self, "table", None) is None:
            return

        nucleus = self.current_nucleus
        if nucleus == "All":
            self.displayed_data = self.data
        else:
            self.displayed_data = [d for d in self.data if d["atom_sym"] == nucleus]

        # Update reference combo for this nucleus
        self.update_reference_combo()

        # Force recalculation with new reference values
        self.recalc()

    def recalc(self):
        """Recalculate chemical shifts and update table + spectrum"""
        # Get merged peak groups
        merged_indices = set()
        for group in self.merged_peaks:
            merged_indices.update(group["indices"])

        # Track which rows to display
        rows_to_display = []

        # Add merged peaks first
        for group in self.merged_peaks:
            # Calculate averages dynamically based on current reference
            total_sigma = 0.0
            total_delta = 0.0
            count = 0

            for atom_idx in group["indices"]:
                item = next(
                    (d for d in self.data if d.get("atom_idx", None) == atom_idx), None
                )
                if item:
                    sigma = item.get("shielding", 0.0)
                    delta = self.delta_ref + (self.sigma_ref - sigma)
                    total_sigma += sigma
                    total_delta += delta
                    count += 1

            if count > 0:
                avg_sigma = total_sigma / count
                avg_delta = total_delta / count

                # Check if any atoms in this group are in displayed_data
                group_items = [
                    item
                    for item in self.displayed_data
                    if item.get("atom_idx", -1) in group["indices"]
                ]
                if group_items:
                    rows_to_display.append(
                        (
                            "merged",
                            {
                                "indices": group["indices"],
                                "avg_sigma": avg_sigma,
                                "avg_delta": avg_delta,
                            },
                            group_items,
                        )
                    )

        # Add individual non-merged items
        for item in self.displayed_data:
            if item.get("atom_idx", -1) not in merged_indices:
                rows_to_display.append(("individual", item, None))

        # Update table
        self.table.setRowCount(len(rows_to_display))

        for r, row_data in enumerate(rows_to_display):
            row_type, data, group_items = row_data

            if row_type == "merged":
                # Display merged group
                indices_str = ", ".join([str(idx) for idx in data["indices"]])
                self.table.setItem(r, 0, QTableWidgetItem(f"[{indices_str}]"))
                self.table.setItem(
                    r,
                    1,
                    QTableWidgetItem(
                        f"{len(data['indices'])}{group_items[0].get('atom_sym', '')}"
                    ),
                )
                self.table.setItem(r, 2, QTableWidgetItem(f"{data['avg_sigma']:.2f}"))
                self.table.setItem(r, 3, QTableWidgetItem(f"{data['avg_delta']:.2f}"))
                # J-Coupling for merged peaks: Show range or 'Complex'
                # Attempt to show formatted list for the first representative or combined
                j_str = self.get_j_coupling_string(data["indices"])
                self.table.setItem(r, 4, QTableWidgetItem(j_str))
            else:
                # Display individual item
                self.table.setItem(
                    r, 0, QTableWidgetItem(str(data.get("atom_idx", "")))
                )
                self.table.setItem(r, 1, QTableWidgetItem(data.get("atom_sym", "")))

                sigma_sample = data.get("shielding", 0.0)
                self.table.setItem(r, 2, QTableWidgetItem(f"{sigma_sample:.2f}"))

                # Chemical shift: δ = δ_ref + (σ_ref - σ_sample)
                delta_sample = self.delta_ref + (self.sigma_ref - sigma_sample)
                self.table.setItem(r, 3, QTableWidgetItem(f"{delta_sample:.2f}"))

                # J-Coupling
                j_str = self.get_j_coupling_string([data.get("atom_idx", -1)])
                self.table.setItem(r, 4, QTableWidgetItem(j_str))

        # Auto enable couplings logic moved here (before plot) to ensure graph is correct immediately
        if getattr(self, "chk_real_spectrum", None) is not None:
            has_relevant_coupling = False
            if self.couplings:
                current_indices = {d["atom_idx"] for d in self.displayed_data}
                # Quick check if any displayed atom is involved in a coupling
                # Optimize: check against set of coupled indices if possible, or loop
                for c in self.couplings:
                    if (
                        c["atom_idx1"] in current_indices
                        or c["atom_idx2"] in current_indices
                    ):
                        has_relevant_coupling = True
                        break

            # Update UI state
            self.chk_real_spectrum.blockSignals(True)

            # Check if reference exists for this nucleus
            nucleus_key = self.get_nucleus_key(self.current_nucleus)
            has_reference = nucleus_key in self.reference_standards

            # Request: For "All", disable coupling (it's confusing/invalid to mix them)
            # ALSO: If reference does not exist, disable (meaningless ppm)
            if self.current_nucleus == "All":
                self.chk_real_spectrum.setChecked(False)
                self.chk_real_spectrum.setEnabled(False)
                self.chk_real_spectrum.setToolTip(
                    "Coupling simulation disabled for 'All' view"
                )
            elif not has_reference:
                self.chk_real_spectrum.setChecked(False)
                self.chk_real_spectrum.setEnabled(False)
                self.chk_real_spectrum.setToolTip(
                    "No reference standard available for this nucleus"
                )
            elif has_relevant_coupling:
                self.chk_real_spectrum.setEnabled(True)
                self.chk_real_spectrum.setToolTip(
                    "Simulate multiplets using J-couplings (requires nmrsim)"
                )
                # Auto-check logic removed from here to prevent forcing it on every redraw
            else:
                self.chk_real_spectrum.setChecked(False)
                self.chk_real_spectrum.setEnabled(False)
                self.chk_real_spectrum.setToolTip(
                    "No coupling information was found for displayed atoms"
                )

            self.chk_real_spectrum.blockSignals(False)

            # Ensure spinboxes match state
            self.toggle_simulation_controls()

        self.plot_spectrum()

    def plot_spectrum(self):
        """NMRのスティックスペクトルを描画"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        # シミュレーションモードなら別メソッドへ
        if (
            getattr(self, "chk_real_spectrum", None) is not None
            and self.chk_real_spectrum.isChecked()
        ):
            self.plot_real_spectrum(ax)
            return

        # 共通ロジックからピークを取得
        self.peaks_metadata = self._get_current_peaks()

        if not self.peaks_metadata:
            ax.text(
                0.5,
                0.5,
                "No data to display",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=14,
            )
            self.canvas.draw_idle()
            return

        # リストから描画用の値を抽出（NameErrorを解決）
        shifts = [p[0] for p in self.peaks_metadata]
        intensities = [p[1] for p in self.peaks_metadata]

        # クリック判定用に保存
        self.current_shifts = shifts

        # X軸の範囲とパディング
        min_shift, max_shift = min(shifts), max(shifts)
        padding = (
            max(0.5, (max_shift - min_shift) * 0.15) if max_shift > min_shift else 5.0
        )

        # 描画
        markerline, stemlines, baseline = ax.stem(
            shifts, intensities, linefmt="b-", markerfmt="None", basefmt="k-"
        )
        stemlines.set_linewidth(2.5)
        stemlines.set_alpha(0.8)
        baseline.set_alpha(0.3)

        # 軸の設定
        if self.chk_auto_x.isChecked():
            ax.set_xlim(max_shift + padding, min_shift - padding)
        else:
            ax.set_xlim(self.spin_x_max.value(), self.spin_x_min.value())

        ax.set_ylim(0, max(intensities) * 1.2)
        ax.set_xlabel("Chemical Shift δ (ppm)", fontsize=10, fontweight="bold")
        ax.set_ylabel(
            f"{self.current_nucleus} Count"
            if self.current_nucleus != "All"
            else "Atom Count",
            fontsize=10,
            fontweight="bold",
        )
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(True, alpha=0.2, linestyle="--")

        # タイトル
        display_nuc = self.get_nucleus_key(self.current_nucleus)
        ax.set_title(
            f"{display_nuc} NMR Stick Spectrum", fontsize=12, fontweight="bold", pad=20
        )

        # 選択中ハイライトの再描画
        if self.selected_peak_indices:
            self.highlight_selected_peaks()

        self.figure.tight_layout()
        self.canvas.draw_idle()

    def highlight_selected_peaks(self):
        """Add red highlights and labels to selected peaks"""
        if not getattr(self, "current_shifts", None):
            # Clear highlights if no selection
            for artist in self.highlight_artists:
                try:
                    artist.remove()
                except Exception as _e:
                    logging.warning("[nmr_analysis.py:1463] silenced: %s", _e)
            self.highlight_artists = []
            self.canvas.draw_idle()
            return

        ax = self.figure.axes[0] if self.figure.axes else None
        if not ax:
            return

        # Clear old highlights and labels
        for artist in self.highlight_artists:
            try:
                artist.remove()
            except Exception as _e:
                logging.warning("[nmr_analysis.py:1477] silenced: %s", _e)
        self.highlight_artists = []

        # Add red highlights and text labels for selected peaks
        for idx in self.selected_peak_indices:
            if idx < len(self.current_shifts):
                shift = self.current_shifts[idx]

                # Draw red line over selected peak (only if not in show_all_mode)
                if not self.show_all_mode:
                    line = ax.axvline(
                        shift,
                        ymin=0,
                        ymax=1,
                        color="red",
                        linewidth=3.5,
                        alpha=0.7,
                        zorder=10,
                    )
                    self.highlight_artists.append(line)

                # Add text label above the peak
                if idx < len(self.peaks_metadata):
                    # Get peak metadata
                    _, _, is_merged, atom_indices = self.peaks_metadata[idx]

                    # Build label text from all atoms in this peak
                    label_parts = []
                    for atom_idx in atom_indices:
                        atom_item = next(
                            (
                                d
                                for d in self.data
                                if d.get("atom_idx", None) == atom_idx
                            ),
                            None,
                        )
                        if atom_item:
                            atom_sym = atom_item.get("atom_sym", "?")
                            label_parts.append(f"{atom_sym}{atom_idx}")

                    label_text = ",".join(label_parts) if label_parts else "?"

                    # Position label above the peak
                    label_color = "black" if self.show_all_mode else "red"
                    text = ax.text(
                        shift,
                        1.1,
                        label_text,
                        ha="center",
                        va="bottom",
                        fontsize=10,
                        fontweight="bold",
                        color=label_color,
                        zorder=11,
                    )
                    self.highlight_artists.append(text)

        self.canvas.draw_idle()

    def _get_current_peaks(self):
        """
        現在の核種フィルタ、マージ設定、リファレンス値に基づき、
        プロットおよび解析用の共通ピークリストを生成する。
        Returns: List of (delta, intensity, is_merged, atom_indices)
        """
        # 原子インデックスからデータへの高速ルックアップ用マップ
        atom_map = {d["atom_idx"]: d for d in self.data}
        # 現在の表示対象原子のインデックスセット
        displayed_indices = {d["atom_idx"] for d in self.displayed_data}

        peaks = []
        processed_indices = set()

        # 1. マージされたグループの処理
        for group in self.merged_peaks:
            indices = group["indices"]
            # グループ内の原子が表示対象（現在の核種）に含まれているかチェック
            if any(idx in displayed_indices for idx in indices):
                total_sigma = 0.0
                valid_count = 0
                for idx in indices:
                    atom = atom_map.get(idx, None)
                    if atom:
                        total_sigma += atom.get("shielding", 0.0)
                        valid_count += 1

                if valid_count > 0:
                    avg_sigma = total_sigma / valid_count
                    avg_delta = self.delta_ref + (self.sigma_ref - avg_sigma)
                    # 強度はマージされた原子数（積分値）
                    peaks.append((avg_delta, float(len(indices)), True, indices))
                    processed_indices.update(indices)

        # 2. 個別（マージされていない）原子の処理
        for item in self.displayed_data:
            idx = item.get("atom_idx", None)
            if idx not in processed_indices:
                sigma = item.get("shielding", 0.0)
                delta = self.delta_ref + (self.sigma_ref - sigma)
                peaks.append((delta, 1.0, False, [idx]))

        # 化学シフトの降順（NMRの慣習）でソートしておくと管理しやすい
        return sorted(peaks, key=lambda x: x[0], reverse=True)

    def export_spectrum(self):
        """Export spectrum to file"""
        current_nucleus = self.current_nucleus
        default_path = get_default_export_path(
            self.file_path, suffix=f"_nmr_{current_nucleus}_spectrum", extension=".png"
        )

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Spectrum",
            default_path,
            "PNG Image (*.png);;PDF (*.pdf);;SVG (*.svg)",
        )

        if filename:
            try:
                self.figure.savefig(filename, dpi=300, bbox_inches="tight")
                if self.parent_dlg and hasattr(self.parent_dlg, "mw"):
                    self.parent_dlg.mw.statusBar().showMessage(
                        f"Spectrum exported to: {os.path.basename(filename)}", 5000
                    )
                else:
                    print(f"Spectrum exported to: {filename}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    def export_spectrum_csv(self):
        """Export spectrum data to CSV (Stick or Real)"""
        if not self.displayed_data:
            QMessageBox.warning(self, "No Data", "No spectrum data to export.")
            return

        default_path = get_default_export_path(
            self.file_path,
            suffix=f"_nmr_{self.current_nucleus}_sticks",
            extension=".csv",
        )
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Spectrum CSV", default_path, "CSV Files (*.csv)"
        )
        if not filename:
            return

        try:
            with open(filename, "w") as f:
                # Check mode
                is_real = (
                    getattr(self, "chk_real_spectrum", None) is not None
                    and self.chk_real_spectrum.isChecked()
                )

                if is_real:
                    # Export XY data from matplotlib line
                    ax = self.figure.axes[0] if self.figure.axes else None
                    if ax and ax.lines:
                        # The last line drawn is usually the spectrum curve or baseline?
                        # In plot_real_spectrum, we draw: ax.plot(x_hz_grid / spectrometer_freq, y_total, 'b-', linewidth=1.2)
                        # We should make sure we grab the right line.
                        # It is the only 'b-' line usually.

                        # Better approach: Re-calculate the data?
                        # Extracting from plot is hacky but consistent with what is seen.
                        # Let's try to extract x,y data from the first non-stem line.

                        line = None
                        for l in ax.lines:
                            # Stems are Line2D but usually handled differently.
                            # nmrsim plot is a standard plot.
                            if l.get_marker() == "None" and l.get_linestyle() == "-":
                                line = l
                                # If we have multiple, the 'blue' one is the spectrum.
                                if l.get_color() == "b":
                                    break

                        if line:
                            x_data = line.get_xdata()
                            y_data = line.get_ydata()

                            f.write("Chemical Shift (ppm),Intensity\n")
                            for x, y in zip(x_data, y_data):
                                f.write(f"{x:.6f},{y:.6f}\n")
                        else:
                            f.write("Error: Could not extract curve data.\n")
                    else:
                        f.write("Error: No plot found.\n")
                else:
                    # Stick Spectrum Export
                    # Export Peak List: Shift, Intensity
                    f.write("Chemical Shift (ppm),Intensity,AtomIndices\n")

                    # We can use peaks_metadata if available, or reconstruct
                    if (
                        getattr(self, "peaks_metadata", None) is not None
                        and self.peaks_metadata
                    ):
                        for (
                            shift,
                            intensity,
                            is_merged,
                            atom_indices,
                        ) in self.peaks_metadata:
                            indices_str = ";".join(str(i) for i in atom_indices)
                            f.write(f"{shift:.6f},{intensity:.4f},{indices_str}\n")
                    else:
                        # Fallback
                        f.write("# No peak data available.\n")

            if self.parent_dlg and hasattr(self.parent_dlg, "mw"):
                self.parent_dlg.mw.statusBar().showMessage(
                    f"Spectrum data exported to: {os.path.basename(filename)}", 5000
                )
            else:
                print(f"Spectrum data exported to: {filename}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    def export_table_csv(self):
        """Export table data to CSV"""
        current_nucleus = self.current_nucleus
        default_path = get_default_export_path(
            self.file_path, suffix=f"_nmr_{current_nucleus}_table", extension=".csv"
        )

        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Table CSV", default_path, "CSV Files (*.csv)"
        )
        if not filename:
            return

        try:
            with open(filename, "w") as f:
                # Headers
                headers = []
                for c in range(self.table.columnCount()):
                    headers.append(self.table.horizontalHeaderItem(c).text())
                f.write(",".join(headers) + "\n")

                # Rows
                for r in range(self.table.rowCount()):
                    cols = []
                    for c in range(self.table.columnCount()):
                        it = self.table.item(r, c)
                        text = it.text() if it else ""
                        # Escape commas if present
                        if "," in text:
                            text = f'"{text}"'
                        cols.append(text)
                    f.write(",".join(cols) + "\n")

            if self.parent_dlg and hasattr(self.parent_dlg, "mw"):
                self.parent_dlg.mw.statusBar().showMessage(
                    f"Table data exported to: {os.path.basename(filename)}", 5000
                )
            else:
                print(f"Table data exported to: {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    def get_j_coupling_string(self, atom_indices):
        """Format J-couplings for a list of atoms (merged or single)"""
        if not self.couplings:
            return ""

        # If merged group, we might have too many couplings.
        # Strategy:
        # 1. Collect all couplings involving ANY atom in atom_indices.
        # 2. Filter out couplings WITHIN the group (intra-group couplings often effectively 0 or infinite depending on equivalence,
        #    but usually we care about couplings to OUTSIDE).
        # 3. Format as "AtomSymIdx=J"
        # 4. If multiple atoms in group have SAME coupling to SAME partner (magnetic equivalence), show once.

        relevant_couplings = []  # (PartnerIdx, J)
        group_set = set(atom_indices)

        # Optimize: Pre-filter couplings? dataset is small usually.
        for c in self.couplings:
            idx1 = c["atom_idx1"]
            idx2 = c["atom_idx2"]
            j_val = abs(c["coupling"])

            if j_val < 0.1:
                continue  # Ignore tiny couplings

            partner = None
            if idx1 in group_set and idx2 not in group_set:
                partner = idx2
            elif idx2 in group_set and idx1 not in group_set:
                partner = idx1

            if partner is not None:
                relevant_couplings.append((partner, j_val))

        if not relevant_couplings:
            return ""

        # If merged group, we might see duplicates: (PartnerA, 7.5) from Atom1, (PartnerA, 7.5) from Atom2...
        # We should average them or check consistency.
        # Group by partner
        partner_map = {}
        for p, j in relevant_couplings:
            if p not in partner_map:
                partner_map[p] = []
            partner_map[p].append(j)

        # Format strings
        parts = []

        # Sort partners by Atom Index for consistency
        sorted_partners = sorted(partner_map.keys())

        for p in sorted_partners:
            # Find symbol for partner
            p_item = next((d for d in self.data if d.get("atom_idx", None) == p), None)
            p_sym = p_item.get("atom_sym", "") if p_item else ""
            p_label = f"{p_sym}{p}"

            j_list = partner_map[p]
            avg_j = sum(j_list) / len(j_list)

            # If standard deviation is high, maybe indicate?
            # For now, just show average.
            parts.append(f"{p_label}={avg_j:.1f}")

        return ", ".join(parts)

    def copy_table(self):
        """Copy table data to clipboard"""
        text = "Idx\tNucleus\tShielding\tShift\tJ-coupling\n"
        for r in range(self.table.rowCount()):
            cols = []
            for c in range(self.table.columnCount()):
                it = self.table.item(r, c)
                cols.append(it.text() if it else "")
            text += "\t".join(cols) + "\n"
        QApplication.clipboard().setText(text)
        if self.parent_dlg and hasattr(self.parent_dlg, "mw"):
            self.parent_dlg.mw.statusBar().showMessage(
                "Table data copied to clipboard!", 5000
            )
        else:
            print("Table data copied to clipboard!")

    def toggle_simulation_controls(self):
        """Enable/Disable simulation spinboxes based on checkbox state"""
        if getattr(self, "chk_real_spectrum", None) is None:
            return

        is_sim_active = (
            self.chk_real_spectrum.isChecked()
            and self.chk_real_spectrum.isEnabled()
            and nmrsim is not None
        )

        if getattr(self, "spin_real_width", None) is not None:
            self.spin_real_width.setEnabled(is_sim_active)
        if getattr(self, "spin_mhz", None) is not None:
            self.spin_mhz.setEnabled(is_sim_active)

    def on_peak_click(self, event):
        """Handle clicking on a peak in the spectrum"""
        click_x = event.xdata
        if click_x is None:
            return

        # Find nearest peak within tolerance (relative to current x-axis range)
        xlim = event.inaxes.get_xlim()
        x_range = abs(xlim[1] - xlim[0])
        tolerance = x_range * 0.01  # 1% of view width

        distances = [abs(shift - click_x) for shift in self.current_shifts]
        min_distance = min(distances)

        if min_distance > tolerance:
            return  # Click too far from any peak

        # Find the clicked peak index
        peak_idx = distances.index(min_distance)

        # Check for Shift or Ctrl key (using Qt modifiers)
        modifiers = QApplication.keyboardModifiers()
        is_multi = bool(
            modifiers
            & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier)
        )

        if peak_idx < len(self.current_shifts):
            if not is_multi:
                # Normal click
                if (
                    len(self.selected_peak_indices) == 1
                    and peak_idx in self.selected_peak_indices
                ):
                    # User clicked the ONLY selected peak -> Toggle OFF (Deselect)
                    self.selected_peak_indices = set()
                else:
                    # New peak, or switching from multi-selection -> Select ONLY this peak
                    self.selected_peak_indices = {peak_idx}
            else:
                # Shift or Ctrl + Click: Toggle selection
                if peak_idx in self.selected_peak_indices:
                    self.selected_peak_indices.remove(peak_idx)
                else:
                    self.selected_peak_indices.add(peak_idx)

            # Force a fresh draw on next poll if desired
            self._last_highlight_atoms = set()

            # Update highlights
            self.highlight_selected_peaks()

            # Update 3D labels for all selected peaks
            self.update_selected_labels()

    def update_selected_labels(self, is_external_sync=False):
        """Update 3D labels and spheres for all selected peaks"""
        # 1. Clear existing labels
        self.clear_atom_labels()

        # 2. Add labels for each selected peak
        all_peak_indices = set()

        for peak_idx in sorted(self.selected_peak_indices):
            if peak_idx < len(self.peaks_metadata):
                # Get metadata for this peak (shift, intensity, is_merged, atom_indices)
                _, _, is_merged, atom_indices = self.peaks_metadata[peak_idx]

                # Add label for each atom in this peak (handles both merged and individual)
                for atom_idx in atom_indices:
                    all_peak_indices.add(atom_idx)
                    # Find atom symbol from data
                    atom_item = next(
                        (d for d in self.data if d.get("atom_idx", None) == atom_idx),
                        None,
                    )
                    if atom_item:
                        atom_sym = atom_item.get("atom_sym", "?")
                        self.add_atom_label(atom_idx, atom_sym)

        # 3. Synchronize with Main Window
        if hasattr(self.parent_dlg, "mw"):
            mw = self.parent_dlg.mw

            # If this is a sync FROM 3D (user clicked in viewer), we should NOT clear the 3D selection!
            # We only clear it if the user clicked in the Graph (internal sync), to replace Green with Yellow.
            if not is_external_sync:
                # Ensure we don't have double spheres (Green + Yellow).
                # We clear the global selection so ONLY our internal Yellow spheres are visible.
                e3d = getattr(mw, "edit_3d_manager", None)
                if e3d:
                    e3d.selected_atoms_3d.clear()

                # CRITICAL: Update our sync tracker so the polling loop knows WE did this
                # and doesn't interpret the empty set as a user unselection.
                self._last_synced_mw_selection = frozenset()

                # Sync to MW if we are the originator (internal sync)
                e3d = getattr(mw, "edit_3d_manager", None)
                if e3d:
                    try:
                        e3d.update_3d_selection_display()
                    except Exception as _e:
                        logging.warning("[nmr_analysis.py:1880] silenced: %s", _e)

            # Draw yellow highlights for NMR selection
            self.draw_custom_nmr_highlights_3d(all_peak_indices)

            # Debug print to confirming highlighting path is taken
            # print(f"Highlighting {len(self.selected_peak_indices)} peaks with {len(atom_coords)} atoms")

            # Render once after all labels added
            v3d = getattr(self.parent_dlg.mw, "view_3d_manager", None)
            if v3d and hasattr(v3d, "plotter"):
                v3d.plotter.render()

    def add_atom_label(self, atom_idx, atom_sym):
        """Add a single atom label to 3D viewer"""
        # Check if parent has plotter
        v3d = (
            getattr(self.parent_dlg.mw, "view_3d_manager", None)
            if hasattr(self.parent_dlg, "mw")
            else None
        )
        if not v3d or not hasattr(v3d, "plotter"):
            return

        # Get coordinates
        coords = self.parent_dlg.parser.data.get("coords", [])
        if not coords or atom_idx >= len(coords):
            return

        try:
            pos = coords[atom_idx]
            label_pos = [pos[0], pos[1], pos[2] + 0.4]  # Offset above atom

            label_name = f"nmr_label_{atom_idx}"
            actor = v3d.plotter.add_point_labels(
                [label_pos],
                [f"{atom_sym}{atom_idx}"],
                font_size=12,
                text_color="cyan",
                point_size=0,
                always_visible=True,
                bold=True,
                name=label_name,
            )
            self._atom_labels.append(actor)
            self._nmr_label_names.append(label_name)
        except Exception as e:
            logging.warning("[nmr_analysis.py:1924] silenced: %s", e)

    def highlight_atom_in_3d(self, atom_idx, atom_sym):
        """Highlight selected atom with a label in 3D viewer (legacy - now uses update_selected_labels)"""
        # This is now handled by update_selected_labels

    def clear_peak_selection(self):
        """Clear all selected peaks and their labels"""
        # Clear selected peaks
        self.selected_peak_indices.clear()

        # Clear highlights
        for artist in self.highlight_artists:
            try:
                artist.remove()
            except Exception as _e:
                logging.warning("[nmr_analysis.py:1941] silenced: %s", _e)
        self.highlight_artists = []

        # Clear 3D labels
        self.clear_atom_labels()

        # [Commented out to avoid doubled spheres]
        if hasattr(self.parent_dlg, "mw"):
            mw = self.parent_dlg.mw
            e3d = getattr(mw, "edit_3d_manager", None)
            if e3d:
                e3d.selected_atoms_3d.clear()
                try:
                    e3d.update_3d_selection_display()
                except Exception as _e:
                    logging.warning("[nmr_analysis.py:1956] silenced: %s", _e)

        # Redraw spectrum
        if getattr(self, "canvas", None) is not None:
            self.canvas.draw_idle()

    def plot_real_spectrum(self, ax):
        """Jカップリングを考慮したスペクトル描画（共通ロジック完全統合版）"""

        # 1. 【共通ロジック】からピーク情報を取得し、既存コードが期待する変数名に代入
        peaks_to_simulate = self._get_current_peaks()

        if not peaks_to_simulate:
            ax.text(
                0.5,
                0.5,
                "No data to simulate",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            self.canvas.draw_idle()
            return

        # 2. クラス変数にも保存（クリック判定やハイライト用）
        self.peaks_metadata = peaks_to_simulate
        self.current_shifts = [p[0] for p in peaks_to_simulate]

        # 3. 核種シンボルの定義
        target_nuc_sym = self.current_nucleus
        if target_nuc_sym == "All":
            if self.displayed_data:
                target_nuc_sym = self.displayed_data[0].get("atom_sym", "H")
            else:
                target_nuc_sym = "H"

        # 4. 周波数と磁気回転比の計算
        base_freq_mhz = (
            self.spin_mhz.value()
            if getattr(self, "spin_mhz", None) is not None
            else 400.0
        )

        ratio = 1.0
        lookup_sym = target_nuc_sym.upper()

        import re

        element_only = re.sub(r"[^A-Z]", "", lookup_sym)

        if lookup_sym in self.GAMMA:
            ratio = abs(self.GAMMA[lookup_sym]) / self.GAMMA["H"]
        elif element_only in self.GAMMA:
            ratio = abs(self.GAMMA[element_only]) / self.GAMMA["H"]

        spectrometer_freq = base_freq_mhz * ratio
        points = 8192
        width_hz = (
            max(0.1, self.spin_real_width.value())
            if getattr(self, "spin_real_width", None) is not None
            else 0.5
        )

        # --- Jカップリングの計算とMultiplet生成 ---
        atom_to_group_size = {
            idx: len(g["indices"]) for g in self.merged_peaks for idx in g["indices"]
        }
        all_multiplets = []

        # 5. ここからループ開始（peaks_to_simulate が定義されているのでエラーになりません）
        for shift, intensity, is_merged, atom_indices in peaks_to_simulate:
            peak_atoms_set = set(atom_indices)
            rep_item = next(
                (d for d in self.data if d.get("atom_idx", None) == atom_indices[0]),
                None,
            )
            peak_nuc = rep_item.get("atom_sym", "").strip().upper() if rep_item else ""

            couplings_list = []

            # カップリング計算ロジック
            if (
                getattr(self, "chk_real_spectrum", None) is not None
                and self.chk_real_spectrum.isChecked()
            ):
                partner_j_values = {}
                partner_multiplicity = {}

                for atom_idx in atom_indices:
                    for c in self.couplings:
                        other = (
                            c["atom_idx2"]
                            if c["atom_idx1"] == atom_idx
                            else (
                                c["atom_idx1"] if c["atom_idx2"] == atom_idx else None
                            )
                        )
                        if other is not None and other not in peak_atoms_set:
                            other_item = next(
                                (d for d in self.data if d["atom_idx"] == other), None
                            )
                            if (
                                other_item
                                and other_item.get("atom_sym", "").strip().upper()
                                == peak_nuc
                            ):
                                J = abs(c["coupling"])
                                if J > 0.01:
                                    pid = next(
                                        (
                                            min(g["indices"])
                                            for g in self.merged_peaks
                                            if other in g["indices"]
                                        ),
                                        other,
                                    )
                                    if pid not in partner_j_values:
                                        partner_j_values[pid] = []
                                    partner_j_values[pid].append(J)
                                    partner_multiplicity[pid] = atom_to_group_size.get(
                                        other, 1
                                    )

                for pid, j_list in partner_j_values.items():
                    n_partner = partner_multiplicity[pid]
                    n_current = len(atom_indices)
                    total_connections = n_current * n_partner
                    avg_J = sum(j_list) / total_connections
                    if avg_J > 0.1:
                        couplings_list.append((avg_J, n_partner))

            if Multiplet:
                try:
                    m = Multiplet(shift * spectrometer_freq, intensity, couplings_list)
                    m.w = width_hz
                    all_multiplets.append(m)
                except Exception as e:
                    logging.error(
                        "NMR: Multiplet creation failed for shift=%.3f: %s", shift, e
                    )

        # --- プロット範囲計算 ---
        if self.chk_auto_x.isChecked():
            shifts_all_sim = [p[0] for p in peaks_to_simulate]
            min_s, max_s = min(shifts_all_sim), max(shifts_all_sim)
            padding = max(0.5, (max_s - min_s) * 0.15) if max_s > min_s else 5.0
            x_limit_max, x_limit_min = max_s + padding, min_s - padding
        else:
            x_limit_max, x_limit_min = self.spin_x_max.value(), self.spin_x_min.value()

        # シミュレーション範囲の計算
        l_low = min(x_limit_min, x_limit_max) * spectrometer_freq
        l_high = max(x_limit_min, x_limit_max) * spectrometer_freq
        if l_low > l_high:
            l_low, l_high = l_high, l_low
        span = l_high - l_low
        if span == 0:
            span = 100
        l_low_sim = l_low - span * 0.1
        l_high_sim = l_high + span * 0.1

        x_hz_grid = np.linspace(l_low_sim, l_high_sim, points)
        y_total = np.zeros_like(x_hz_grid)
        gamma = width_hz / 2.0

        # --- nmrsim計算またはフォールバック ---
        nmrsim_success = False
        if all_multiplets and Spectrum:
            try:
                spec = Spectrum(all_multiplets)
                x_sim, y_sim = spec.lineshape(points=points)
                y_total = np.interp(x_hz_grid, x_sim, y_sim, left=0, right=0)
                if np.max(y_total) >= 1e-9:
                    nmrsim_success = True
            except Exception as e:
                logging.error("NMR: nmrsim Spectrum simulation failed: %s", e)

        if not nmrsim_success:
            all_peaks = []
            for p in peaks_to_simulate:
                freq = p[0] * spectrometer_freq
                intensity = p[1]
                all_peaks.append((freq, intensity))

            if all_peaks:
                vs = np.array([p[0] for p in all_peaks])
                is_ = np.array([p[1] for p in all_peaks])
                chunk_size = 100
                for i in range(0, len(vs), chunk_size):
                    v_chunk = vs[i : i + chunk_size]
                    i_chunk = is_[i : i + chunk_size]
                    for v, I in zip(v_chunk, i_chunk):
                        y_total += I * (
                            gamma / (np.pi * ((x_hz_grid - v) ** 2 + gamma**2))
                        )

        # 高さの正規化
        max_p = max(p[1] for p in peaks_to_simulate) if peaks_to_simulate else 1.0
        if np.max(y_total) > 0:
            y_total = y_total / np.max(y_total) * max_p

        # プロット
        ax.plot(x_hz_grid / spectrometer_freq, y_total, "b-", linewidth=1.2)
        if np.max(y_total) > 0:
            ax.set_ylim(0, np.max(y_total) * 1.3)

        # スタイル設定
        ax.set_xlim(x_limit_max, x_limit_min)
        ax.set_xlabel("Chemical Shift δ (ppm)", fontsize=10, fontweight="bold")
        ax.set_ylabel("Intensity", fontsize=10, fontweight="bold")
        ax.tick_params(axis="both", which="major", labelsize=8)

        # タイトル設定
        current_nucleus_title = self.current_nucleus
        display_nuc = self.get_nucleus_key(current_nucleus_title)
        ref_name_title = (
            self.combo_ref.currentText()
            if getattr(self, "combo_ref", None) is not None
            else "Custom"
        )

        ax.set_title(
            f"{display_nuc} NMR Spectrum", fontsize=12, fontweight="bold", pad=20
        )
        ref_text = f"Ref: {ref_name_title} (δ_ref = {self.delta_ref:.2f} ppm, σ_ref = {self.sigma_ref:.1f} ppm)"
        ax.text(
            0.5,
            1.02,
            ref_text,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=9,
        )

        ax.yaxis.set_visible(True)
        ax.grid(True, alpha=0.25, linestyle="--", axis="x", linewidth=0.8)
        ax.grid(True, alpha=0.15, linestyle=":", axis="y", linewidth=0.5)

        self.figure.tight_layout()

        # --- 修正3: イベント接続とハイライト更新 ---
        # クリックイベントを再接続
        # self.canvas.mpl_connect('button_press_event', self.on_peak_click)

        # 選択済みのピークがあればハイライト（赤い線やラベル）を再描画
        if self.selected_peak_indices:
            self.highlight_selected_peaks()

        self.canvas.draw_idle()

    def toggle_all_labels(self):
        """Toggle showing all atom labels on the spectrum graph"""
        show_all = self.chk_show_all_labels.isChecked()

        if show_all:
            # Enable show all mode (labels without red highlights)
            self.show_all_mode = True
            self.selected_peak_indices.clear()

            # Select all peaks to show labels
            # Use peaks_metadata length if available, otherwise displayed_data as fallback
            count = (
                len(self.peaks_metadata)
                if getattr(self, "peaks_metadata", None) is not None
                else len(self.displayed_data)
            )
            for i in range(count):
                self.selected_peak_indices.add(i)

            # Update graph with labels only (no red highlights)
            self.highlight_selected_peaks()

            # For "Show All", we also want to update the 3D view to reflect "All"
            # or at least clear the specific "red" selection we had.
            # If we want to show labels for ALL atoms in 3D:
            self.update_selected_labels()
        else:
            # Disable show all mode
            self.show_all_mode = False
            # Clear all selections (graph and 3D)
            self.clear_peak_selection()

    def clear_atom_labels(self):
        """Remove all atom labels and custom selection spheres from 3D viewer"""
        v3d = (
            getattr(self.parent_dlg.mw, "view_3d_manager", None)
            if hasattr(self.parent_dlg, "mw")
            else None
        )
        if not v3d or not hasattr(v3d, "plotter"):
            return

        plotter = v3d.plotter

        # 1. Clear custom NMR selection spheres by name (most reliable in PyVista)
        try:
            plotter.remove_actor("nmr_selection_highlights")
        except Exception as _e:
            logging.warning("[nmr_analysis.py:2182] silenced: %s", _e)

        # 2. Clear labels by tracked name
        if getattr(self, "_nmr_label_names", None) is not None:
            for name in self._nmr_label_names:
                try:
                    plotter.remove_actor(name)
                except Exception as _e:
                    logging.warning("[nmr_analysis.py:2190] silenced: %s", _e)
            self._nmr_label_names = []

        # 3. Fallback: Clear labels by list reference
        for actor in self._atom_labels:
            try:
                plotter.remove_actor(actor)
            except Exception as _e:
                logging.warning("[nmr_analysis.py:2198] silenced: %s", _e)
        self._atom_labels = []

        # 4. Clean up private spheres actor list
        if getattr(self, "_nmr_sphere_actors", None) is not None:
            for actor in self._nmr_sphere_actors:
                try:
                    plotter.remove_actor(actor)
                except Exception as _e:
                    logging.warning("[nmr_analysis.py:2207] silenced: %s", _e)
            self._nmr_sphere_actors = []

        try:
            plotter.render()
        except Exception as _e:
            logging.warning("[nmr_analysis.py:2213] silenced: %s", _e)

    def draw_custom_nmr_highlights_3d(self, atom_indices):
        """Draw yellow highlight spheres for selected atoms in 3D viewer"""
        mw = self.parent_dlg.mw if hasattr(self.parent_dlg, "mw") else None
        v3d = getattr(mw, "view_3d_manager", None) if mw else None
        if not v3d or not hasattr(v3d, "plotter"):
            return

        plotter = v3d.plotter

        # ALWAYS clear existing custom highlights first to prevent stacking/phantom spheres
        try:
            plotter.remove_actor("nmr_selection_highlights")
        except Exception as _e:
            logging.warning("[nmr_analysis.py:2228] silenced: %s", _e)

        # Clear tracker list to prevent stale references
        self._nmr_sphere_actors = []

        # If no indices provided, just render the cleared state and return
        if not atom_indices or not hasattr(v3d, "atom_positions_3d"):
            try:
                plotter.render()
            except Exception as _e:
                logging.warning("[nmr_analysis.py:2238] silenced: %s", _e)
            return

        indices = list(atom_indices)
        valid_indices = [
            i for i in indices if i < len(mw.view_3d_manager.atom_positions_3d)
        ]

        if not valid_indices:
            return

        try:
            # Get positions
            selected_positions = mw.view_3d_manager.atom_positions_3d[valid_indices]

            # Highlight sphere size: 40% (1.4x) relative to VDW radii per user request
            radii = []
            for i in valid_indices:
                # Find atom symbol from parser data
                atom_item = next(
                    (d for d in self.data if i == d.get("atom_idx", None)), None
                )
                sym = atom_item.get("atom_sym", "C") if atom_item else "C"
                # Use 1.4x scaling factor (40% larger)
                r = VDW_RADII.get(sym, 0.4) * 1.4
                radii.append(r)

            # Create glyphs for highlights
            highlight_source = pv.PolyData(selected_positions)
            highlight_source["radii"] = np.array(radii)

            highlight_glyphs = highlight_source.glyph(
                scale="radii",
                geom=pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16),
                orient=False,
            )

            # Add to plotter and track actor
            actor = plotter.add_mesh(
                highlight_glyphs,
                color="yellow",
                opacity=0.3,
                name="nmr_selection_highlights",
            )
            self._nmr_sphere_actors.append(actor)
            plotter.render()

        except Exception as e:
            logging.warning("[nmr_analysis.py:2282] silenced: %s", e)

    def closeEvent(self, event):
        """Clean up labels when dialog closes"""
        self.save_settings()
        self.clear_atom_labels()
        super().closeEvent(event)

    def keyPressEvent(self, event):
        """Intercept Enter/Return keys to prevent dialog reset or closure"""
        if event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            # Focus on parent to clear focus from spinboxes but don't close
            self.setFocus()
            event.accept()
            return
        super().keyPressEvent(event)

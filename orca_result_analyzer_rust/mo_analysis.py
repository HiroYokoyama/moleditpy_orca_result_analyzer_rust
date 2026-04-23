import os
import tempfile
import numpy as np
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QTreeWidget,
    QTreeWidgetItem,
    QAbstractItemView,
    QMessageBox,
    QFileDialog,
    QProgressDialog,
    QHeaderView,
    QGroupBox,
    QSpinBox,
    QDoubleSpinBox,
    QFormLayout,
    QTreeWidgetItemIterator,
    QApplication,
    QColorDialog,
    QInputDialog,
    QComboBox,
    QCheckBox,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QBrush
import json
from .utils import get_default_export_path
import logging

try:
    from .mo_engine import BasisSetEngine, CalcWorker
    from .vis import CubeVisualizer
except ImportError:
    try:
        from mo_engine import BasisSetEngine, CalcWorker
        from vis import CubeVisualizer
    except Exception:
        BasisSetEngine = None
        CalcWorker = None
        CubeVisualizer = None

try:
    from .energy_diag import EnergyDiagramDialog
except ImportError:
    EnergyDiagramDialog = None


class MODialog(QDialog):
    def __init__(self, parent, mo_data, result_dir=None):
        super().__init__(parent)
        self.mw = None
        if hasattr(parent, "mw"):
            self.mw = parent.mw
        elif hasattr(parent, "context"):
            self.mw = parent.context.get_main_window()

        self.mo_data = mo_data
        self.setWindowTitle("MO Analysis & Visualization")
        self.resize(550, 750)
        self.mos = mo_data
        self.parent_dlg = parent
        self.last_cube_path = None
        self.generation_queue = []  # Init queue
        self.energy_dlg = None  # Track Energy Diagram
        self.setup_ui()

    def get_cube_path(self, display_id):
        if not hasattr(self.parent_dlg, "parser") or not self.parent_dlg.parser:
            return None

        parser = self.parent_dlg.parser
        if hasattr(parser, "filename") and parser.filename:
            fpath = parser.filename
            base_dir = os.path.dirname(fpath)
            filename_base = os.path.splitext(os.path.basename(fpath))[0]
            out_dir = os.path.join(base_dir, f"{filename_base}_cubes")

            # Just return the path, don't create dir here (lazy creation requested)
            # Sanitize display ID (which might contain "1 (a)")
            # Example: "1 (a)" -> "1_a"
            safe_id = (
                str(display_id)
                .replace(" ", "_")
                .replace("(", "")
                .replace(")", "")
                .replace(":", "")
            )
            # Ensure no double underscores
            while "__" in safe_id:
                safe_id = safe_id.replace("__", "_")

            # User requested prefix: benzene-ene_MO_... (Using filename base)
            return os.path.join(out_dir, f"{filename_base}_MO_{safe_id}.cube")
        return None

    def setup_ui(self):
        # Use simpler Vertical Layout to fill the window
        layout = QVBoxLayout(self)

        # 1. Visualization Settings (Grouping "Grid" vs "Calc")
        vis_grp = QGroupBox("Visualization Controls")
        vis_layout = QVBoxLayout(vis_grp)

        # --- Presets ---
        h_preset = QHBoxLayout()
        h_preset.addWidget(QLabel("Preset:"))
        self.combo_presets = QComboBox()
        self.combo_presets.addItem("Default")
        h_preset.addWidget(self.combo_presets)

        btn_save_preset = QPushButton("Save")
        btn_save_preset.setFixedWidth(50)
        btn_save_preset.clicked.connect(self.save_preset)
        h_preset.addWidget(btn_save_preset)

        btn_del_preset = QPushButton("Del")
        btn_del_preset.setFixedWidth(40)
        btn_del_preset.clicked.connect(self.delete_preset)
        h_preset.addWidget(btn_del_preset)
        vis_layout.addLayout(h_preset)

        # --- Style & Colors ---
        h_style = QHBoxLayout()
        h_style.addWidget(QLabel("Style:"))
        self.combo_style = QComboBox()
        self.combo_style.addItems(["Surface", "Wireframe", "Points"])
        h_style.addWidget(self.combo_style)

        self.check_smooth = QCheckBox("Smooth Shading")
        self.check_smooth.setChecked(True)
        h_style.addWidget(self.check_smooth)
        vis_layout.addLayout(h_style)

        h_colors = QHBoxLayout()
        self.btn_color_p = QPushButton("Pos (+)")
        self.btn_color_p.setStyleSheet(
            "background-color: red; color: white; font-weight: bold;"
        )
        self.btn_color_p.clicked.connect(lambda: self.pick_color("p"))
        h_colors.addWidget(self.btn_color_p)

        self.btn_color_n = QPushButton("Neg (-)")
        self.btn_color_n.setStyleSheet(
            "background-color: blue; color: white; font-weight: bold;"
        )
        self.btn_color_n.clicked.connect(lambda: self.pick_color("n"))
        h_colors.addWidget(self.btn_color_n)
        vis_layout.addLayout(h_colors)

        # Isovalue / Opacity
        h_iso = QHBoxLayout()
        h_iso.addWidget(QLabel("Isovalue:"))
        self.spin_iso = QDoubleSpinBox()
        self.spin_iso.setRange(0.001, 1.0)
        self.spin_iso.setSingleStep(0.005)
        self.spin_iso.setDecimals(3)
        self.spin_iso.setValue(0.02)
        h_iso.addWidget(self.spin_iso)

        h_iso.addWidget(QLabel("Opacity:"))
        self.spin_opacity = QDoubleSpinBox()
        self.spin_opacity.setRange(0.0, 1.0)
        self.spin_opacity.setSingleStep(0.1)
        self.spin_opacity.setValue(0.5)
        h_iso.addWidget(self.spin_opacity)
        vis_layout.addLayout(h_iso)

        # Connect changes to update view immediately if exists
        self.combo_style.currentTextChanged.connect(self.update_vis_only)
        self.check_smooth.toggled.connect(self.update_vis_only)
        self.spin_iso.valueChanged.connect(self.update_vis_only)
        self.spin_opacity.valueChanged.connect(self.update_vis_only)

        # Load Settings must be called after UI init
        self.load_settings()
        self.combo_presets.currentTextChanged.connect(self.apply_preset)

        # Calculation Settings
        calc_grp = QGroupBox("Calculation Parameters")
        calc_layout = QFormLayout(calc_grp)

        self.spin_pts = QSpinBox()
        self.spin_pts.setRange(10, 200)
        self.spin_pts.setValue(40)
        self.spin_pts.setSuffix(" pts")
        calc_layout.addRow("Grid Resolution (x,y,z):", self.spin_pts)

        self.spin_margin = QDoubleSpinBox()
        self.spin_margin.setRange(1.0, 15.0)
        self.spin_margin.setValue(4.0)
        self.spin_margin.setSuffix(" Bohr")  # Correct unit
        calc_layout.addRow("Calc Boundary Margin:", self.spin_margin)

        vis_layout.addWidget(calc_grp)

        # Warning Layout (Label + Copy Button)
        warn_layout = QHBoxLayout()

        self.lbl_warning = QLabel("")
        self.lbl_warning.setStyleSheet(
            "color: #e65100; font-weight: bold; font-size: 9pt;"
        )
        self.lbl_warning.setWordWrap(True)
        self.lbl_warning.setVisible(False)
        warn_layout.addWidget(self.lbl_warning)

        self.btn_copy_input = QPushButton("Copy Input")
        self.btn_copy_input.setFixedWidth(80)
        self.btn_copy_input.setVisible(False)
        self.btn_copy_input.clicked.connect(self.copy_orca_input)
        warn_layout.addWidget(self.btn_copy_input)

        vis_layout.addLayout(warn_layout)

        layout.addWidget(vis_grp)

        # 2. MO List
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["MO", "Label", "Occ", "E (eV)", "E (Eh)"])
        # Columns: Let them fill
        self.tree.setColumnWidth(0, 120)
        self.tree.setColumnWidth(1, 100)  # Label
        self.tree.setColumnWidth(2, 80)
        self.tree.setColumnWidth(3, 80)
        # Last column stretches
        header = self.tree.header()
        if header:
            header.setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)

        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.tree.itemDoubleClicked.connect(self.on_double_click)
        # Connect currentItemChanged for both keyboard navigation and single click
        self.tree.currentItemChanged.connect(self.on_item_changed)
        # Keep selection changed only for button state specific logic if needed,
        # but currentItemChanged covers the "primary" selection change.
        self.tree.itemSelectionChanged.connect(self.on_selection_changed)
        self.tree.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(self.tree)

        # 3. Action Buttons
        btn_layout = QHBoxLayout()
        self.btn_vis = QPushButton("Visualize Selected")
        self.btn_vis.setStyleSheet("font-weight: bold; background-color: #d0f0c0;")
        self.btn_vis.clicked.connect(self.visualize_selected_mos)
        self.btn_vis.setEnabled(False)  # Default disabled until selection
        btn_layout.addWidget(self.btn_vis)

        # Connect selection
        self.tree.itemSelectionChanged.connect(self.on_selection_changed)

        # MO Diagram Button
        btn_diag = QPushButton("Show MO Diagram")
        btn_diag.clicked.connect(self.show_mo_diagram)
        btn_layout.addWidget(btn_diag)

        btn_csv = QPushButton("Export CSV")
        btn_csv.clicked.connect(self.export_csv)
        btn_layout.addWidget(btn_csv)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btn_layout.addWidget(btn_close)

        # Add Stretch to push buttons to right or keep centered?
        # Standard dialogs often have buttons on right.
        # But centering is fine.

        # IMPORTANT: Add the button layout to the main layout!
        layout.addLayout(btn_layout)
        # If we want interaction with 3D view, we need Non-Modal (show()).

        # For now, we keep it Modal, but 'Visualize' updates the background window?
        # If modal, user can't rotate 3D view easily without closing dialog.
        # Let's add a "Apply/Update" button or just let it update.

        # To make it better:
        # We'll just show info on right side.
        # info_panel = QWidget()
        # info_layout = QVBoxLayout(info_panel)
        # info_layout.addWidget(QLabel("<b>Selected MO Info</b>"))
        # self.lbl_info = QLabel("Select an MO to view details.")
        # self.lbl_info.setWordWrap(True)
        # info_layout.addWidget(self.lbl_info)
        # info_layout.addStretch()
        # main_layout.addWidget(info_panel)

        # Populate
        self.normalize_and_populate()

    def normalize_and_populate(self):
        self.tree.clear()  # Verify clear first
        self.mo_list = []

        # Sort keys logic:
        # Keys are likely strings "0_alpha", "1_alpha" etc. or integers (if old/other format)
        # We want to sort by Spin then Index, or just by Energy if available?
        # Standard: Sort by Spin (Alpha first), then Index.

        keys = []
        if isinstance(self.mos, dict):
            raw_keys = list(self.mos.keys())
            # Check key type
            # print(f"DEBUG Keys: {raw_keys[:10]}") # Debugging
            if raw_keys and isinstance(raw_keys[0], str) and "_" in raw_keys[0]:
                # Sort by spin, then index -> NOW Sort by Index, then Spin (interleaved)
                def sort_key(k):
                    if "_" in str(k):
                        idx_s, spin_s = str(k).split("_")

                        if spin_s == "alpha":
                            pass
                        else:
                            pass  # beta comes first in sorted list -> appears later in reverse iteration -> appears LOWER in UI?

                        # So we need CUSTOM sort where Beta < Alpha.
                        s_priority = 0 if spin_s == "beta" else 1

                        return (int(idx_s), s_priority)
                    return (0, 0)

                keys = sorted(raw_keys, key=sort_key)
            else:
                # Fallback numeric sort
                keys = sorted(
                    raw_keys, key=lambda x: int(x) if str(x).isdigit() else str(x)
                )

            for k in keys:
                val = self.mos[k]
                if "id" not in val:
                    val["id"] = k
                # Store the Access Key used in self.mos dictionary
                val["_access_key"] = k
                self.mo_list.append(val)
        elif isinstance(self.mos, list):
            self.mo_list = self.mos
            # Add index as access key
            for i, mo in enumerate(self.mo_list):
                mo["_access_key"] = i

        # Separate HOMO/LUMO logic for each spin
        spin_mos = {"restricted": [], "alpha": [], "beta": []}
        for mo in self.mo_list:
            s = mo.get("spin", "restricted")
            if s not in spin_mos:
                s = "restricted"  # fallback using restricted logic for unknown strings if any
            spin_mos[s].append(mo)

        spin_homo_idx = {}
        for s in spin_mos:
            n_occ = 0
            for mo in spin_mos[s]:
                occ = mo.get("occ", mo.get("occupation", 0.0))
                if occ > 0.1:
                    n_occ += 1
            spin_homo_idx[s] = n_occ - 1

        # Display order: Reversed (High Energy Top)
        # We need to sort global list carefully if mixing spins?
        # Or just show them in the order we sorted above (Alpha then Beta)
        # Our self.mo_list is already sorted by Spin, Index.
        # But we iterate backwards! so Beta (high idx -> low idx) -> Alpha (high idx -> low idx)

        count = len(self.mo_list)
        for i in range(count - 1, -1, -1):
            mo = self.mo_list[i]

            # Internal unique key
            access_key = mo.get("_access_key", None)

            # Display ID
            mo_idx_val = mo.get("id", mo.get("index", i))
            spin = mo.get("spin", "restricted")

            # Determine if HOMO/LUMO for this specific spin channel
            # We need the index OF THIS MO within its spin group
            # Since self.mo_list is sorted, we can find it?
            # Or just use the ID? Usually ID corresponds to index.
            # Let's use the ID if available, assuming 0-based index.

            local_idx = -1
            try:
                local_idx = int(mo_idx_val)
            except Exception as _e:
                logging.warning("[mo_analysis.py:348] silenced: %s", _e)

            if spin in spin_homo_idx:
                h = spin_homo_idx[spin]

            display_id = str(mo_idx_val)

            label_id = display_id
            if spin == "alpha":
                label_id += " (a)"
            elif spin == "beta":
                label_id += " (b)"

            # Label Separation
            homo_lumo_label = ""
            if spin in spin_homo_idx:
                diff = local_idx - h

                if diff == 0:
                    homo_lumo_label = "HOMO"
                elif diff == 1:
                    homo_lumo_label = "LUMO"
                elif diff > 1 and diff < 100:
                    homo_lumo_label = f"LUMO+{diff - 1}"
                elif diff < 0 and diff > -100:
                    homo_lumo_label = f"HOMO{diff}"

            occ = mo.get("occ", mo.get("occupation", 0.0))
            e_eh = mo.get("energy_eh", mo.get("energy", None))
            e_ev = mo.get("energy_ev", None)

            if e_eh is None and e_ev is not None:
                e_eh = e_ev / 27.2114
            elif e_ev is None and e_eh is not None:
                e_ev = e_eh * 27.2114
            if e_eh is None:
                e_eh = 0.0
            if e_ev is None:
                e_ev = 0.0

            # Check if file exists to highlight
            bg_color = None
            try:
                # Use label_id as the display_id to match logic in visualize_current_mo
                path = self.get_cube_path(label_id)
                if path and os.path.exists(path):
                    bg_color = QColor(240, 255, 240)  # Light Green
            except Exception as _e:
                logging.warning("[mo_analysis.py:388] silenced: %s", _e)

            item = QTreeWidgetItem(
                [label_id, homo_lumo_label, f"{occ:.2f}", f"{e_ev:.3f}", f"{e_eh:.5f}"]
            )
            if bg_color:
                for c in range(5):
                    item.setBackground(c, QBrush(bg_color))

            # Store the unique lookup key
            item.setData(0, Qt.ItemDataRole.UserRole, access_key)

            # Text Color based on Occ
            if occ > 1.9:
                item.setForeground(0, QColor("blue"))
            elif occ > 0.1:
                item.setForeground(0, QColor("green"))
            else:
                item.setForeground(0, QColor("gray"))

            self.tree.addTopLevelItem(item)

        # Scroll to HOMO of first available spin branch
        # Since we just fill them, scrolling to center might be ambiguous.
        # Let's scroll to the first 'Occupied' item found from top of tree (High energy)
        # That would be the HOMO of the last spin block added (e.g. Beta HOMO if present, or Alpha HOMO).
        # Actually tree has High Energy at top. So First Occupied Item is HOMO.

        iterator = QTreeWidgetItemIterator(self.tree)
        while iterator.value():
            item = iterator.value()
            # Check column 1 for Label "HOMO"
            if "HOMO" in item.text(1):
                self.tree.setCurrentItem(item)
                self.tree.scrollToItem(
                    item, QAbstractItemView.ScrollHint.PositionAtCenter
                )
                break
            iterator += 1

    def on_double_click(self, item, col):
        # Double click always tries to visualize (generate if needed)
        self.visualize_selected_mos()

    def on_item_changed(self, current, previous):
        # Single click or keyboard change
        if not current:
            return

        # 1. Update Buttons (reuse logic)
        # on_selection_changed handles button state but relies on "selectedItems"
        # which might not be updated yet? usually it is.
        # But let's check for Auto-Load

        try:
            # Logic: If generated file exists, load it.
            display_id = current.text(0)
            # Need to clean display_id to match file system (e.g. "17 (a)" -> "17_a")
            # Use self.get_cube_path which handles this
            path = self.get_cube_path(display_id)
            if path and os.path.exists(path):
                self.show_cube(path)
        except Exception as _e:
            logging.warning("[mo_analysis.py:445] silenced: %s", _e)

    def on_selection_changed(self):
        items = self.tree.selectedItems()
        has_coeffs = False
        if items:
            try:
                # Use UserRole data for robust lookup
                key = items[0].data(0, Qt.ItemDataRole.UserRole)

                # Check coeffs in parser data
                if hasattr(self.parent_dlg, "parser") and self.parent_dlg.parser:
                    if (
                        hasattr(self.parent_dlg.parser, "data")
                        and self.parent_dlg.parser.data
                    ):
                        if "mo_coeffs" in self.parent_dlg.parser.data:
                            # Use key directly
                            if key in self.parent_dlg.parser.data["mo_coeffs"]:
                                has_coeffs = True
            except Exception as _e:
                logging.warning("[mo_analysis.py:462] silenced: %s", _e)

        self.btn_vis.setEnabled(has_coeffs)
        if items and not has_coeffs:
            self.btn_vis.setToolTip(
                f"No coefficients available for MO {items[0].text(0)}"
            )
            self.lbl_warning.setText(
                "<b>Coefficients Missing.</b> Required Input:"
                "<pre style='margin-top:0; margin-bottom:0;'>%output<br>  Print[P_Basis] 2<br>  Print[P_Mos] 1<br>end</pre>"
            )
            self.lbl_warning.setVisible(True)
            self.btn_copy_input.setVisible(True)
        else:
            self.btn_vis.setToolTip("")
            self.lbl_warning.setVisible(False)
            self.btn_copy_input.setVisible(False)

    def copy_orca_input(self):
        text = "%output\n  Print[P_Basis] 2\n  Print[P_Mos] 1\nend"
        QApplication.clipboard().setText(text)
        if self.mw and hasattr(self.mw, "statusBar"):
            sb = self.mw.statusBar()
            if sb:
                sb.showMessage("ORCA Input block copied to clipboard.", 5000)
            else:
                print("ORCA Input block copied to clipboard.")
        else:
            print("ORCA Input block copied to clipboard.")

    def get_engine(self):
        if not BasisSetEngine:
            QMessageBox.critical(self, "Error", "BasisSetEngine not available")
            return None
        try:
            # Safety checks
            if not hasattr(self.parent_dlg, "parser") or not self.parent_dlg.parser:
                return None
            if (
                not hasattr(self.parent_dlg.parser, "data")
                or not self.parent_dlg.parser.data
            ):
                return None

            shells = self.parent_dlg.parser.data.get("basis_set_shells", [])

            # Convert lists to numpy arrays for the engine AND map keys
            # IMPORTANT: Engine expects Bohr for centers (matching internal grid)
            # Parser provides Angstroms.
            BOHR_TO_ANG = 0.529177249

            clean_shells = []
            for s in shells:
                center_ang = np.array(s.get("origin", s.get("center", [0, 0, 0])))
                # Convert Angstrom -> Bohr
                center_bohr = center_ang / BOHR_TO_ANG

                d = {
                    "type": s.get("l", s.get("type", 0)),
                    "center": center_bohr,
                    "exps": np.array(s["exps"]),
                    "coeffs": np.array(s["coeffs"]),
                }
                clean_shells.append(d)

            return BasisSetEngine(clean_shells)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Engine Init Failed: {e}")
            return None

    def visualize_selected_mos(self):
        # Batch generation for selected items
        selected = self.tree.selectedItems()
        if not selected:
            return

        # Build Queue
        self.generation_queue = []
        for item in selected:
            key = item.data(0, Qt.ItemDataRole.UserRole)
            if key is not None:  # Ensure key is valid
                self.generation_queue.append(key)

        if not self.generation_queue:
            return

        # Start Batch
        self.process_generation_queue()

    def process_generation_queue(self):
        if not getattr(self, "generation_queue", None):
            # Done
            if (
                getattr(self, "progress_dialog", None) is not None
                and self.progress_dialog
            ):
                self.progress_dialog.close()
                self.progress_dialog = None
            return

        # Get next key
        key = self.generation_queue.pop(0)

        # Check if worker running
        # Actually we force sequential.

        self._generate_single_mo(key)

    def _generate_single_mo(self, key):
        # ... (Logic from old visualize_current_mo) ...
        # Re-fetch data using key
        # Key might be index (int) or dict key (str)
        # Find the correct MO in self.mo_list
        mo_data = None
        for mo in self.mo_list:
            if mo.get("_access_key", None) == key:
                mo_data = mo
                break

        if not mo_data:
            print(f"MO not found for key: {key}")
            self.process_generation_queue()  # Skip invalid key
            return

        # Get Display ID
        # Need to get the display_id from the original item for consistent caching/display
        # Find the item in the tree by key
        display_id = None
        it = QTreeWidgetItemIterator(self.tree)
        while it.value():
            item = it.value()
            if item.data(0, Qt.ItemDataRole.UserRole) == key:
                display_id = item.text(0)
                break
            it += 1

        if display_id is None:
            # Fallback if item not found (shouldn't happen if key is from selected item)
            display_id = str(mo_data.get("id", key))
            spin = mo_data.get("spin", "")
            if spin == "alpha":
                display_id += " (a)"
            elif spin == "beta":
                display_id += " (b)"

        # Safety checks for parser data
        if not hasattr(self.parent_dlg, "parser") or not self.parent_dlg.parser:
            QMessageBox.warning(self, "Error", "Parser not available")
            self.process_generation_queue()
            return
        if (
            not hasattr(self.parent_dlg.parser, "data")
            or not self.parent_dlg.parser.data
        ):
            QMessageBox.warning(self, "Error", "No parser data available")
            self.process_generation_queue()
            return

        coeffs_map = self.parent_dlg.parser.data.get("mo_coeffs", {})
        mo_data_coeffs = coeffs_map.get(key, None)  # This is the actual coeffs data
        if not mo_data_coeffs:
            QMessageBox.warning(self, "Error", f"No coefficients for MO {display_id}")
            self.process_generation_queue()
            return

        engine = self.get_engine()
        if not engine:
            self.process_generation_queue()  # Skip
            return

        if engine.n_basis == 0:
            QMessageBox.warning(
                self,
                "Error",
                "Basis Set Engine initialized with 0 basis functions. Check parser.",
            )
            self.process_generation_queue()
            return

        # Prepare Data
        raw_coeffs = [x["coeff"] for x in mo_data_coeffs["coeffs"]]
        dense_vec = np.zeros(engine.n_basis)
        n = min(len(raw_coeffs), engine.n_basis)
        dense_vec[:n] = raw_coeffs[:n]

        if np.sum(np.abs(dense_vec)) < 1e-9:
            QMessageBox.warning(
                self,
                "Data Error",
                f"MO coefficients for {display_id} are all zero/empty! Cube will be empty.",
            )
            self.process_generation_queue()
            return

        out_path = self.get_cube_path(display_id)
        if not out_path:
            out_path = os.path.join(tempfile.gettempdir(), f"orca_mo_{display_id}.cube")
        else:
            # Lazy directory creation: create it now if it doesn't exist
            out_dir = os.path.dirname(out_path)
            if not os.path.exists(out_dir):
                try:
                    os.makedirs(out_dir)
                except Exception as _e:
                    logging.warning("[mo_analysis.py:644] silenced: %s", _e)

        self.last_cube_path = out_path

        if os.path.exists(out_path):
            self.show_cube(out_path)
            # Highlight
            it = QTreeWidgetItemIterator(self.tree)
            bg = QBrush(QColor(240, 255, 240))
            while it.value():
                item = it.value()
                if item.data(0, Qt.ItemDataRole.UserRole) == key:
                    for c in range(5):
                        item.setBackground(c, bg)
                    break
                it += 1

            # Next!
            self.process_generation_queue()
            return

        # Start Worker
        # Manage Progress Dialog (Shared)
        if getattr(self, "progress_dialog", None) is None:
            self.progress_dialog = QProgressDialog(
                "Generating Cubes...", "Cancel", 0, 100, self
            )
            self.progress_dialog.setWindowModality(Qt.WindowModality.WindowModal)
            self.progress_dialog.setAutoClose(False)  # Keep open for batch
            self.progress_dialog.show()

        self.progress_dialog.setLabelText(f"Generating {display_id}...")
        self.progress_dialog.setValue(0)

        self.worker = CalcWorker(
            engine,
            key,
            self.spin_pts.value(),
            self.spin_margin.value(),
            self.parent_dlg.parser.data["atoms"],
            self.parent_dlg.parser.data["coords"],
            dense_vec,
            out_path,
        )

        self.worker.progress_sig.connect(self.progress_dialog.setValue)

        def on_finished(success, res):
            if success:
                self.show_cube(res)
                # Highlight
                it = QTreeWidgetItemIterator(self.tree)
                bg = QBrush(QColor(240, 255, 240))
                while it.value():
                    item = it.value()
                    if item.data(0, Qt.ItemDataRole.UserRole) == key:
                        for c in range(5):
                            item.setBackground(c, bg)
                        break
                    it += 1

                if (
                    getattr(self, "diag_dlg", None) is not None
                    and self.diag_dlg
                    and self.diag_dlg.isVisible()
                ):
                    self.diag_dlg.status_label.setText(
                        f"Generated: {os.path.basename(res)}"
                    )
            else:
                # If one fails, maybe continue?
                # Or stop? let's continue but warn?
                print(f"Failed: {res}")
                QMessageBox.warning(
                    self, "Generation Failed", f"Failed to generate cube:\n{res}"
                )
                if (
                    getattr(self, "diag_dlg", None) is not None
                    and self.diag_dlg
                    and self.diag_dlg.isVisible()
                ):
                    self.diag_dlg.status_label.setText("Generation Failed")

            # Trigger Next
            self.process_generation_queue()

        self.worker.finished_sig.connect(on_finished)
        self.worker.start()

    def pick_color(self, which):
        current_col = QColor("red") if which == "p" else QColor("blue")
        # Try to parse from button style
        try:
            style = (
                self.btn_color_p.styleSheet()
                if which == "p"
                else self.btn_color_n.styleSheet()
            )
            if "background-color:" in style:
                c_str = style.split("background-color:")[1].split(";")[0].strip()
                current_col = QColor(c_str)
        except Exception as _e:
            logging.warning("[mo_analysis.py:719] silenced: %s", _e)

        col = QColorDialog.getColor(current_col, self, "Select Color")
        if col.isValid():
            hex_c = col.name()
            # Determine contrasting text color
            # Simple brightness check
            brightness = (col.red() * 299 + col.green() * 587 + col.blue() * 114) / 1000
            text_c = "black" if brightness > 128 else "white"

            style_sheet = (
                f"background-color: {hex_c}; color: {text_c}; font-weight: bold;"
            )
            if which == "p":
                self.btn_color_p.setStyleSheet(style_sheet)
            else:
                self.btn_color_n.setStyleSheet(style_sheet)
            self.update_vis_only()

    def load_settings(self):
        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")
        self.presets = {
            "Default": {
                "iso": 0.02,
                "opacity": 0.5,
                "style": "Surface",
                "color_p": "#ff0000",
                "color_n": "#0000ff",
                "smooth_shading": True,
            }
        }

        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r") as f:
                    all_settings = json.load(f)
                    mo_settings = all_settings.get("mo_settings", {})

                    # Load presets
                    saved_presets = mo_settings.get("presets", {})
                    for name, data in saved_presets.items():
                        self.presets[name] = data

                    # Last used
                    last_preset = mo_settings.get("last_preset", "Default")

                    # Populate combo
                    self.combo_presets.blockSignals(True)
                    self.combo_presets.clear()
                    self.combo_presets.addItems(list(self.presets.keys()))

                    if last_preset in self.presets:
                        self.combo_presets.setCurrentText(last_preset)
                        self.apply_preset(last_preset)
                    else:
                        self.combo_presets.setCurrentText("Default")
                        self.apply_preset("Default")

                    self.combo_presets.blockSignals(False)
            except Exception as e:
                print(f"Error loading settings: {e}")

    def save_settings(self):
        # Save current presets and selection
        all_settings = {}
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r") as f:
                    all_settings = json.load(f)
            except Exception as _e:
                logging.warning("[mo_analysis.py:779] silenced: %s", _e)

        mo_settings = {
            "presets": {k: v for k, v in self.presets.items() if k != "Default"},
            "last_preset": self.combo_presets.currentText(),
            "smooth_shading": self.check_smooth.isChecked(),
        }
        all_settings["mo_settings"] = mo_settings

        try:
            with open(self.settings_file, "w") as f:
                json.dump(all_settings, f, indent=2)
        except Exception as e:
            print(f"Error saving settings: {e}")

    def save_preset(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:")
        if not ok or not name:
            return

        # Get current state
        data = {
            "iso": self.spin_iso.value(),
            "opacity": self.spin_opacity.value(),
            "style": self.combo_style.currentText(),
            "color_p": self.get_color_hex("p"),
            "color_n": self.get_color_hex("n"),
            "smooth_shading": self.check_smooth.isChecked(),
        }

        self.presets[name] = data

        # Update combo
        self.combo_presets.currentText()
        self.combo_presets.blockSignals(True)
        self.combo_presets.clear()
        self.combo_presets.addItems(list(self.presets.keys()))
        self.combo_presets.setCurrentText(name)  # switch to new
        self.combo_presets.blockSignals(False)

        self.save_settings()

    def delete_preset(self):
        curr = self.combo_presets.currentText()
        if curr == "Default":
            QMessageBox.warning(self, "Error", "Cannot delete Default preset.")
            return

        if curr in self.presets:
            del self.presets[curr]

        self.combo_presets.blockSignals(True)
        self.combo_presets.clear()
        self.combo_presets.addItems(list(self.presets.keys()))
        self.combo_presets.setCurrentText("Default")
        self.apply_preset("Default")
        self.combo_presets.blockSignals(False)

        self.save_settings()

    def apply_preset(self, name):
        if name not in self.presets:
            return
        data = self.presets[name]

        # Block signals to avoid partial updates? No, updates are fine.
        self.spin_iso.setValue(data.get("iso", 0.02))
        self.spin_opacity.setValue(data.get("opacity", 0.5))

        style = data.get("style", "Surface")
        idx = self.combo_style.findText(style)
        if idx >= 0:
            self.combo_style.setCurrentIndex(idx)

        cp = data.get("color_p", "#ff0000")
        cn = data.get("color_n", "#0000ff")
        self.set_btn_color(self.btn_color_p, cp)
        self.set_btn_color(self.btn_color_n, cn)

        self.check_smooth.setChecked(data.get("smooth_shading", True))

        self.update_vis_only()
        self.save_settings()  # Save last used

    def get_color_hex(self, which):
        # Extract from stylesheet
        btn = self.btn_color_p if which == "p" else self.btn_color_n
        style = btn.styleSheet()
        if "background-color:" in style:
            return style.split("background-color:")[1].split(";")[0].strip()
        return "#ff0000" if which == "p" else "#0000ff"

    def set_btn_color(self, btn, hex_c):
        col = QColor(hex_c)
        brightness = (col.red() * 299 + col.green() * 587 + col.blue() * 114) / 1000
        text_c = "black" if brightness > 128 else "white"
        btn.setStyleSheet(
            f"background-color: {hex_c}; color: {text_c}; font-weight: bold;"
        )

    def update_vis_only(self):
        if self.last_cube_path and os.path.exists(self.last_cube_path):
            self.show_cube(self.last_cube_path)

    def show_cube(self, path):
        if not CubeVisualizer:
            print("Warning: CubeVisualizer module not loaded.")
            QMessageBox.warning(
                self,
                "Visualizer Error",
                "CubeVisualizer module not loaded.\nCheck if 'pyvista' is installed.",
            )
            return

        # Access Main Window
        mw = self.mw
        if not mw:
            return

        cp = self.get_color_hex("p")
        cn = self.get_color_hex("n")
        style = self.combo_style.currentText().lower()  # surface, wireframe, points

        vis = CubeVisualizer(mw)
        if vis.load_file(path):
            vis.show_iso(
                self.spin_iso.value(),
                opacity=self.spin_opacity.value(),
                color_p=cp,
                color_n=cn,
                style=style,
                smooth_shading=self.check_smooth.isChecked(),
            )
            mw.plotter.render()

    def closeEvent(self, event):
        """Clean up 3D actors when closing"""
        # Clean up tracked sub-dialogs
        if getattr(self, "energy_dlg", None) is not None and self.energy_dlg:
            try:
                self.energy_dlg.close()
            except Exception as _e:
                logging.warning("[mo_analysis.py:905] silenced: %s", _e)
            self.energy_dlg = None

        if hasattr(self.parent_dlg, "mw"):
            plotter = self.parent_dlg.mw.plotter
            plotter.remove_actor("mo_iso_p")
            plotter.remove_actor("mo_iso_n")
            plotter.render()
        elif hasattr(self.parent_dlg, "context"):
            plotter = self.parent_dlg.context.get_main_window().plotter
            plotter.remove_actor("mo_iso_p")
            plotter.remove_actor("mo_iso_n")
            plotter.render()

        super().closeEvent(event)

    def export_csv(self):
        default_path = get_default_export_path(
            self.parent_dlg.file_path, suffix="_mo_list", extension=".csv"
        )
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export MO Data", default_path, "CSV Files (*.csv)"
        )
        if not filename:
            return

        try:
            import csv

            with open(filename, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["ID", "Occupation", "Energy (eV)", "Energy (Eh)"])

                # Iterate tree items to respect current sort/display?
                # Using the data list is safer and complete.
                # However, the tree has formatted text. Let's use the list self.mo_list
                # But self.mo_list is raw. Let's replicate the display logic or use tree items.

                # Using tree items ensures we export what is seen (including calculated Energy if missing)
                it = QTreeWidgetItemIterator(self.tree)
                while it.value():
                    item = it.value()
                    # ID, Occ, eV, Eh
                    row = [item.text(0), item.text(1), item.text(2), item.text(3)]
                    writer.writerow(row)
                    it += 1
            # print(f"Data exported to {filename}")
            # QMessageBox.information(self, "Success", f"Data exported to {filename}")
            if self.mw and hasattr(self.mw, "statusBar"):
                sb = self.mw.statusBar()
                if sb:
                    sb.showMessage(f"Data exported to {filename}", 5000)
                else:
                    print(f"Data exported to {filename}")
            else:
                print(f"Data exported to {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to export CSV: {e}")

    def show_mo_diagram(self):
        if not EnergyDiagramDialog:
            return

        # Prepare Data for Diagram
        # 1. Separate by spin
        mos_alpha = []
        mos_beta = []
        mos_restr = []

        for mo in self.mo_list:
            s = mo.get("spin", "restricted")
            e_eh = mo.get("energy_eh", mo.get("energy", None))
            e_ev = mo.get("energy_ev", None)
            # Ensure e_eh is set (EnergyDiagram usually expects Hartree if unit selection is implemented there, defaults to "Ha")
            if e_eh is None and e_ev is not None:
                e_eh = e_ev / 27.2114
            if e_eh is None:
                e_eh = 0.0

            occ = mo.get("occ", mo.get("occupation", 0.0))

            item = {"e": e_eh, "occ": occ}

            if s == "alpha":
                mos_alpha.append(item)
            elif s == "beta":
                mos_beta.append(item)
            else:
                mos_restr.append(item)

        # 2. Determine Type
        is_uhf = len(mos_beta) > 0

        diag_data = {}
        if is_uhf:
            diag_data["type"] = "UHF"
            # Lists must be sorted by energy usually? They are sorted in list by spin/index.
            # Assuming 'mo_list' was sorted by index effectively
            # Let's ensure sort by index or energy
            # Since self.mo_list keys were keys...
            # Let's rely on standard order.

            es_a = [m["e"] for m in mos_alpha]
            noccs_a = [m["occ"] for m in mos_alpha]
            es_b = [m["e"] for m in mos_beta]
            noccs_b = [m["occ"] for m in mos_beta]

            diag_data["energies"] = [es_a, es_b]
            diag_data["occupations"] = [noccs_a, noccs_b]
        else:
            diag_data["type"] = "RHF"
            es = [m["e"] for m in mos_restr]
            noccs = [m["occ"] for m in mos_restr]
            if not es and mos_alpha:  # Fallback if labeled alpha but no beta
                es = [m["e"] for m in mos_alpha]
                noccs = [m["occ"] for m in mos_alpha]

            diag_data["energies"] = es
            diag_data["occupations"] = noccs

        # Determine Result Dir
        res_dir = None
        if hasattr(self.parent_dlg, "parser") and self.parent_dlg.parser:
            if (
                hasattr(self.parent_dlg.parser, "filename")
                and self.parent_dlg.parser.filename
            ):
                res_dir = os.path.dirname(self.parent_dlg.parser.filename)
                # Check for dedicated cube dir?
                base = os.path.splitext(
                    os.path.basename(self.parent_dlg.parser.filename)
                )[0]
                cube_dir = os.path.join(res_dir, f"{base}_cubes")
                if os.path.exists(cube_dir):
                    res_dir = cube_dir

        # Make modeless
        if getattr(self, "energy_dlg", None) is not None and self.energy_dlg:
            self.energy_dlg.close()

        self.energy_dlg = EnergyDiagramDialog(
            diag_data, parent=self, result_dir=res_dir
        )
        self.energy_dlg.show()

    def load_file_by_path(self, path):
        """Called from Diagram to load valid existing file"""
        if os.path.exists(path):
            self.show_cube(path)

            # Highlight in tree if possible
            # Need to reverse map path to key? Or user can just view it.
            # Best effort to highlight
            # We don't have the key here easily unless we parse filename.
            # But the visualization is what matters.

    def generate_specific_orbital(self, index, label, spin_suffix=""):
        """Called from Diagram to generate cube"""
        # We need to map index -> Key.
        # self.mo_list has keys.
        # Index is 0-based.
        # But separate by spin?
        # Diagram index is index within spin channel.
        # We need to find the MO with that index and spin.

        found_mo = None
        target_spin = "beta" if spin_suffix == "_B" else "alpha"  # or restricted
        # Diagram treats 'restricted' as alpha usually (col 0).

        # Refined Logic:
        # iterate self.mo_list, count indices for that spin
        # 0-based index means nth orbital of that spin.

        # Check type
        is_uhf = False
        for mo in self.mo_list:
            if mo.get("spin", None) == "beta":
                is_uhf = True
                break

        if not is_uhf:
            target_spin = "restricted"

        # If we asked for Alpha but it's restricted, map alpha->restricted
        if target_spin == "alpha" and not is_uhf:
            target_spin = "restricted"

        curr_idx = 0
        mo_key = None

        # Sort self.mo_list by ID to match diagram order?
        # self.mo_list is appended in loop.
        # normalize_and_populate sorts keys.
        # We need to trust the sort order is: Alpha 0..N, Beta 0..N ??
        # Or keys are arbitrary.
        # Let's re-sort to be safe.

        sorted_mos = sorted(
            self.mo_list, key=lambda x: (x.get("spin", ""), int(x.get("id", -1)))
        )

        for mo in sorted_mos:
            s = mo.get("spin", "restricted")
            if s != target_spin:
                continue

            if curr_idx == index:
                mo_key = mo.get("_access_key", None)
                found_mo = mo
                break
            curr_idx += 1

        if found_mo:
            # Trigger Visualization
            # We need to select it in the tree?
            # Or just call worker directly.
            # Visualizing relies on selection in tree for `visualize_current_mo`.
            # Let's Select it in tree.

            # Find item with data key
            it = QTreeWidgetItemIterator(self.tree)
            while it.value():
                item = it.value()
                if item.data(0, Qt.ItemDataRole.UserRole) == mo_key:
                    self.tree.setCurrentItem(item)
                    self.visualize_selected_mos()
                    break
                it += 1
        else:
            QMessageBox.warning(
                self, "Error", f"Could not find MO for Index {index} ({target_spin})"
            )

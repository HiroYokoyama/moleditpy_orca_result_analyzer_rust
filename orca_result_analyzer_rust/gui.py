import os
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QWidget,
    QGridLayout,
    QMessageBox,
    QMenuBar,
    QFileDialog,
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import QSize

try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from rdkit.Chem import rdDetermineBonds
except ImportError:
    Chem = None
    Point3D = None
    rdDetermineBonds = None

# Imported Modules for Analysis
from .mo_analysis import MODialog
from .freq_analysis import FrequencyDialog
from .traj_analysis import TrajectoryResultDialog
from .force_analysis import ForceViewerDialog
from .charge_analysis import ChargeDialog
from .dipole_analysis import DipoleDialog
from .nmr_analysis import NMRDialog
from .tddft_analysis import TDDFTDialog
from .thermal_analysis import ThermalTableDialog
from .scf_analysis import SCFTraceDialog

from . import PLUGIN_VERSION
import logging
# from .logger import Logger


class OrcaResultAnalyzerDialog(QDialog):
    def __init__(self, parent, parser, file_path, context=None):
        super().__init__(parent)
        self.mw = parent
        self.parser = parser
        self.file_path = file_path
        self.context = context

        self.setWindowTitle(f"ORCA Result Analyzer (v{PLUGIN_VERSION})")
        self.resize(450, 600)

        # self.logger = Logger.get_logger("OrcaResultAnalyzerDialog")
        self.init_ui()

        # Update main window title to reflect ORCA result
        if hasattr(self.mw, "init_manager"):
            self.mw.init_manager.current_file_path = self.file_path
        if hasattr(self.mw, "state_manager") and hasattr(
            self.mw.state_manager, "update_window_title"
        ):
            self.mw.state_manager.update_window_title()

    def get_icon(self, name):
        """Helper to load icon from icon directory"""
        icon_path = os.path.join(os.path.dirname(__file__), "icon", name)
        if os.path.exists(icon_path):
            return QIcon(icon_path)
        return QIcon()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # Menu Bar (added as widget since QDialog doesn't have native menu bar)
        menu_bar = QMenuBar(self)
        menu_bar.setStyleSheet(
            "QMenuBar { background-color: #f0f0f0; padding: 5px; font-size: 11pt; }"
        )
        layout.addWidget(menu_bar)

        # File Menu
        file_menu = menu_bar.addMenu("&File")

        open_action = QAction("&Open File...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.setIcon(self.get_icon("menu_open.svg"))
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        reload_action = QAction("&Reload File", self)
        reload_action.setShortcut("Ctrl+R")
        reload_action.setIcon(self.get_icon("menu_reload.svg"))
        reload_action.triggered.connect(self.reload_file)
        file_menu.addAction(reload_action)

        file_menu.addSeparator()

        close_action = QAction("&Close", self)
        close_action.setShortcut("Ctrl+W")
        close_action.setIcon(self.get_icon("menu_close.svg"))
        close_action.triggered.connect(self.close)
        file_menu.addAction(close_action)

        # Current File Display with Open Button
        file_frame = QWidget()
        file_frame.setStyleSheet("""
            QWidget {
                background-color: #e8f4f8;
                border: 2px solid #0066cc;
                border-radius: 5px;
                padding: 10px;
            }
        """)
        file_frame_layout = QHBoxLayout(file_frame)
        file_frame_layout.setContentsMargins(10, 10, 10, 10)

        # File info section
        file_info_layout = QVBoxLayout()
        lbl_current = QLabel("<b>Current File:</b>")
        lbl_current.setStyleSheet(
            "font-size: 10pt; background: transparent; border: none; padding: 0;"
        )
        file_info_layout.addWidget(lbl_current)

        self.lbl_file_path = QLabel(self.file_path)
        self.lbl_file_path.setStyleSheet(
            "color: #0066cc; font-size: 9pt; background: transparent; border: none; padding: 0;"
        )
        self.lbl_file_path.setWordWrap(True)
        file_info_layout.addWidget(self.lbl_file_path)

        # Updated Time Label
        self.lbl_updated = QLabel("Updated: ---")
        self.lbl_updated.setStyleSheet(
            "color: #555; font-size: 9pt; background: transparent; border: none; padding: 0;"
        )
        file_info_layout.addWidget(self.lbl_updated)

        # ORCA Version Label
        version_str = (
            self.parser.data.get("version", "Unknown") if self.parser else "Unknown"
        )
        self.lbl_version = QLabel(f"ORCA Version: {version_str}")
        self.lbl_version.setStyleSheet(
            "color: #555; font-size: 9pt; background: transparent; border: none; padding: 0;"
        )
        file_info_layout.addWidget(self.lbl_version)

        # Initial call to update labels
        self.update_file_info_labels()

        file_frame_layout.addLayout(file_info_layout, 1)

        # Buttons Layout (Open + Reload)
        btns_top_layout = QVBoxLayout()

        # Large Open File Button
        btn_open_large = QPushButton("Open File...")
        btn_open_large.setIcon(self.get_icon("menu_open.svg"))
        btn_open_large.setStyleSheet("""
            QPushButton {
                background-color: #0066cc;
                color: white;
                font-size: 10pt;
                font-weight: bold;
                padding: 8px 15px;
                border-radius: 5px;
                border: none;
                text-align: left;
            }
            QPushButton:hover {
                background-color: #0052a3;
            }
            QPushButton:pressed {
                background-color: #003d7a;
            }
        """)
        btn_open_large.clicked.connect(self.open_file)
        btns_top_layout.addWidget(btn_open_large)

        # Reload Button
        btn_reload = QPushButton("Reload")
        btn_reload.setIcon(self.get_icon("menu_reload.svg"))
        btn_reload.setStyleSheet("""
            QPushButton {
                background-color: #28a745;
                color: white;
                font-size: 10pt;
                font-weight: bold;
                padding: 8px 15px;
                border-radius: 5px;
                border: none;
                text-align: left;
            }
            QPushButton:hover {
                background-color: #218838;
            }
            QPushButton:pressed {
                background-color: #1e7e34;
            }
        """)
        btn_reload.clicked.connect(self.reload_file)
        btns_top_layout.addWidget(btn_reload)

        file_frame_layout.addLayout(btns_top_layout)

        layout.addWidget(file_frame)

        # Grid for buttons
        grid = QGridLayout()
        grid.setSpacing(12)
        layout.addLayout(grid)

        # Button style for all analysis buttons
        button_style = """
            QPushButton {
                background-color: #ffffff;
                border: 2px solid #cccccc;
                border-radius: 8px;
                padding: 12px 15px;
                font-size: 10pt;
                font-weight: bold;
                text-align: left;
                min-height: 40px;
            }
            QPushButton:hover {
                background-color: #f0f8ff;
                border: 2px solid #0066cc;
            }
            QPushButton:pressed {
                background-color: #e0e8ff;
            }
            QPushButton:disabled {
                background-color: #f0f0f0;
                border: 2px solid #dddddd;
                color: #999999;
            }
        """

        icon_size = QSize(32, 32)

        # Row 0: Electronic Structure
        self.btn_scf = QPushButton("SCF Trace")
        self.btn_scf.setIcon(self.get_icon("icon_scf.svg"))
        self.btn_scf.setIconSize(icon_size)
        self.btn_scf.setStyleSheet(button_style)
        self.btn_scf.clicked.connect(self.show_scf_trace)
        grid.addWidget(self.btn_scf, 0, 0)

        self.btn_mo = QPushButton("MO Analysis")
        self.btn_mo.setIcon(self.get_icon("icon_mo.svg"))
        self.btn_mo.setIconSize(icon_size)
        self.btn_mo.setStyleSheet(button_style)
        self.btn_mo.clicked.connect(self.show_mo_analyzer)
        grid.addWidget(self.btn_mo, 0, 1)

        # Row 1: Geometry Trajectory
        self.btn_traj = QPushButton("Optimization / Scan")
        self.btn_traj.setIcon(self.get_icon("icon_traj.svg"))
        self.btn_traj.setIconSize(icon_size)
        self.btn_traj.setStyleSheet(button_style)
        self.btn_traj.clicked.connect(self.show_trajectory)
        grid.addWidget(self.btn_traj, 1, 0)

        self.btn_forces = QPushButton("Forces")
        self.btn_forces.setIcon(self.get_icon("icon_forces.svg"))
        self.btn_forces.setIconSize(icon_size)
        self.btn_forces.setStyleSheet(button_style)
        self.btn_forces.clicked.connect(self.show_forces)
        grid.addWidget(self.btn_forces, 1, 1)

        # Row 2: Atomic Properties
        self.btn_charge = QPushButton("Atomic Charges")
        self.btn_charge.setIcon(self.get_icon("icon_charge.svg"))
        self.btn_charge.setIconSize(icon_size)
        self.btn_charge.setStyleSheet(button_style)
        self.btn_charge.clicked.connect(self.show_charges)
        grid.addWidget(self.btn_charge, 2, 0)

        self.btn_dipole = QPushButton("Dipole Moment")
        self.btn_dipole.setIcon(self.get_icon("icon_dipole.svg"))
        self.btn_dipole.setIconSize(icon_size)
        self.btn_dipole.setStyleSheet(button_style)
        self.btn_dipole.clicked.connect(self.show_dipole)
        grid.addWidget(self.btn_dipole, 2, 1)

        # Row 3: Vibrations & Thermodynamics
        self.btn_freq = QPushButton("Frequencies")
        self.btn_freq.setIcon(self.get_icon("icon_freq.svg"))
        self.btn_freq.setIconSize(icon_size)
        self.btn_freq.setStyleSheet(button_style)
        self.btn_freq.clicked.connect(self.show_freq)
        grid.addWidget(self.btn_freq, 3, 0)

        self.btn_therm = QPushButton("Thermochemistry")
        self.btn_therm.setIcon(self.get_icon("icon_therm.svg"))
        self.btn_therm.setIconSize(icon_size)
        self.btn_therm.setStyleSheet(button_style)
        self.btn_therm.clicked.connect(self.show_thermal)
        grid.addWidget(self.btn_therm, 3, 1)

        # Row 4: Advanced Spectroscopy
        self.btn_tddft = QPushButton("TDDFT")
        self.btn_tddft.setIcon(self.get_icon("icon_tddft.svg"))
        self.btn_tddft.setIconSize(icon_size)
        self.btn_tddft.setStyleSheet(button_style)
        self.btn_tddft.clicked.connect(self.show_tddft)
        grid.addWidget(self.btn_tddft, 4, 0)

        self.btn_nmr = QPushButton("NMR")
        self.btn_nmr.setIcon(self.get_icon("icon_nmr.svg"))
        self.btn_nmr.setIconSize(icon_size)
        self.btn_nmr.setStyleSheet(button_style)
        self.btn_nmr.clicked.connect(self.show_nmr)
        grid.addWidget(self.btn_nmr, 4, 1)

        layout.addStretch()

        # Close button only (Open is now in menu)
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_layout.addWidget(btn_close)

        layout.addLayout(btn_layout)

        self.update_button_states()

    def close_all_sub_dialogs(self):
        """Close all tracked analysis dialogs to prevent orphaned windows."""
        dialog_attrs = [
            "mo_dlg",
            "freq_dlg",
            "traj_dlg",
            "forces_dlg",
            "thermal_dlg",
            "tddft_dlg",
            "dipole_dlg",
            "charges_dlg",
            "nmr_dlg",
            "scf_dlg",
        ]
        for attr in dialog_attrs:
            if getattr(self, attr, None) is not None:
                dlg = getattr(self, attr)
                if dlg is not None:
                    try:
                        dlg.close()
                    except Exception as _e:
                        logging.warning("[gui.py:319] silenced: %s", _e)
                setattr(self, attr, None)

    def closeEvent(self, event):
        """Ensure all sub-dialogs close when the main analyzer window is closed."""
        self.close_all_sub_dialogs()
        super().closeEvent(event)

    def open_file(self):
        # Get last directory from current file
        start_dir = os.path.dirname(self.file_path) if self.file_path else ""

        path, _ = QFileDialog.getOpenFileName(
            self, "Open ORCA Output", start_dir, "ORCA Output (*.out)"
        )
        if not path:
            return

        self.load_file(path)

    def load_file(self, path):
        """Load a file and update UI"""
        # Close existing dialogs to prevent confusion
        self.close_all_sub_dialogs()

        try:
            content = ""
            encodings = ["utf-8", "utf-16", "latin-1", "cp1252"]
            found = False
            for enc in encodings:
                try:
                    with open(path, "r", encoding=enc) as f:
                        content = f.read()
                    found = True
                    break
                except UnicodeError:
                    continue

            if not found:
                with open(path, "r", encoding="utf-8", errors="replace") as f:
                    content = f.read()

            from .parser import OrcaParser

            new_parser = OrcaParser()
            new_parser.load_from_memory(content, path)

            # --- Auto-load NEB Trajectory if present ---
            # Try to get explicit filename from parser first
            parsed_trj = new_parser.data.get("neb_trj_file", None)
            base_dir = os.path.dirname(path)

            potential_paths = []
            if parsed_trj:
                potential_paths.append(os.path.join(base_dir, parsed_trj))

            # Fallback: Standard naming
            base, ext = os.path.splitext(path)
            potential_paths.append(base + "_MEP_trj.xyz")

            trj_path = None
            for p in potential_paths:
                if os.path.exists(p):
                    trj_path = p
                    break

            if trj_path:
                try:
                    with open(
                        trj_path, "r", encoding="utf-8", errors="replace"
                    ) as ftrj:
                        trj_content = ftrj.read()
                    trj_steps = new_parser.parse_xyz_content(trj_content)
                    if trj_steps:
                        # Prioritize detailed TRJ structure over summary
                        new_parser.data["scan_steps"] = trj_steps

                        # Update main structure with the last frame of the trajectory
                        if trj_steps[-1].get("atoms", None) and trj_steps[-1].get(
                            "coords", None
                        ):
                            new_parser.data["atoms"] = trj_steps[-1]["atoms"]
                            new_parser.data["coords"] = trj_steps[-1]["coords"]

                        # print(f"Loaded NEB Trajectory: {len(trj_steps)} frames from {os.path.basename(trj_path)}")
                        self.mw.statusBar().showMessage(
                            f"Loaded NEB Trajectory from {os.path.basename(trj_path)}",
                            5000,
                        )
                except Exception as e:
                    # print(f"Failed to load associated TRJ: {e}")
                    logging.warning("[gui.py:404] silenced: %s", e)
            else:
                pass

            self.parser = new_parser
            self.file_path = path

            # --- Sync with Main Window Title ---
            if hasattr(self.mw, "init_manager"):
                self.mw.init_manager.current_file_path = path
            if hasattr(self.mw, "state_manager") and hasattr(
                self.mw.state_manager, "update_window_title"
            ):
                self.mw.state_manager.update_window_title()

            # Update File Info Labels
            self.update_file_info_labels()

            # Auto-load 3D structure
            self.load_structure_3d()
            self.update_button_states()

            # QMessageBox.information(self, "Loaded", f"Successfully loaded:\n{os.path.basename(path)}")
            # print(f"Successfully loaded: {os.path.basename(path)}")
            self.mw.statusBar().showMessage(
                f"Successfully loaded: {os.path.basename(path)}", 5000
            )

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file:\n{e}")

    def update_file_info_labels(self):
        """Update the file info labels (Path, Time, Versions)"""
        if getattr(self, "lbl_file_path", None) is None:
            return

        self.lbl_file_path.setText(self.file_path)

        # Updated Time
        mtime_str = "---"
        if self.file_path and os.path.exists(self.file_path):
            try:
                from datetime import datetime

                dt = datetime.fromtimestamp(os.path.getmtime(self.file_path))
                mtime_str = dt.strftime("%Y-%m-%d %H:%M:%S")
            except Exception as _e:
                logging.warning("[gui.py:446] silenced: %s", _e)

        if getattr(self, "lbl_updated", None) is not None:
            self.lbl_updated.setText(f"Updated: {mtime_str}")

        # ORCA Version
        v = self.parser.data.get("version", "Unknown") if self.parser else "Unknown"
        if getattr(self, "lbl_version", None) is not None:
            self.lbl_version.setText(f"ORCA Version: {v}")

    def reload_file(self):
        if self.file_path and os.path.exists(self.file_path):
            self.load_file(self.file_path)
        else:
            QMessageBox.warning(
                self, "Error", "No file currently loaded or file not found."
            )

    def update_button_states(self):
        data = self.parser.data

        # Enable MO button if MO coefficients OR orbital energies exist
        # Strict check: must be non-empty dict or list
        mo_coeffs = data.get("mo_coeffs", {})
        orb_energies = data.get("orbital_energies", [])

        has_mo = bool(mo_coeffs) or (bool(orb_energies) and len(orb_energies) > 0)
        self.btn_mo.setEnabled(has_mo)
        tooltip = ""
        if not has_mo:
            tooltip = "No MO data found"
        elif orb_energies and not mo_coeffs:
            tooltip = "Orbital energies available (no coefficients for visualization)"
        self.btn_mo.setToolTip(tooltip)

        freqs = data.get("frequencies", [])
        has_freq = bool(freqs and len(freqs) > 0)
        self.btn_freq.setEnabled(has_freq)
        self.btn_freq.setToolTip("" if has_freq else "No frequency data found")

        # Scan / Optimization
        scan_steps = data.get("scan_steps", [])
        has_scan = bool(scan_steps and len(scan_steps) > 0)
        self.btn_traj.setEnabled(has_scan)
        self.btn_traj.setToolTip("" if has_scan else "No trajectory/scan steps found")

        # Enable Forces button if gradients exist or scan exists (BUT NOT FOR NEB SUMMARY)
        grads = data.get("gradients", [])
        scan_steps = data.get("scan_steps", [])

        # Check if scan steps are purely NEB images (summary only, no force data)
        # Note: If we auto-loaded XYZ, we now HAVE structure steps (type='neb_step' from xyz parser or similar)
        # xyz parser sets type='neb_step', summary parser sets type='neb_image'
        # If it's 'neb_step', it might have structure but NO FORCES.
        # Actually standard XYZ doesn't have forces.
        # So we should likely still keep Forces disabled if 'neb_step' unless we parsed gradients separately.

        is_neb_summary = False
        if scan_steps and isinstance(scan_steps, list):
            first_type = scan_steps[0].get("type", None)
            if first_type in [
                "neb_image",
                "neb_step",
            ]:  # Both summary and XYZ steps usually lack Forces in standard ORCA XYZ
                is_neb_summary = True

        has_forces = bool(grads)
        if scan_steps and not is_neb_summary:
            has_forces = True

        self.btn_forces.setEnabled(has_forces)
        tooltip = ""
        if not has_forces:
            tooltip = "No gradients or optimization trajectory found"
        elif grads:
            tooltip = "View Forces (Gradients)"
        elif scan_steps:
            tooltip = "View trajectory steps"
        self.btn_forces.setToolTip(tooltip)

        bool(
            data.get("thermal")
            or (data.get("frequencies", None) and "thermo" in str(data))
        )
        self.btn_therm.setEnabled(bool(data.get("thermal")))

        has_tddft = bool(data.get("tddft"))
        self.btn_tddft.setEnabled(has_tddft)

        has_dipole = bool(data.get("dipoles"))
        self.btn_dipole.setEnabled(has_dipole)

        has_charges = bool(data.get("charges"))
        self.btn_charge.setEnabled(has_charges)

        has_nmr = bool(data.get("nmr_shielding"))
        self.btn_nmr.setEnabled(has_nmr)

        # SCF Trace
        has_scf = bool(data.get("scf_traces"))
        self.btn_scf.setEnabled(has_scf)
        self.btn_scf.setToolTip("" if has_scf else "No SCF iteration data found")

    def load_structure_3d(self):
        # Helper to ensure 3D structure is loaded
        atoms = self.parser.data.get("atoms", [])
        coords = self.parser.data.get("coords", [])

        if not atoms or not coords:
            return

        if not Chem:
            return

        try:
            mol = Chem.RWMol()
            conf = Chem.Conformer(len(atoms))

            for i, sym in enumerate(atoms):
                idx = mol.AddAtom(Chem.Atom(sym))
                x, y, z = coords[i]
                conf.SetAtomPosition(idx, Point3D(x, y, z))

            mol.AddConformer(conf)

            # Bond determination must run on the RWMol (mutable) before GetMol().
            # DetermineBonds on a read-only Mol silently fails.
            # Skip if any animation is currently playing (traj or freq) — the result
            # would be immediately overwritten and the cost is wasted.
            _traj_playing = (
                getattr(self, "traj_dlg", None) is not None
                and getattr(self.traj_dlg, "is_playing", False)
            )
            _freq_playing = (
                getattr(self, "freq_dlg", None) is not None
                and getattr(self.freq_dlg, "is_playing", False)
            )
            if rdDetermineBonds and not (_traj_playing or _freq_playing):
                try:
                    charge = self.parser.data.get("charge", 0)
                    rdDetermineBonds.DetermineConnectivity(mol)
                    rdDetermineBonds.DetermineBondOrders(mol, charge=charge)
                except Exception:
                    pass  # Non-fatal; some charge states are unsupported

            final_mol = mol.GetMol()

            # Set as current molecule for export functionality
            if hasattr(self.mw, "current_mol"):
                self.mw.current_mol = final_mol

            if hasattr(self.mw, "view_3d_manager") and hasattr(
                self.mw.view_3d_manager, "draw_molecule_3d"
            ):
                self.mw.view_3d_manager.draw_molecule_3d(final_mol)

            # Sync with main window features
            if hasattr(self.mw, "ui_manager"):
                # Treat ORCA molecules as XYZ-derived (fixed geometry)
                self.mw.is_xyz_derived = True

                # Enter 3D mode which enables export/analysis buttons and hides 2D panel
                if hasattr(self.mw.ui_manager, "_enter_3d_viewer_ui_mode"):
                    self.mw.ui_manager._enter_3d_viewer_ui_mode()
                elif hasattr(self.mw.ui_manager, "_enable_3d_features"):
                    self.mw.ui_manager._enable_3d_features(True)
                    if hasattr(self.mw.ui_manager, "minimize_2d_panel"):
                        self.mw.ui_manager.minimize_2d_panel()

            # Reset 3D camera to fit molecule
            if hasattr(self.mw, "view_3d_manager") and hasattr(
                self.mw.view_3d_manager, "plotter"
            ):
                try:
                    self.mw.view_3d_manager.plotter.reset_camera()
                    self.mw.view_3d_manager.plotter.render()
                except Exception as _e:
                    logging.warning("[gui.py:592] silenced: %s", _e)
        except Exception as e:
            # self.logger.error(f"Error loading 3D: {e}")
            # print(f"Error loading 3D: {e}")
            import traceback

            traceback.print_exc()

    def show_mo_analyzer(self):
        self.load_structure_3d()
        mo_coeffs = self.parser.data.get("mo_coeffs", None)
        orb_energies = self.parser.data.get("orbital_energies", None)

        data_to_show = mo_coeffs if mo_coeffs else orb_energies

        if not data_to_show:
            QMessageBox.warning(
                self, "No Data", "No Molecular Orbital coefficients or energies found."
            )
            return

        if getattr(self, "mo_dlg", None) is not None and self.mo_dlg is not None:
            self.mo_dlg.close()

        self.mo_dlg = MODialog(self, data_to_show)
        self.mo_dlg.show()

    def show_freq(self):
        self.load_structure_3d()  # Reset to final structure before opening
        freqs = self.parser.data.get("frequencies", [])
        if not freqs:
            QMessageBox.warning(self, "No Data", "No frequency data found.")
            return
        atoms = self.parser.data.get("atoms", [])
        coords = self.parser.data.get("coords", [])

        # Ensure only one instance is open
        if getattr(self, "freq_dlg", None) is not None and self.freq_dlg is not None:
            self.freq_dlg.close()

        self.freq_dlg = FrequencyDialog(self.mw, freqs, atoms, coords)
        self.freq_dlg.show()

    def show_trajectory(self):
        data = self.parser.data.get("scan_steps", [])
        if not data:
            QMessageBox.warning(
                self, "No Info", "No trajectory steps (Optimization / Scan) found."
            )
            return
        if getattr(self, "traj_dlg", None) is not None and self.traj_dlg is not None:
            self.traj_dlg.close()
        charge = self.parser.data.get("charge", 0)
        base_dir = os.path.dirname(self.file_path) if self.file_path else None
        parsed_trj = self.parser.data.get("neb_trj_file", None)
        self.traj_dlg = TrajectoryResultDialog(
            self.mw,
            data,
            charge=charge,
            title="Trajectory / NEB Analysis",
            base_dir=base_dir,
            output_path=self.file_path,
            predicted_trj=parsed_trj,
        )
        self.traj_dlg.show()

    def show_forces(self):
        grads = self.parser.data.get("gradients", [])
        has_scan = bool(self.parser.data.get("scan_steps"))

        if not grads and not has_scan:
            QMessageBox.warning(
                self, "No Info", "No cartesian gradients or optimization steps found."
            )
            return

        if (
            getattr(self, "forces_dlg", None) is not None
            and self.forces_dlg is not None
        ):
            self.forces_dlg.close()
        # Open Force Viewer with trajectory support
        self.forces_dlg = ForceViewerDialog(self, grads, parser=self.parser)
        self.forces_dlg.show()

    def show_thermal(self):
        self.load_structure_3d()
        data = self.parser.data.get("thermal", {})
        if not data:
            QMessageBox.warning(self, "No Info", "No thermochemistry section found.")
            return
        if (
            getattr(self, "thermal_dlg", None) is not None
            and self.thermal_dlg is not None
        ):
            self.thermal_dlg.close()
        self.thermal_dlg = ThermalTableDialog(self, data)
        self.thermal_dlg.show()

    def show_tddft(self):
        self.load_structure_3d()
        excitations = self.parser.data.get("tddft", [])
        if not excitations:
            QMessageBox.warning(
                self, "No Analysis", "No TDDFT/TDA excitation energies found."
            )
            return

        if getattr(self, "tddft_dlg", None) is not None and self.tddft_dlg is not None:
            self.tddft_dlg.close()

        self.tddft_dlg = TDDFTDialog(self, excitations)
        self.tddft_dlg.show()

    def show_dipole(self):
        self.load_structure_3d()
        d = self.parser.data.get("dipoles", None)
        if not d:
            QMessageBox.warning(self, "No Info", "No dipole moment found.")
            return
        if (
            getattr(self, "dipole_dlg", None) is not None
            and self.dipole_dlg is not None
        ):
            self.dipole_dlg.close()
        self.dipole_dlg = DipoleDialog(self, d)
        self.dipole_dlg.show()

    def show_charges(self):
        self.load_structure_3d()
        charges = self.parser.data.get("charges", {})
        if not charges:
            QMessageBox.warning(self, "No Info", "No atomic charges found.")
            return
        if (
            getattr(self, "charges_dlg", None) is not None
            and self.charges_dlg is not None
        ):
            self.charges_dlg.close()
        self.charges_dlg = ChargeDialog(self, charges)
        self.charges_dlg.show()

    def show_nmr(self):
        self.load_structure_3d()
        data = self.parser.data.get("nmr_shielding", [])
        couplings = self.parser.data.get("nmr_couplings", [])
        if not data:
            QMessageBox.warning(
                self, "No Info", "No NMR chemical shielding data found."
            )
            return
        if getattr(self, "nmr_dlg", None) is not None and self.nmr_dlg is not None:
            self.nmr_dlg.close()
        self.nmr_dlg = NMRDialog(
            self, data, couplings=couplings, file_path=self.file_path
        )
        self.nmr_dlg.show()

    def show_scf_trace(self):
        self.load_structure_3d()
        data = self.parser.data.get("scf_traces", [])
        if not data:
            QMessageBox.warning(self, "No Info", "No SCF energy trace data found.")
            return
        if getattr(self, "scf_dlg", None) is not None and self.scf_dlg is not None:
            self.scf_dlg.close()
        self.scf_dlg = SCFTraceDialog(self, data)
        self.scf_dlg.show()

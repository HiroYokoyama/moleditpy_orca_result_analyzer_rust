from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QRadioButton,
    QDoubleSpinBox,
    QCheckBox,
    QPushButton,
    QFileDialog,
    QMessageBox,
    QGroupBox,
    QComboBox,
)
import os
import json

try:
    from .spectrum_widget import SpectrumWidget
except ImportError:
    try:
        from spectrum_widget import SpectrumWidget
    except Exception:
        SpectrumWidget = None
from .utils import get_default_export_path
import logging


class TDDFTDialog(QDialog):
    def __init__(self, parent, excitations):
        super().__init__(parent)
        self.setWindowTitle("TDDFT Spectrum")
        self.resize(700, 700)
        self.excitations = excitations
        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")
        self.prev_sigma_unit_idx = 0  # 0: cm-1, 1: eV

        main_layout = QVBoxLayout(self)

        # 0. Spectrum Display (Existing Widget)
        if SpectrumWidget:
            # Data format for widget: List of dicts.
            # Our excitations have 'energy_nm' and 'osc'.
            self.spectrum = SpectrumWidget(
                self.excitations,
                x_key="energy_nm",
                y_key="osc",
                x_unit="Wavelength (nm)",
                sigma=3000.0,
            )
            self.spectrum.show_legend = False
            main_layout.addWidget(self.spectrum)
            self.spectrum.clicked.connect(self.on_spectrum_click)
        else:
            main_layout.addWidget(QLabel("SpectrumWidget not available."))
            self.spectrum = None

        # 1. Spectrum Settings Section
        settings_group = QGroupBox("Spectrum Settings")
        settings_vbox = QVBoxLayout(settings_group)

        # Row 1: Type & Gauge
        row1_layout = QHBoxLayout()
        row1_layout.addWidget(QLabel("Spectrum Type:"))
        self.radio_abs = QRadioButton("Absorption")
        self.radio_abs.setChecked(True)
        self.radio_abs.toggled.connect(self.switch_spectrum_type)
        row1_layout.addWidget(self.radio_abs)

        self.radio_cd = QRadioButton("CD")
        self.radio_cd.toggled.connect(self.switch_spectrum_type)
        row1_layout.addWidget(self.radio_cd)

        row1_layout.addWidget(QLabel(" | Gauge:"))
        self.combo_gauge = QComboBox()
        self.combo_gauge.addItems(["Length (Electric)", "Velocity"])
        self.combo_gauge.currentIndexChanged.connect(self.switch_spectrum_type)
        self.combo_gauge.setEnabled(True)
        self.combo_gauge.setFixedWidth(130)
        row1_layout.addWidget(self.combo_gauge)
        row1_layout.addStretch()
        settings_vbox.addLayout(row1_layout)

        # Row 2: Sigma & Sticks
        row2_layout = QHBoxLayout()
        row2_layout.addWidget(QLabel("Broadening (FWHM):"))
        # Default Broadening
        self.spin_sigma = QDoubleSpinBox()
        self.spin_sigma.setRange(0.1, 10000.0)
        self.spin_sigma.setValue(3000.0)
        self.spin_sigma.setSingleStep(10.0)
        self.spin_sigma.valueChanged.connect(self.update_spectrum_sigma)
        row2_layout.addWidget(self.spin_sigma)

        self.combo_sigma_unit = QComboBox()
        self.combo_sigma_unit.addItems(["cm⁻¹", "eV"])
        self.combo_sigma_unit.setCurrentIndex(0)  # Default to cm-1
        self.combo_sigma_unit.currentIndexChanged.connect(self.on_sigma_unit_changed)
        row2_layout.addWidget(self.combo_sigma_unit)

        row2_layout.addSpacing(20)
        self.chk_sticks = QCheckBox("Show Transitions (Sticks)")
        self.chk_sticks.setChecked(True)
        if self.spectrum:
            self.chk_sticks.stateChanged.connect(self.spectrum.set_sticks)
        row2_layout.addWidget(self.chk_sticks)

        self.chk_markers = QCheckBox("Markers")
        self.chk_markers.setChecked(True)
        if self.spectrum:
            self.chk_markers.stateChanged.connect(self.spectrum.set_markers)
        row2_layout.addWidget(self.chk_markers)

        row2_layout.addStretch()

        # Reset Button
        self.btn_reset = QPushButton("Reset Defaults")
        self.btn_reset.clicked.connect(self.reset_defaults)
        row2_layout.addWidget(self.btn_reset)

        settings_vbox.addLayout(row2_layout)

        # Row 3: Physical Broadening
        row3_layout = QHBoxLayout()
        self.chk_physical = QCheckBox("Physical (Area-preserving) Broadening")
        self.chk_physical.setChecked(True)
        self.chk_physical.toggled.connect(self.switch_spectrum_type)
        row3_layout.addWidget(self.chk_physical)
        row3_layout.addStretch()
        settings_vbox.addLayout(row3_layout)

        main_layout.addWidget(settings_group)

        # 2. Axis Controls Section
        axis_group = QGroupBox("Axis Scale")
        axis_layout = QVBoxLayout(axis_group)

        # X-Range
        x_row = QHBoxLayout()
        self.chk_auto_x = QCheckBox("Auto Wavelength Range")
        self.chk_auto_x.setChecked(True)
        self.chk_auto_x.stateChanged.connect(self.toggle_auto_x)
        x_row.addWidget(self.chk_auto_x)

        x_row.addWidget(QLabel("Range (nm):"))
        self.spin_x_min = QDoubleSpinBox()
        self.spin_x_min.setRange(0, 50000)
        self.spin_x_min.setValue(100)
        self.spin_x_min.setEnabled(False)
        self.spin_x_min.setFixedWidth(80)
        self.spin_x_min.valueChanged.connect(self.update_x_range)
        x_row.addWidget(self.spin_x_min)

        x_row.addWidget(QLabel("-"))
        self.spin_x_max = QDoubleSpinBox()
        self.spin_x_max.setRange(0, 50000)
        self.spin_x_max.setValue(800)
        self.spin_x_max.setEnabled(False)
        self.spin_x_max.setFixedWidth(80)
        self.spin_x_max.valueChanged.connect(self.update_x_range)
        x_row.addWidget(self.spin_x_max)
        x_row.addStretch()
        axis_layout.addLayout(x_row)

        # Y-Range
        y_row = QHBoxLayout()
        self.chk_auto_y = QCheckBox("Auto Intensity Range  ")
        self.chk_auto_y.setChecked(True)
        self.chk_auto_y.stateChanged.connect(self.toggle_auto_y)
        y_row.addWidget(self.chk_auto_y)

        y_row.addWidget(QLabel("Range:"))
        self.spin_y_min = QDoubleSpinBox()
        self.spin_y_min.setRange(-10000, 10000)
        self.spin_y_min.setValue(0)
        self.spin_y_min.setEnabled(False)
        self.spin_y_min.setFixedWidth(80)
        self.spin_y_min.valueChanged.connect(self.update_range)
        y_row.addWidget(self.spin_y_min)

        y_row.addWidget(QLabel("-"))
        self.spin_y_max = QDoubleSpinBox()
        self.spin_y_max.setRange(-10000, 10000)
        self.spin_y_max.setValue(1.0)
        self.spin_y_max.setEnabled(False)
        self.spin_y_max.setFixedWidth(80)
        self.spin_y_max.valueChanged.connect(self.update_range)
        y_row.addWidget(self.spin_y_max)
        y_row.addStretch()
        axis_layout.addLayout(y_row)

        main_layout.addWidget(axis_group)

        # 3. Actions Section
        action_layout = QHBoxLayout()

        self.btn_png = QPushButton("Export Image (PNG)")
        self.btn_png.clicked.connect(self.save_png)
        action_layout.addWidget(self.btn_png)

        self.btn_csv = QPushButton("Export Data (CSV)")
        self.btn_csv.clicked.connect(self.save_csv)
        action_layout.addWidget(self.btn_csv)

        self.btn_sticks = QPushButton("Export Sticks (CSV)")
        self.btn_sticks.clicked.connect(self.save_sticks)
        action_layout.addWidget(self.btn_sticks)

        self.btn_report = QPushButton("Export Full Data (.txt)")
        self.btn_report.clicked.connect(self.save_orca_report)
        action_layout.addWidget(self.btn_report)

        action_layout.addStretch()

        self.btn_close = QPushButton("Close")
        self.btn_close.setFixedWidth(100)
        self.btn_close.clicked.connect(self.close)
        action_layout.addWidget(self.btn_close)

        main_layout.addLayout(action_layout)
        self.load_settings()
        self.switch_spectrum_type()

    def closeEvent(self, event):
        self.save_settings()
        super().closeEvent(event)

    def load_settings(self):
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r") as f:
                    all_settings = json.load(f)

                settings = all_settings.get("tddft_settings", {})

                # Block signals during loading to prevent unintended conversion logic
                self.spin_sigma.blockSignals(True)
                self.combo_sigma_unit.blockSignals(True)

                if "sigma" in settings:
                    self.spin_sigma.setValue(float(settings["sigma"]))

                if "sigma_unit_idx" in settings:
                    idx = int(settings["sigma_unit_idx"])
                    self.combo_sigma_unit.setCurrentIndex(idx)
                    self.prev_sigma_unit_idx = idx
                    # Set initial step size
                    self.spin_sigma.setSingleStep(0.1 if idx == 1 else 10.0)

                if "show_sticks" in settings:
                    self.chk_sticks.setChecked(bool(settings["show_sticks"]))

                if "physical" in settings:
                    self.chk_physical.setChecked(bool(settings["physical"]))

                self.spin_sigma.blockSignals(False)
                self.combo_sigma_unit.blockSignals(False)

            except Exception as e:
                self.spin_sigma.blockSignals(False)
                self.combo_sigma_unit.blockSignals(False)
                print(f"Error loading TDDFT settings: {e}")

    def update_spectrum_sigma(self):
        if not self.spectrum:
            return
        val = self.spin_sigma.value()  # This is FWHM
        unit_idx = self.combo_sigma_unit.currentIndex()  # 0: cm-1, 1: eV

        # Convert val to cm-1 if eV
        # 1 eV = 8065.544 cm-1
        fwhm_cm = val
        if unit_idx == 1:  # eV
            fwhm_cm = val * 8065.544

        # Convert FWHM to Sigma (FWHM = 2 * sqrt(2 * ln 2) * sigma approx 2.355 * sigma)
        self.spectrum.sigma = fwhm_cm / 2.355
        self.spectrum.broaden_in_energy = True  # Both cm-1 and eV are Energy
        self.spectrum.update()

    def on_sigma_unit_changed(self):
        curr_idx = self.combo_sigma_unit.currentIndex()
        if curr_idx == self.prev_sigma_unit_idx:
            return

        val = self.spin_sigma.value()
        # 1 eV = 8065.544 cm-1
        if self.prev_sigma_unit_idx == 0 and curr_idx == 1:  # cm-1 to eV
            new_val = val / 8065.544
        elif self.prev_sigma_unit_idx == 1 and curr_idx == 0:  # eV to cm-1
            new_val = val * 8065.544
        else:
            new_val = val

        self.spin_sigma.blockSignals(True)
        self.spin_sigma.setValue(new_val)
        # Update step size: 0.1 for eV, 10.0 for cm-1
        self.spin_sigma.setSingleStep(0.1 if curr_idx == 1 else 10.0)
        self.spin_sigma.blockSignals(False)

        self.prev_sigma_unit_idx = curr_idx
        self.update_spectrum_sigma()
        self.switch_spectrum_type()

    def save_settings(self):
        all_settings = {}
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r") as f:
                    all_settings = json.load(f)
            except Exception as _e:
                logging.warning("[tddft_analysis.py:291] silenced: %s", _e)

        tddft_settings = {
            "sigma": self.spin_sigma.value(),
            "sigma_unit_idx": self.combo_sigma_unit.currentIndex(),
            "show_sticks": self.chk_sticks.isChecked(),
            "physical": self.chk_physical.isChecked(),
        }

        all_settings["tddft_settings"] = tddft_settings

        try:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.settings_file), exist_ok=True)
            with open(self.settings_file, "w") as f:
                json.dump(all_settings, f, indent=2)

            # log to status bar if possible
            if hasattr(self.parent(), "mw") and self.parent().mw:
                self.parent().mw.statusBar().showMessage("TDDFT settings saved.", 3000)
        except Exception as e:
            print(f"Error saving TDDFT settings: {e}")

    def toggle_auto_y(self):
        is_auto = self.chk_auto_y.isChecked()
        self.spin_y_min.setEnabled(not is_auto)
        self.spin_y_max.setEnabled(not is_auto)
        if is_auto and self.spectrum:
            self.spectrum.set_auto_range()
        else:
            self.update_range()

    def update_range(self):
        if self.chk_auto_y.isChecked() or not self.spectrum:
            return
        ymin = self.spin_y_min.value()
        ymax = self.spin_y_max.value()
        self.spectrum.set_y_range(ymin, ymax)

    def toggle_auto_x(self):
        is_auto = self.chk_auto_x.isChecked()
        self.spin_x_min.setEnabled(not is_auto)
        self.spin_x_max.setEnabled(not is_auto)
        if is_auto and self.spectrum:
            self.spectrum.set_auto_x_range()
        else:
            self.update_x_range()

    def update_x_range(self):
        if self.chk_auto_x.isChecked() or not self.spectrum:
            return
        xmin = self.spin_x_min.value()
        xmax = self.spin_x_max.value()
        self.spectrum.set_x_range(xmin, xmax)

    def save_png(self):
        if not self.spectrum:
            return
        default_path = get_default_export_path(
            self.parent().file_path, suffix="_tddft", extension=".png"
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Graph", default_path, "Images (*.png)"
        )
        if path:
            self.spectrum.save_png(path)

    def save_csv(self):
        if not self.spectrum:
            return
        default_path = get_default_export_path(
            self.parent().file_path, suffix="_tddft", extension=".csv"
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Data", default_path, "CSV Files (*.csv)"
        )
        if path:
            success = self.spectrum.save_csv(path)
            if success:
                # print(f"Data saved to {path}")
                # QMessageBox.information(self, "Saved", f"Data saved to:\n{path}")
                if hasattr(self.parent(), "mw") and self.parent().mw:
                    self.parent().mw.statusBar().showMessage(
                        f"Data saved to {path}", 5000
                    )
            else:
                QMessageBox.warning(self, "Error", "Failed to save CSV.")

    def save_sticks(self):
        if not self.spectrum:
            return
        default_path = get_default_export_path(
            self.parent().file_path, suffix="_tddft_sticks", extension=".csv"
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Sticks", default_path, "CSV Files (*.csv)"
        )
        if path:
            success = self.spectrum.save_sticks_csv(path)
            if success:
                # print(f"Stick data saved to {path}")
                # QMessageBox.information(self, "Exported", f"Stick data saved to:\n{path}")
                if hasattr(self.parent(), "mw") and self.parent().mw:
                    self.parent().mw.statusBar().showMessage(
                        f"Stick data saved to {path}", 5000
                    )
            else:
                QMessageBox.warning(self, "Error", "Failed to export stick data.")

    def save_orca_report(self):
        if not self.excitations:
            QMessageBox.warning(self, "No Data", "No excitation data to export.")
            return

        default_path = get_default_export_path(
            self.parent().file_path, suffix="_tddft_report", extension=".txt"
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Report", default_path, "Text Files (*.txt)"
        )
        if not path:
            return

        import datetime

        try:
            with open(path, "w", encoding="utf-8") as f:
                # --- Header ---
                f.write("=" * 80 + "\n")
                f.write("             ORCA TD-DFT / TDA EXCITATION SPECTRUM ANALYSIS\n")
                f.write("=" * 80 + "\n")
                f.write(
                    f"Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                )
                f.write("-" * 80 + "\n")
                f.write("UNITS & DEFINITIONS:\n")
                f.write("  Energy      : eV, nm, cm^-1\n")
                f.write("  f (Osc)     : Oscillator Strength (Dimensionless)\n")
                f.write("  R (Rot)     : Rotatory Strength   (10^-40 esu^2 cm^2)\n")
                f.write("  Gauge       : L = Length (Dipole), V = Velocity\n")
                f.write("-" * 80 + "\n\n")

                for item in self.excitations:
                    state = item.get("state", "?")

                    # Skip Ground State if labeled '0'
                    if str(state) == "0":
                        continue

                    # --- Prepare Energy Values ---
                    energy_ev = item.get("energy_ev", 0.0)
                    energy_nm = item.get("energy_nm", 0.0)
                    energy_cm = item.get("energy_cm", 0.0)

                    # Auto-fill missing energy units
                    if energy_cm == 0.0 and energy_nm > 0:
                        energy_cm = 10000000.0 / energy_nm
                    if energy_ev == 0.0 and energy_nm > 0:
                        energy_ev = 1239.84193 / energy_nm

                    # --- Prepare Strength Values (Length vs Velocity) ---
                    # Oscillator Strength
                    # 変更: デフォルトを 0.0 ではなく None にして区別する
                    f_len = item.get("osc_len", item.get("osc", None))
                    f_vel = item.get("osc_vel", None)

                    # Rotatory Strength
                    r_len = item.get("rot_len", item.get("rotatory_strength", None))
                    r_vel = item.get("rot_vel", None)

                    # Helper for formatting (None -> "N/A")
                    def fmt_val(v):
                        return f"{v:>9.6f}" if v is not None else "      N/A"

                    # --- Writing the Block ---
                    f.write(
                        f"STATE {state:>3}  --------------------------------------------------------\n"
                    )
                    f.write(
                        f"  Energy : {energy_ev:>8.4f} eV  |  {energy_nm:>8.2f} nm  |  {energy_cm:>8.1f} cm^-1\n"
                    )

                    f.write(f"  f (Osc): L={fmt_val(f_len)}  |  V={fmt_val(f_vel)}\n")

                    # Rotatory Strength (CD計算データがある場合のみ行自体を表示)
                    if r_len is not None or r_vel is not None:
                        f.write(
                            f"  R (Rot): L={fmt_val(r_len)}  |  V={fmt_val(r_vel)}\n"
                        )

                    f.write("\n  Major Transitions:\n")

                    # Transition Details
                    transitions = item.get("transitions", [])
                    if isinstance(transitions, list):
                        if not transitions:
                            f.write("    (No transition details available)\n")
                        for trans in transitions:
                            # Clean up formatting if needed
                            f.write(f"    {trans}\n")
                    else:
                        f.write(f"    {transitions}\n")

                    f.write("\n")  # Empty line between states

            if hasattr(self.parent(), "mw") and self.parent().mw:
                self.parent().mw.statusBar().showMessage(
                    f"Report saved to {path}", 5000
                )
            else:
                QMessageBox.information(self, "Exported", f"Report saved to:\n{path}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save report:\n{e}")

    def reset_defaults(self):
        self.spin_sigma.blockSignals(True)
        self.combo_sigma_unit.blockSignals(True)
        self.chk_sticks.blockSignals(True)
        self.chk_physical.blockSignals(True)

        self.spin_sigma.setValue(3000.0)
        self.spin_sigma.setRange(0.1, 10000.0)
        self.spin_sigma.setSingleStep(10.0)
        self.combo_sigma_unit.setCurrentIndex(0)  # cm-1
        self.prev_sigma_unit_idx = 0
        self.chk_sticks.setChecked(True)
        self.chk_physical.setChecked(True)

        self.spin_sigma.blockSignals(False)
        self.combo_sigma_unit.blockSignals(False)
        self.chk_sticks.blockSignals(False)
        self.chk_physical.blockSignals(False)

        self.switch_spectrum_type()

    def switch_spectrum_type(self):
        """
        Calculates processed_y (Intensity or Area) based on Physical/Gauge settings.
        Correctly handles CD energy dependence (R * Energy).
        """
        if not self.excitations:
            return
        if not self.spectrum:
            return

        is_physical = self.chk_physical.isChecked()
        is_cd_mode = self.radio_cd.isChecked()  # True if CD
        gauge_mode = self.combo_gauge.currentIndex()  # 0: Length, 1: Velocity

        target_key = ""
        y_unit_main = ""
        y_unit_sticks = ""

        # --- 1. Key Selection & Units ---
        sample = self.excitations[0] if self.excitations else {}

        if not is_cd_mode:  # Absorption
            if gauge_mode == 1:  # Velocity
                target_key = "osc_vel" if "osc_vel" in sample else "osc"
            else:  # Length
                target_key = "osc_len" if "osc_len" in sample else "osc"

            if is_physical:
                y_unit_main = "ε (M⁻¹ cm⁻¹)"
                y_unit_sticks = "f (osc)"
            else:
                y_unit_main = "Intensity (Arb.)"
                y_unit_sticks = "f (osc)"

        else:  # CD
            if gauge_mode == 1:  # Velocity
                target_key = "rot_vel" if "rot_vel" in sample else "rotatory_strength"
            else:  # Length
                target_key = "rot_len" if "rot_len" in sample else "rotatory_strength"

            if is_physical:
                y_unit_main = "Δε (M⁻¹ cm⁻¹)"
                y_unit_sticks = "R (10⁻⁴⁰ cgs)"
            else:
                y_unit_main = "Intensity (Arb.)"
                y_unit_sticks = "R (10⁻⁴⁰ cgs)"

        # --- 2. Constants ---
        # Absorption: Integral(epsilon) dE = 2.315e8 * f
        ABS_FACTOR = 2.315e8

        # CD: Integral(DeltaEps) dE ~ R * Energy.
        # Factor derived from R(cgs) relation.
        # R is usually 10^-40 cgs. The factor to map to DeltaEps integral is approx 0.04355
        CD_FACTOR = 0.04355

        # --- 3. Process Data ---
        processed_data = []

        for item in self.excitations:
            new_item = item.copy()
            val = item.get(target_key, 0.0)
            if val is None:
                val = 0.0

            if is_physical:
                if not is_cd_mode:
                    # Absorption: Area is proportional to f (constant over energy)
                    new_item["processed_y"] = val * ABS_FACTOR
                else:
                    # CD: Area is proportional to R * Energy(wavenumber)
                    # We need the energy in cm-1 for correct scaling
                    e_cm = item.get("energy_cm", 0.0)
                    if not isinstance(e_cm, (int, float)):
                        e_cm = 0.0

                    # Use nm to calculate cm-1 if missing
                    if e_cm <= 0:
                        e_nm = item.get("energy_nm", 0.0)
                        if isinstance(e_nm, (int, float)) and e_nm > 0:
                            e_cm = 1e7 / e_nm
                        else:
                            e_cm = 0.0

                    # 修正点: R値にエネルギーを掛けてスケーリング
                    new_item["processed_y"] = val * e_cm * CD_FACTOR
            else:
                # Non-physical (Height = Value)
                new_item["processed_y"] = val

            processed_data.append(new_item)

        # --- 4. Update Widget ---
        self.spectrum.data_list = processed_data
        self.spectrum.y_key = "processed_y"
        self.spectrum.y_key_sticks = target_key  # Sticks always show raw f or R

        self.spectrum.y_unit = y_unit_main
        self.spectrum.y_unit_sticks = y_unit_sticks

        self.spectrum.normalization_mode = "area" if is_physical else "height"
        self.spectrum.broaden_in_energy = (
            True  # Always broaden in Energy for physical correctness
        )

        self.spectrum.set_dual_axis(is_physical)
        self.update_spectrum_sigma()
        self.spectrum.update()

    def on_spectrum_click(self, item):
        """Handle click on spectrum peak"""
        if item is None:
            return
        state = item.get("state", "?")
        energy_ev = item.get("energy_ev", 0.0)
        energy_nm = item.get("energy_nm", 0.0)

        is_cd = self.radio_cd.isChecked()
        gauge_mode = self.combo_gauge.currentIndex()  # 0: Length, 1: Velocity

        msg = f"Transition to Excited State {state}\n"
        msg += f"Energy: {energy_ev:.4f} eV ({energy_nm:.2f} nm)\n"
        msg += "-" * 30 + "\n"

        # Show Oscillator Strengths
        osc_len = item.get("osc_len", item.get("osc", 0.0))
        osc_vel = item.get("osc_vel", 0.0)

        # Show Rotatory Strengths
        rot_len = item.get("rot_len", item.get("rotatory_strength", 0.0))
        rot_vel = item.get("rot_vel", 0.0)

        if not is_cd:
            msg += "Oscillator Strength (f):\n"
            msg += (
                f"  - Length   : {osc_len:.6f}"
                + (" *" if gauge_mode == 0 else "")
                + "\n"
            )
            msg += (
                f"  - Velocity : {osc_vel:.6f}"
                + (" *" if gauge_mode == 1 else "")
                + "\n"
            )
        else:
            msg += "Rotatory Strength (R) [10⁻⁴⁰ esu² cm²]:\n"
            msg += (
                f"  - Length   : {rot_len:.6f}"
                + (" *" if gauge_mode == 0 else "")
                + "\n"
            )
            msg += (
                f"  - Velocity : {rot_vel:.6f}"
                + (" *" if gauge_mode == 1 else "")
                + "\n"
            )

        msg += "-" * 30 + "\n"

        # Use instantiated QMessageBox for selectable text
        from PyQt6.QtCore import Qt

        msg_box = QMessageBox(self)
        msg_box.setWindowTitle(f"Transition Details - State {state}")

        # Add scroll area if message is too long? QMessageBox doesn't support scrolling easily.
        # But for 20 lines it's fine.
        # Let's increase the limit or make it nicer.

        if "transitions" in item:
            trans_data = item["transitions"]
            if trans_data and len(trans_data) > 0:
                msg += "\nOrbital Contributions:\n"
                # Show ALL transitions
                msg += "\n".join(trans_data)
            else:
                msg += f"\n(No detailed orbital contributions found - data type: {type(trans_data)}, len: {len(trans_data) if trans_data else 0})"
        else:
            msg += "\n(No 'transitions' key in item)"

        msg_box.setText(msg)
        msg_box.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        msg_box.exec()

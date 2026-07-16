from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QPushButton,
    QApplication,
    QCheckBox,
    QAbstractItemView,
    QFileDialog,
)
import os
import json
import csv
import logging


class ThermalTableDialog(QDialog):
    def __init__(self, parent, data):
        super().__init__(parent)
        self.setWindowTitle("Thermochemistry")
        self.resize(550, 460)
        self.resize(550, 460)
        self.data = data
        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")

        layout = QVBoxLayout(self)

        # Checkbox for details
        self.chk_details = QCheckBox("Show Detailed Values")
        self.chk_details.setChecked(False)  # Reduced by default
        self.chk_details.stateChanged.connect(self.update_table)
        layout.addWidget(self.chk_details)

        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Property", "Value"])
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        layout.addWidget(self.table)

        btn_copy = QPushButton("Copy to Clipboard")
        btn_copy.clicked.connect(self.copy_table)
        layout.addWidget(btn_copy)

        btn_csv = QPushButton("Export to CSV")
        btn_csv.clicked.connect(self.export_csv)
        layout.addWidget(btn_csv)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close)

        self.update_table()
        self.load_settings()

    def closeEvent(self, event):
        self.save_settings()
        super().closeEvent(event)

    def load_settings(self):
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r", encoding="utf-8") as f:
                    all_settings = json.load(f)

                settings = all_settings.get("thermal_settings", {})
                if "show_details" in settings:
                    self.chk_details.setChecked(bool(settings["show_details"]))

            except Exception as e:
                logging.warning("Error loading thermal settings: %s", e)

    def save_settings(self):
        all_settings = {}
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r", encoding="utf-8") as f:
                    all_settings = json.load(f)
            except Exception as _e:
                logging.warning("silenced: %s", _e)

        thermal_settings = {"show_details": self.chk_details.isChecked()}

        all_settings["thermal_settings"] = thermal_settings

        try:
            with open(self.settings_file, "w", encoding="utf-8") as f:
                json.dump(all_settings, f, indent=2)
        except Exception as e:
            logging.warning("Error saving thermal settings: %s", e)

    def update_table(self):
        show_details = self.chk_details.isChecked()
        data = self.data

        # Extract values
        temp = data.get("temperature", None)
        e_el = data.get("electronic_energy", None)
        zpe = data.get("zpe", None)
        u = data.get("thermal_energy", None)
        h = data.get("enthalpy", None)
        h_corr = data.get("enthalpy_corr", None)
        s = data.get("entropy", None)
        g = data.get("gibbs", None)
        g_corr = data.get("gibbs_corr", None)

        # Detailed
        corr_vib = data.get("corr_vib", None)
        corr_rot = data.get("corr_rot", None)
        corr_trans = data.get("corr_trans", None)
        corr_thermal_total = data.get("corr_thermal_total", None)
        corr_zpe = data.get("corr_zpe", None)
        corr_total = data.get("corr_total", None)

        s_el = data.get("s_el", None)
        s_vib = data.get("s_vib", None)
        s_rot = data.get("s_rot", None)
        s_trans = data.get("s_trans", None)

        items = []
        if e_el is not None:
            items.append(("Electronic Energy (SP)", f"{e_el:.8f} Eh"))

        if temp is not None:
            items.append(("Temperature (K)", f"{temp:.2f} K"))

        # Basic Items
        items.extend(
            [
                ("Zero Point Energy", zpe),
                ("Total Thermal Energy (U)", u),
                ("Total Enthalpy (H)", h),
                ("Enthalpy Correction (H - E_el)", h_corr),
                ("Entropy Term (T*S)", s),
                ("Gibbs Free Energy (G)", g),
                ("Gibbs Correction (G - E_el)", g_corr),
            ]
        )

        # Imaginary Frequencies
        imag_count = data.get("imaginary_freq_count", None)
        if imag_count is not None:
            items.append(("Imaginary Frequencies", str(imag_count)))

        if show_details:
            # Add Separator and Detailed Items
            items.extend(
                [
                    ("", ""),  # Separator
                    ("--- Detailed Corrections ---", ""),
                    ("Thermal Vib Correction", corr_vib),
                    ("Thermal Rot Correction", corr_rot),
                    ("Thermal Trans Correction", corr_trans),
                    ("Total Thermal Correction", corr_thermal_total),
                    ("Non-thermal (ZPE) Correction", corr_zpe),
                    ("Total Correction (ZPE + Thermal)", corr_total),
                    ("", ""),
                    ("--- Entropy Breakdown ---", ""),
                    ("Electronic Entropy (T*S_el)", s_el),
                    ("Vibrational Entropy (T*S_vib)", s_vib),
                    ("Rotational Entropy (T*S_rot)", s_rot),
                    ("Translational Entropy (T*S_trans)", s_trans),
                ]
            )

        self.table.setRowCount(len(items))

        for i, (name, val) in enumerate(items):
            self.table.setItem(i, 0, QTableWidgetItem(name))
            if isinstance(val, str):
                self.table.setItem(i, 1, QTableWidgetItem(val))
            elif val is not None:
                self.table.setItem(i, 1, QTableWidgetItem(f"{val:.8f} Eh"))
            else:
                self.table.setItem(i, 1, QTableWidgetItem("-"))

    def copy_table(self):
        text = ""
        for r in range(self.table.rowCount()):
            p = self.table.item(r, 0).text()
            v = self.table.item(r, 1).text()
            text += f"{p}\t{v}\n"
        QApplication.clipboard().setText(text)

    def export_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if path:
            try:
                with open(path, "w", newline="", encoding="utf-8") as f:
                    writer = csv.writer(f)
                    writer.writerow(["Property", "Value"])
                    for r in range(self.table.rowCount()):
                        p = self.table.item(r, 0).text()
                        v = self.table.item(r, 1).text()
                        writer.writerow([p, v])
                if self.parent() and self.parent().context:
                    self.parent().context.show_status_message(
                        f"Data exported to {path}", 5000
                    )
            except Exception:
                logging.debug("Thermochemistry export failed", exc_info=True)

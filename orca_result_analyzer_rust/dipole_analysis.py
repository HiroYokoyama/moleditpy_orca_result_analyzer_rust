import numpy as np
import os
import json
import pyvista as pv
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QCheckBox,
    QDoubleSpinBox,
    QColorDialog,
    QSpinBox,
    QGroupBox,
)
from PyQt6.QtGui import QColor
import logging


class DipoleDialog(QDialog):
    def __init__(self, parent_dlg, dipole_data):
        super().__init__(parent_dlg)
        self.setWindowTitle("Dipole Moment")
        self.resize(350, 300)
        self.parent_dlg = parent_dlg
        self.dipole_data = dipole_data
        self.arrow_actor = None
        self.arrow_actor = None
        self.arrow_color = "cyan"
        self.arrow_res = 20
        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")

        main_layout = QVBoxLayout(self)

        # 1. Info Section
        info_group = QGroupBox("Dipole Data")
        info_layout = QVBoxLayout(info_group)
        vec = dipole_data.get("vector", [0.0, 0.0, 0.0])
        mag = dipole_data.get("magnitude", 0.0)

        info_text = (
            f"Vector (X, Y, Z):<br>"
            f"<b>{vec[0]:.4f}, {vec[1]:.4f}, {vec[2]:.4f}</b> Debye<br><br>"
            f"Magnitude: <b>{mag:.4f} Debye</b>"
        )
        lbl_info = QLabel(info_text)
        lbl_info.setWordWrap(True)
        info_layout.addWidget(lbl_info)
        main_layout.addWidget(info_group)

        # 2. 3D Appearance Section
        view_group = QGroupBox("3D Visualization")
        view_layout = QVBoxLayout(view_group)

        self.chk_show = QCheckBox("Show in 3D Viewer")
        self.chk_show.setChecked(False)
        self.chk_show.stateChanged.connect(self.update_view)
        view_layout.addWidget(self.chk_show)

        self.chk_reverse = QCheckBox("Reverse Vector")
        self.chk_reverse.setChecked(True)
        self.chk_reverse.stateChanged.connect(self.update_view)
        view_layout.addWidget(self.chk_reverse)

        # Scale and Resolution in one row
        res_row = QHBoxLayout()
        res_row.addWidget(QLabel("Scale:"))
        self.spin_scale = QDoubleSpinBox()
        self.spin_scale.setRange(0.1, 10.0)
        self.spin_scale.setSingleStep(0.1)
        self.spin_scale.setValue(2.0)
        self.spin_scale.valueChanged.connect(self.update_view)
        res_row.addWidget(self.spin_scale)

        res_row.addWidget(QLabel(" Res:"))
        self.spin_res = QSpinBox()
        self.spin_res.setRange(3, 100)
        self.spin_res.setValue(20)
        self.spin_res.valueChanged.connect(self.on_res_changed)
        res_row.addWidget(self.spin_res)
        view_layout.addLayout(res_row)

        # Color picker row
        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color:"))
        self.btn_color = QPushButton()
        self.btn_color.setFixedWidth(60)
        self.btn_color.setStyleSheet(
            f"background-color: {self.arrow_color}; border: 1px solid gray; height: 20px;"
        )
        self.btn_color.clicked.connect(self.pick_color)
        color_row.addWidget(self.btn_color)
        color_row.addStretch()
        view_layout.addLayout(color_row)

        main_layout.addWidget(view_group)

        main_layout.addStretch()

        # 3. Actions
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        main_layout.addWidget(btn_close)

        self.load_settings()

    def update_view(self):
        # Clear old
        if self.arrow_actor:
            try:
                self.parent_dlg.mw.plotter.remove_actor(self.arrow_actor)
            except Exception as _e:
                logging.warning("[dipole_analysis.py:99] silenced: %s", _e)
            self.arrow_actor = None

        if not self.chk_show.isChecked():
            self.parent_dlg.mw.plotter.render()
            return

        try:
            mw = self.parent_dlg.mw
            # Center of mass calculation
            coords = self.parent_dlg.parser.data.get("coords", [])
            if coords:
                center = np.mean(coords, axis=0)
            else:
                center = np.array([0.0, 0.0, 0.0])

            vec = np.array(self.dipole_data.get("vector", [0.0, 0.0, 0.0]))
            mag = np.linalg.norm(vec)

            if mag < 1e-6:
                return

            direction = vec / mag
            if self.chk_reverse.isChecked():
                direction = -direction

            scale = self.spin_scale.value()
            length = mag * scale

            # Add Arrow using add_mesh and pv.Arrow for robustness
            arrow = pv.Arrow(
                start=center,
                direction=direction,
                scale=length,
                shaft_resolution=self.arrow_res,
                tip_resolution=self.arrow_res,
            )
            self.arrow_actor = mw.plotter.add_mesh(
                arrow, color=self.arrow_color, name="dipole_vector"
            )
            mw.plotter.render()

        except Exception as e:
            print(f"Error drawing dipole: {e}")

    def pick_color(self):
        color = QColorDialog.getColor(
            QColor(self.arrow_color), self, "Select Arrow Color"
        )
        if color.isValid():
            self.arrow_color = color.name()
            self.btn_color.setStyleSheet(f"background-color: {self.arrow_color};")
            self.update_view()

    def on_res_changed(self, val):
        self.arrow_res = val
        self.update_view()

    def closeEvent(self, event):
        if self.arrow_actor:
            try:
                self.parent_dlg.mw.plotter.remove_actor(self.arrow_actor)
                self.parent_dlg.mw.plotter.render()
            except Exception as _e:
                logging.warning("[dipole_analysis.py:152] silenced: %s", _e)
        # Clean up reference in parent
        # Clean up reference in parent
        if hasattr(self.parent_dlg, "dipole_dlg"):
            self.parent_dlg.dipole_dlg = None

        self.save_settings()
        event.accept()

    def load_settings(self):
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r") as f:
                    all_settings = json.load(f)

                settings = all_settings.get("dipole_settings", {})

                # if "scale" in settings:
                #    self.spin_scale.setValue(float(settings["scale"]))

                if "res" in settings:
                    self.spin_res.setValue(int(settings["res"]))

                if "color" in settings:
                    self.arrow_color = settings["color"]
                    self.btn_color.setStyleSheet(
                        f"background-color: {self.arrow_color}; border: 1px solid gray; height: 20px;"
                    )

                if "show" in settings:
                    self.chk_show.setChecked(bool(settings["show"]))

                if "reverse" in settings:
                    self.chk_reverse.setChecked(bool(settings["reverse"]))

            except Exception as e:
                print(f"Error loading dipole settings: {e}")

    def save_settings(self):
        all_settings = {}
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r") as f:
                    all_settings = json.load(f)
            except Exception:
                pass  # settings file may be empty or corrupt; start fresh

        dipole_settings = {
            # "scale": self.spin_scale.value(),
            "res": self.spin_res.value(),
            "color": self.arrow_color,
            "show": self.chk_show.isChecked(),
            "reverse": self.chk_reverse.isChecked(),
        }

        all_settings["dipole_settings"] = dipole_settings

        try:
            with open(self.settings_file, "w") as f:
                json.dump(all_settings, f, indent=2)
        except Exception as e:
            print(f"Error saving dipole settings: {e}")

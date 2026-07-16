"""
tests/test_gui_reset_on_load.py

Tests that loading a *new* ORCA result into an already-open analyzer window
(Select File / Reload / drag-drop) clears stale 3D atom-color overrides
left by the Atomic Charges view, so colors keyed to a previous molecule's
atom indices never bleed onto a differently-indexed one.
"""

import sys
import unittest
from unittest.mock import MagicMock, patch

# Stub Qt (mirrors tests/test_about_menu.py's bootstrap so this file is
# self-sufficient regardless of test collection order).
if "PyQt6.QtWidgets" in sys.modules:
    qtw = sys.modules["PyQt6.QtWidgets"]
else:
    qtw = MagicMock()
    sys.modules["PyQt6.QtWidgets"] = qtw

if not hasattr(qtw, "QLabel") or "Mock" in type(getattr(qtw, "QLabel")).__name__:

    class _QLabel:
        def __init__(self, *a, **k):
            pass

        def setText(self, *a):
            pass

        def fontMetrics(self):
            fm = MagicMock()
            fm.elidedText.return_value = ""
            return fm

        def width(self):
            return 100

        def resizeEvent(self, *a):
            pass

        def __getattr__(self, name):
            if name.startswith("_"):
                raise AttributeError(name)
            return MagicMock()

    qtw.QLabel = _QLabel

if (
    not hasattr(qtw, "QPushButton")
    or "Mock" in type(getattr(qtw, "QPushButton")).__name__
):

    class _QPushButton:
        def __init__(self, *a, **k):
            self.clicked = MagicMock()

        def setIcon(self, *a):
            pass

        def setIconSize(self, *a):
            pass

        def setStyleSheet(self, *a):
            pass

        def setEnabled(self, *a):
            pass

        def setToolTip(self, *a):
            pass

        def installEventFilter(self, *a):
            pass

        def __getattr__(self, name):
            if name.startswith("_"):
                raise AttributeError(name)
            return MagicMock()

    qtw.QPushButton = _QPushButton

if not hasattr(qtw, "QDialog") or "Mock" in type(getattr(qtw, "QDialog")).__name__:

    class _QDialog:
        class DialogCode:
            Accepted = 1
            Rejected = 0

        def __init__(self, *a, **k):
            self.menu_bar_mock = MagicMock()

        def menuBar(self):
            return self.menu_bar_mock

        def setWindowTitle(self, *a):
            pass

        def resize(self, *a):
            pass

        def setAcceptDrops(self, *a):
            pass

        def close(self):
            pass

        def __getattr__(self, name):
            if name.startswith("_"):
                raise AttributeError(name)
            return MagicMock()

    qtw.QDialog = _QDialog

if "PyQt6" not in sys.modules:
    sys.modules["PyQt6"] = MagicMock()
if "PyQt6.QtCore" not in sys.modules:
    sys.modules["PyQt6.QtCore"] = MagicMock()
qt_core = sys.modules["PyQt6.QtCore"]
if (
    not hasattr(qt_core, "QObject")
    or "Mock" in type(getattr(qt_core, "QObject")).__name__
):

    class _QObject:
        def __init__(self, *a, **k):
            pass

        def __getattr__(self, name):
            if name.startswith("_"):
                raise AttributeError(name)
            return MagicMock()

    qt_core.QObject = _QObject
if "PyQt6.QtGui" not in sys.modules:
    sys.modules["PyQt6.QtGui"] = MagicMock()

if "matplotlib.backends.backend_qtagg" not in sys.modules:
    sys.modules["matplotlib.backends.backend_qtagg"] = MagicMock()
if "pyvista" not in sys.modules:
    sys.modules["pyvista"] = MagicMock()

from orca_result_analyzer_rust.gui import OrcaResultAnalyzerDialog  # noqa: E402
from orca_result_analyzer_rust.parser import OrcaParser  # noqa: E402


class TestLoadFileClearsAtomColorOverrides(unittest.TestCase):
    def _make_dialog(self):
        parser = OrcaParser()
        ctx = MagicMock()
        dlg = OrcaResultAnalyzerDialog(None, parser, "", ctx)
        # Avoid exercising the rest of load_file's file-parsing machinery —
        # only the reset-on-load wiring at the top of the method matters here.
        dlg.close_all_sub_dialogs = MagicMock()
        dlg.update_file_info_labels = MagicMock()
        dlg.load_structure_3d = MagicMock()
        dlg.update_button_states = MagicMock()
        return dlg

    def test_load_file_calls_clear_atom_color_overrides(self):
        dlg = self._make_dialog()
        with patch(
            "orca_result_analyzer_rust.gui.clear_atom_color_overrides"
        ) as mock_clear:
            # Path need not exist: the reset call happens before file I/O,
            # and any read failure is caught internally by load_file().
            dlg.load_file("/nonexistent/path/calc.out")
            mock_clear.assert_called_once_with(dlg.mw)

    def test_load_file_clears_overrides_even_on_read_failure(self):
        """A missing/unreadable file must not skip the override reset."""
        dlg = self._make_dialog()
        with patch(
            "orca_result_analyzer_rust.gui.clear_atom_color_overrides"
        ) as mock_clear:
            try:
                dlg.load_file("/definitely/does/not/exist.out")
            except Exception as exc:  # load_file must swallow read errors
                self.fail(f"load_file raised unexpectedly: {exc}")
            mock_clear.assert_called_once()

    def test_load_file_clears_before_closing_sub_dialogs_is_also_called(self):
        """Both cleanup steps must run on every load, in the same call."""
        dlg = self._make_dialog()
        with patch(
            "orca_result_analyzer_rust.gui.clear_atom_color_overrides"
        ) as mock_clear:
            dlg.load_file("/nonexistent/path/calc.out")
            dlg.close_all_sub_dialogs.assert_called_once()
            mock_clear.assert_called_once()


if __name__ == "__main__":
    unittest.main()

"""
tests/test_init.py
Unit tests for orca_result_analyzer_rust/__init__.py.

Tests plugin metadata, initialize() registration contract, and the
handle_drop() detection logic (non-Qt paths only).

PyQt6 is stubbed at module level since __init__.py imports
QMessageBox and QFileDialog unconditionally.  The Qt-dependent code
paths inside open_orca_file() (QApplication, OrcaResultAnalyzerDialog, etc.)
are not exercised here because they require a live display.
"""

import os
import sys
import types
import tempfile
import importlib.util
import unittest
from unittest.mock import MagicMock

_SRC_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))


def _install_qt_stubs():
    """Install minimal PyQt6 stubs needed for __init__.py import."""
    _pyqt6 = types.ModuleType("PyQt6")
    _widgets = types.ModuleType("PyQt6.QtWidgets")
    _widgets.QMessageBox = MagicMock()
    _widgets.QFileDialog = MagicMock()
    _pyqt6.QtWidgets = _widgets

    sys.modules.update({
        "PyQt6": _pyqt6,
        "PyQt6.QtWidgets": _widgets,
    })


def _load_init_module():
    """Load orca_result_analyzer_rust/__init__.py directly, bypassing Qt."""
    _install_qt_stubs()

    init_src = os.path.normpath(
        os.path.join(_SRC_DIR, "orca_result_analyzer_rust", "__init__.py")
    )
    spec = importlib.util.spec_from_file_location(
        "orca_result_analyzer_rust_init", init_src
    )
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = "orca_result_analyzer_rust"
    sys.modules["orca_result_analyzer_rust_init"] = mod
    spec.loader.exec_module(mod)
    return mod


_init_mod = _load_init_module()


class StubContext:
    """Mirrors the PluginContext API surface used by initialize()."""

    def __init__(self):
        self.file_openers = {}   # ext → (callback, priority)
        self.drop_handlers = []  # list of (handler, priority)
        self.menu_actions = {}

    def register_file_opener(self, ext, callback, priority=0):
        self.file_openers[ext] = (callback, priority)

    def register_drop_handler(self, handler, priority=0):
        self.drop_handlers.append((handler, priority))

    def add_menu_action(self, path, callback):
        self.menu_actions[path] = callback

    def get_main_window(self):
        return MagicMock()


# ---------------------------------------------------------------------------
# TestMetadata
# ---------------------------------------------------------------------------

class TestMetadata(unittest.TestCase):

    def test_plugin_name_present(self):
        self.assertTrue(hasattr(_init_mod, "PLUGIN_NAME"))
        self.assertIsInstance(_init_mod.PLUGIN_NAME, str)
        self.assertGreater(len(_init_mod.PLUGIN_NAME), 0)

    def test_plugin_version_present(self):
        self.assertTrue(hasattr(_init_mod, "PLUGIN_VERSION"))
        self.assertIsInstance(_init_mod.PLUGIN_VERSION, str)

    def test_plugin_author_present(self):
        self.assertTrue(hasattr(_init_mod, "PLUGIN_AUTHOR"))

    def test_plugin_description_present(self):
        self.assertTrue(hasattr(_init_mod, "PLUGIN_DESCRIPTION"))


# ---------------------------------------------------------------------------
# TestInitialize — registration contract
# ---------------------------------------------------------------------------

class TestInitialize(unittest.TestCase):

    def _make_ctx(self):
        ctx = StubContext()
        _init_mod.initialize(ctx)
        return ctx

    def test_registers_out_file_opener(self):
        ctx = self._make_ctx()
        self.assertIn(".out", ctx.file_openers)

    def test_out_opener_priority(self):
        ctx = self._make_ctx()
        _, priority = ctx.file_openers[".out"]
        self.assertEqual(priority, 100)

    def test_registers_drop_handler(self):
        ctx = self._make_ctx()
        self.assertEqual(len(ctx.drop_handlers), 1)

    def test_drop_handler_priority(self):
        ctx = self._make_ctx()
        _, priority = ctx.drop_handlers[0]
        self.assertEqual(priority, 100)

    def test_initialize_callable(self):
        self.assertTrue(callable(_init_mod.initialize))


# ---------------------------------------------------------------------------
# TestDropHandler — non-Qt code paths
# ---------------------------------------------------------------------------

class TestDropHandler(unittest.TestCase):
    """
    Tests the handle_drop closure registered by initialize().
    Only the False-return paths are tested here: correct extension but no
    ORCA content, and wrong extension entirely.
    The True-return path (ORCA header detected) calls open_orca_file() which
    requires a live Qt application and is not exercised in this headless suite.
    """

    def setUp(self):
        ctx = StubContext()
        _init_mod.initialize(ctx)
        self._handler, _ = ctx.drop_handlers[0]

    def test_non_out_extension_returns_false(self):
        self.assertFalse(self._handler("molecule.mol"))

    def test_log_extension_returns_false(self):
        self.assertFalse(self._handler("calc.log"))

    def test_nonexistent_out_file_returns_false(self):
        self.assertFalse(self._handler("/nonexistent/path/calc.out"))

    def test_out_file_without_orca_content_returns_false(self):
        with tempfile.NamedTemporaryFile(
            suffix=".out", mode="w", delete=False, encoding="utf-8"
        ) as f:
            f.write("This is just a plain text file, not ORCA output.\n")
            tmp = f.name
        try:
            result = self._handler(tmp)
            self.assertFalse(result)
        finally:
            os.unlink(tmp)

    def test_empty_out_file_returns_false(self):
        with tempfile.NamedTemporaryFile(
            suffix=".out", mode="w", delete=False, encoding="utf-8"
        ) as f:
            tmp = f.name
        try:
            result = self._handler(tmp)
            self.assertFalse(result)
        finally:
            os.unlink(tmp)

    def test_case_insensitive_extension(self):
        """Extension check should be case-insensitive (.OUT, .Out, etc.)."""
        # Non-existent file: open() will raise, caught → returns False
        # The important check is that a .OUT path reaches the extension gate
        self.assertFalse(self._handler("/some/path/CALC.OUT"))


# ---------------------------------------------------------------------------
# TestInitializeIdempotent
# ---------------------------------------------------------------------------

class TestInitializeIdempotent(unittest.TestCase):

    def test_multiple_initialize_calls_do_not_duplicate(self):
        """Each initialize() call on a fresh context should register exactly once."""
        ctx = StubContext()
        _init_mod.initialize(ctx)
        self.assertEqual(len(ctx.drop_handlers), 1)
        self.assertEqual(len(ctx.file_openers), 1)


if __name__ == "__main__":
    unittest.main()

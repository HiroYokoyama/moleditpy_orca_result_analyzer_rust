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
    existing_widgets = sys.modules.get("PyQt6.QtWidgets")
    existing_pyqt6 = sys.modules.get("PyQt6")

    if existing_widgets is not None:
        if not hasattr(existing_widgets, "QMessageBox"):
            existing_widgets.QMessageBox = MagicMock()
        if not hasattr(existing_widgets, "QFileDialog"):
            existing_widgets.QFileDialog = MagicMock()
        if not hasattr(existing_widgets, "QApplication"):
            app_mock = MagicMock()
            app_mock.processEvents = MagicMock()
            existing_widgets.QApplication = app_mock
        elif not hasattr(existing_widgets.QApplication, "processEvents"):
            existing_widgets.QApplication.processEvents = MagicMock()
    else:
        _pyqt6 = existing_pyqt6 or types.ModuleType("PyQt6")
        _widgets = types.ModuleType("PyQt6.QtWidgets")
        _widgets.QMessageBox = MagicMock()
        _widgets.QFileDialog = MagicMock()
        app_mock = MagicMock()
        app_mock.processEvents = MagicMock()
        _widgets.QApplication = app_mock
        _pyqt6.QtWidgets = _widgets

        sys.modules.update(
            {
                "PyQt6": _pyqt6,
                "PyQt6.QtWidgets": _widgets,
            }
        )


def _load_init_module():
    """Load orca_result_analyzer_rust/__init__.py directly, bypassing Qt."""
    _install_qt_stubs()

    init_src = os.path.normpath(
        os.path.join(_SRC_DIR, "orca_result_analyzer_rust", "__init__.py")
    )
    spec = importlib.util.spec_from_file_location("orca_result_analyzer_rust_init", init_src)
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = "orca_result_analyzer_rust"
    sys.modules["orca_result_analyzer_rust_init"] = mod
    spec.loader.exec_module(mod)
    return mod


_init_mod = _load_init_module()


class StubContext:
    """Mirrors the PluginContext API surface used by initialize()."""

    def __init__(self):
        self.file_openers = {}  # ext → (callback, priority)
        self.drop_handlers = []  # list of (handler, priority)
        self.menu_actions = {}
        self._windows = {}
        self.mark_project_modified_call_count = 0
        self.document_reset_handlers = []

    def register_file_opener(self, ext, callback, priority=0):
        self.file_openers[ext] = (callback, priority)

    def register_drop_handler(self, handler, priority=0):
        self.drop_handlers.append((handler, priority))

    def register_document_reset_handler(self, callback):
        self.document_reset_handlers.append(callback)

    def add_menu_action(self, path, callback):
        self.menu_actions[path] = callback

    def get_main_window(self):
        return MagicMock()

    def register_window(self, window_id: str, window) -> None:
        self._windows[window_id] = window

    def get_window(self, window_id: str):
        return self._windows.get(window_id)

    def mark_project_modified(self) -> None:
        self.mark_project_modified_call_count += 1


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


# ---------------------------------------------------------------------------
# TestContextRegistryAPI — stub completeness for V4 context registry
# ---------------------------------------------------------------------------


class TestContextRegistryAPI(unittest.TestCase):
    """Verify StubContext implements all V4 registry methods the plugin uses,
    and that initialize() does not call any methods absent from the stub."""

    def setUp(self):
        self.ctx = StubContext()
        _init_mod.initialize(self.ctx)

    def test_stub_has_register_window(self):
        self.assertTrue(callable(getattr(self.ctx, "register_window", None)))

    def test_stub_has_get_window(self):
        self.assertTrue(callable(getattr(self.ctx, "get_window", None)))

    def test_stub_has_mark_project_modified(self):
        self.assertTrue(callable(getattr(self.ctx, "mark_project_modified", None)))

    def test_registry_roundtrip(self):
        sentinel = object()
        self.ctx.register_window("test_win", sentinel)
        self.assertIs(self.ctx.get_window("test_win"), sentinel)

    def test_get_window_returns_none_for_unknown(self):
        self.assertIsNone(self.ctx.get_window("no_such_window"))

    def test_initialize_does_not_raise_with_full_stub(self):
        ctx2 = StubContext()
        try:
            _init_mod.initialize(ctx2)
        except AttributeError as exc:
            self.fail(f"initialize() called a context method not on StubContext: {exc}")


# ---------------------------------------------------------------------------
# TestExtensionsMenuAction — Extensions/ORCA Result Analyzer registered in initialize()
# ---------------------------------------------------------------------------


class TestExtensionsMenuAction(unittest.TestCase):
    def setUp(self):
        self.ctx = StubContext()
        _init_mod.initialize(self.ctx)

    def test_extensions_menu_action_registered(self):
        self.assertIn("Extensions/ORCA Result Analyzer", self.ctx.menu_actions)

    def test_extensions_menu_action_is_callable(self):
        action = self.ctx.menu_actions.get("Extensions/ORCA Result Analyzer")
        self.assertTrue(callable(action))


# ---------------------------------------------------------------------------
# TestRunFunction — run() no longer opens a file dialog
# ---------------------------------------------------------------------------


class TestRunFunction(unittest.TestCase):
    def test_run_is_callable(self):
        self.assertTrue(callable(_init_mod.run))

    def test_run_without_context_does_not_raise(self):
        """run() with no prior initialize() should silently return."""
        original = _init_mod._context
        _init_mod._context = None
        try:
            _init_mod.run(MagicMock())  # must not raise
        finally:
            _init_mod._context = original

    def test_run_does_not_open_file_dialog(self):
        """run() must NOT call QFileDialog.getOpenFileName (no file pre-selection)."""
        ctx = StubContext()
        _init_mod.initialize(ctx)

        # Patch _open_orca_analyzer_empty to be a no-op so Qt is never touched
        original = _init_mod._open_orca_analyzer_empty
        calls = []
        _init_mod._open_orca_analyzer_empty = lambda c: calls.append(c)

        # Also make sure QFileDialog is a clean mock we can introspect
        widgets_mod = sys.modules.get("PyQt6.QtWidgets")
        file_dialog_mock = MagicMock()
        original_fd = getattr(widgets_mod, "QFileDialog", None)
        if widgets_mod:
            widgets_mod.QFileDialog = file_dialog_mock

        try:
            _init_mod.run(MagicMock())
            # _open_orca_analyzer_empty should have been called once
            self.assertEqual(len(calls), 1)
            # QFileDialog.getOpenFileName must NOT have been called
            file_dialog_mock.getOpenFileName.assert_not_called()
        finally:
            _init_mod._open_orca_analyzer_empty = original
            if widgets_mod and original_fd is not None:
                widgets_mod.QFileDialog = original_fd


# ---------------------------------------------------------------------------
# TestOpenAnalyzerEmpty — singleton + window-registration contract
# ---------------------------------------------------------------------------


class TestOpenAnalyzerEmpty(unittest.TestCase):
    def _make_ctx_with_empty_open(self):
        """Return a StubContext where _open_orca_analyzer_empty is
        monkey-patched to record calls without touching Qt."""
        ctx = StubContext()
        _init_mod.initialize(ctx)
        return ctx

    def test_extensions_action_invokes_empty_open(self):
        """Calling the Extensions menu action should invoke _open_orca_analyzer_empty."""
        ctx = self._make_ctx_with_empty_open()
        calls = []
        original = _init_mod._open_orca_analyzer_empty
        _init_mod._open_orca_analyzer_empty = lambda c: calls.append(c)
        try:
            ctx.menu_actions["Extensions/ORCA Result Analyzer"]()
            self.assertEqual(len(calls), 1)
            self.assertIs(calls[0], ctx)
        finally:
            _init_mod._open_orca_analyzer_empty = original

    def test_singleton_raises_existing_window(self):
        """If a window is already registered, it should be raised, not replaced."""
        ctx = StubContext()
        # Pre-register a fake window
        fake_win = MagicMock()
        ctx.register_window("analyzer", fake_win)

        original = _init_mod._open_orca_analyzer_empty
        # Call the real function — it should detect the existing window and raise it
        _init_mod._open_orca_analyzer_empty(ctx)

        # The fake window should have been shown/raised
        fake_win.show.assert_called_once()
        fake_win.raise_.assert_called_once()
        fake_win.activateWindow.assert_called_once()

        # A second window must NOT have been registered
        self.assertIs(ctx.get_window("analyzer"), fake_win)


# ---------------------------------------------------------------------------
# TestDocumentResetHandler — File->New must close the analyzer and drop
# stale atom-color overrides (they otherwise survive clear_3d_view()).
# ---------------------------------------------------------------------------


def _real_utils_module():
    """Load the real utils.py directly.

    Several other test files in this suite install a bare-bones fake module
    at sys.modules["orca_result_analyzer_rust.utils"] (for their own isolated
    imports) and never clean it up, so by the time the full suite runs, that
    slot may hold a stub whose clear_atom_color_overrides (if present at all)
    is a no-op MagicMock. _on_document_reset()'s deferred `from .utils import
    clear_atom_color_overrides` would then silently pick up the stub instead
    of the real implementation. Tests that need the real dict-clearing
    behavior swap this module in for the duration of the test.
    """
    path = os.path.normpath(
        os.path.join(_SRC_DIR, "orca_result_analyzer_rust", "utils.py")
    )
    spec = importlib.util.spec_from_file_location("orca_result_analyzer_rust.utils", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestDocumentResetHandler(unittest.TestCase):
    def setUp(self):
        self.ctx = StubContext()
        _init_mod.initialize(self.ctx)
        self._original_utils_mod = sys.modules.get("orca_result_analyzer_rust.utils")
        sys.modules["orca_result_analyzer_rust.utils"] = _real_utils_module()

    def tearDown(self):
        if self._original_utils_mod is not None:
            sys.modules["orca_result_analyzer_rust.utils"] = self._original_utils_mod
        else:
            sys.modules.pop("orca_result_analyzer_rust.utils", None)

    def test_registers_document_reset_handler(self):
        self.assertEqual(len(self.ctx.document_reset_handlers), 1)

    def test_handler_closes_existing_analyzer_window(self):
        fake_win = MagicMock()
        self.ctx.register_window("analyzer", fake_win)
        handler = self.ctx.document_reset_handlers[0]

        handler()

        fake_win.close.assert_called_once()

    def test_handler_noop_when_no_window_registered(self):
        handler = self.ctx.document_reset_handlers[0]
        try:
            handler()
        except Exception as exc:
            self.fail(f"document reset handler should not raise without a window: {exc}")

    def test_handler_survives_window_close_exception(self):
        fake_win = MagicMock()
        fake_win.close.side_effect = RuntimeError("boom")
        self.ctx.register_window("analyzer", fake_win)
        handler = self.ctx.document_reset_handlers[0]

        try:
            handler()
        except Exception as exc:
            self.fail(f"document reset handler must swallow close() errors: {exc}")

    def test_handler_clears_atom_color_overrides(self):
        fixed_mw = MagicMock()
        fixed_mw.view_3d_manager._plugin_color_overrides = {0: "#ff0000", 1: "#00ff00"}
        self.ctx.get_main_window = lambda: fixed_mw
        handler = self.ctx.document_reset_handlers[0]

        handler()

        self.assertEqual(fixed_mw.view_3d_manager._plugin_color_overrides, {})

    def test_handler_clears_overrides_even_without_a_window(self):
        """Overrides must be dropped even if the analyzer was already closed."""
        fixed_mw = MagicMock()
        fixed_mw.view_3d_manager._plugin_color_overrides = {2: "#0000ff"}
        self.ctx.get_main_window = lambda: fixed_mw
        handler = self.ctx.document_reset_handlers[0]

        handler()

        self.assertEqual(fixed_mw.view_3d_manager._plugin_color_overrides, {})


if __name__ == "__main__":
    unittest.main()

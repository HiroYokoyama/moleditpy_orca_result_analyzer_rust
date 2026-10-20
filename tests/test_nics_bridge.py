"""
tests/test_nics_bridge.py
Unit tests for orca_result_analyzer_rust/nics_bridge.py, the optional hand-off to
the separately installed ORCA NICS Analyzer plugin.

The module is loaded straight from its file: it is pure Python and imports Qt
only inside the menu-action fallback, which is exercised with a scoped stub.
"""

import os
import sys
import types
import importlib.util
import unittest
from unittest.mock import MagicMock, patch

_SRC = os.path.normpath(
    os.path.join(
        os.path.dirname(__file__), "..", "orca_result_analyzer_rust", "nics_bridge.py"
    )
)


def _load():
    spec = importlib.util.spec_from_file_location("orca_nics_bridge_mod", _SRC)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["orca_nics_bridge_mod"] = mod
    spec.loader.exec_module(mod)
    return mod


bridge = _load()


def _fake_nics_module(name="orca_nics_analyzer_fake"):
    mod = types.ModuleType(name)
    mod.PLUGIN_NAME = "ORCA NICS Analyzer"
    mod.PLUGIN_VERSION = "0.4.3"
    mod._context = MagicMock()
    mod._open_file = MagicMock()
    return mod


def _host(plugins=None):
    """Main window whose plugin manager lists *plugins* (or has none)."""
    mw = MagicMock()
    if plugins is None:
        mw.plugin_manager = None
    else:
        mw.plugin_manager.plugins = plugins
    return mw


class TestFindNicsModule(unittest.TestCase):
    def test_no_plugin_manager_means_not_installed(self):
        self.assertIsNone(bridge.find_nics_module(_host()))
        self.assertFalse(bridge.nics_analyzer_available(_host()))

    def test_none_main_window_is_tolerated(self):
        self.assertIsNone(bridge.find_nics_module(None))

    def test_found_by_plugin_name(self):
        mod = _fake_nics_module()
        host = _host([{"name": "Whatever", "module": mod}])
        self.assertIs(bridge.find_nics_module(host), mod)
        self.assertTrue(bridge.nics_analyzer_available(host))

    def test_found_by_registry_name(self):
        mod = types.ModuleType("plain")
        host = _host([{"name": "orca_nics_analyzer", "module": mod}])
        self.assertIs(bridge.find_nics_module(host), mod)

    def test_found_by_module_name(self):
        mod = types.ModuleType("plugins.orca_nics_analyzer")
        host = _host([{"name": "Something", "module": mod}])
        self.assertIs(bridge.find_nics_module(host), mod)

    def test_other_plugins_are_ignored(self):
        other = types.ModuleType("cif_viewer")
        other.PLUGIN_NAME = "CIF Viewer"
        host = _host([{"name": "CIF Viewer", "module": other}, "not-a-dict"])
        self.assertIsNone(bridge.find_nics_module(host))

    def test_falls_back_to_sys_modules_when_the_registry_is_unusable(self):
        mod = _fake_nics_module("orca_nics_analyzer_in_sys_modules")
        host = MagicMock()
        host.plugin_manager.plugins = 5  # not iterable as a plugin list
        with patch.dict(sys.modules, {mod.__name__: mod}):
            self.assertIs(bridge.find_nics_module(host), mod)

    def test_an_uninitialised_module_is_not_accepted_from_sys_modules(self):
        mod = _fake_nics_module("orca_nics_analyzer_uninitialised")
        mod._context = None
        with patch.dict(sys.modules, {mod.__name__: mod}):
            self.assertIsNone(bridge.find_nics_module(_host()))


class TestOpenNicsAnalyzer(unittest.TestCase):
    def test_missing_plugin_reports_how_to_install_it(self):
        ok, message = bridge.open_nics_analyzer(_host(), "job.out")
        self.assertFalse(ok)
        self.assertIn("not installed", message)

    def test_the_current_file_is_handed_over(self):
        mod = _fake_nics_module()
        host = _host([{"name": "nics", "module": mod}])
        ok, message = bridge.open_nics_analyzer(host, "job.out")
        self.assertTrue(ok)
        self.assertEqual(message, "")
        mod._open_file.assert_called_once_with("job.out", mod._context)

    def test_a_public_opener_wins_over_the_private_one(self):
        mod = _fake_nics_module()
        mod.open_file = MagicMock()
        host = _host([{"name": "nics", "module": mod}])
        bridge.open_nics_analyzer(host, "job.out")
        mod.open_file.assert_called_once_with("job.out", mod._context)
        mod._open_file.assert_not_called()

    def test_a_failing_hand_off_is_reported_not_raised(self):
        mod = _fake_nics_module()
        mod._open_file.side_effect = ValueError("bad ghost block")
        host = _host([{"name": "nics", "module": mod}])
        ok, message = bridge.open_nics_analyzer(host, "job.out")
        self.assertFalse(ok)
        self.assertIn("bad ghost block", message)

    def test_without_a_file_the_plugins_own_menu_entry_is_fired(self):
        mod = _fake_nics_module()
        host = _host([{"name": "nics", "module": mod}])
        action = MagicMock()
        action.text.return_value = "ORCA NICS Analyzer..."
        other = MagicMock()
        other.text.return_value = "Plugin Manager"
        host.findChildren.return_value = [other, action]
        with self._qtgui():
            ok, _ = bridge.open_nics_analyzer(host, None)
        self.assertTrue(ok)
        action.trigger.assert_called_once_with()
        other.trigger.assert_not_called()
        mod._open_file.assert_not_called()

    def test_the_analyzers_own_nics_entry_is_never_fired(self):
        # findChildren is recursive and reaches the analyzer's own
        # "NICS Analysis..." action, which must not re-trigger itself.
        mod = _fake_nics_module()
        host = _host([{"name": "nics", "module": mod}])
        own = MagicMock()
        own.text.return_value = "NICS Analysis..."
        plugin = MagicMock()
        plugin.text.return_value = "ORCA NICS Analyzer..."
        host.findChildren.return_value = [own, plugin]
        with self._qtgui():
            ok, _ = bridge.open_nics_analyzer(host, None, exclude_action=own)
        self.assertTrue(ok)
        own.trigger.assert_not_called()
        plugin.trigger.assert_called_once_with()

    def test_no_menu_entry_and_no_file_reports_failure(self):
        mod = _fake_nics_module()
        host = _host([{"name": "nics", "module": mod}])
        host.findChildren.return_value = []
        with self._qtgui():
            ok, message = bridge.open_nics_analyzer(host, "")
        self.assertFalse(ok)
        self.assertIn("could not be opened", message)

    def _qtgui(self):
        """Scoped PySide6.QtGui so the lazy QAction import resolves."""
        pkg = types.ModuleType("PySide6")
        pkg.__path__ = []
        gui = types.ModuleType("PySide6.QtGui")
        gui.QAction = type("QAction", (), {})
        pkg.QtGui = gui
        return patch.dict(sys.modules, {"PySide6": pkg, "PySide6.QtGui": gui})


if __name__ == "__main__":
    unittest.main()

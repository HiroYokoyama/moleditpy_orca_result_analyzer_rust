"""
tests/test_gui_nics_menu.py
Covers the Analysis > NICS Analysis entry, which is only meaningful while the
separate ORCA NICS Analyzer plugin is installed: it must hide itself when that
plugin is absent and hand the currently loaded file over when it is present.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(__file__))
import gui_harness  # noqa: E402

G = gui_harness.load_isolated("gui")


class _NicsCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.path = os.path.join(self._tmp.name, "job.out")
        with open(self.path, "w", encoding="utf-8") as fh:
            fh.write("ORCA output\n")

        self.parser = MagicMock()
        self.parser.data = {}
        self.mw = MagicMock()
        self.context = MagicMock()
        self.context.get_main_window.return_value = self.mw
        with patch.object(G, "nics_analyzer_available", return_value=True):
            self.dlg = G.OrcaResultAnalyzerDialog(
                self.mw, self.parser, self.path, context=self.context
            )


class TestVisibility(_NicsCase):
    def test_the_entry_is_hidden_while_the_plugin_is_absent(self):
        self.dlg.nics_action = MagicMock()
        with patch.object(G, "nics_analyzer_available", return_value=False) as avail:
            self.dlg._refresh_nics_action()
        avail.assert_called_once_with(self.mw)
        self.dlg.nics_action.setVisible.assert_called_once_with(False)

    def test_the_entry_shows_once_the_plugin_is_installed(self):
        self.dlg.nics_action = MagicMock()
        with patch.object(G, "nics_analyzer_available", return_value=True):
            self.dlg._refresh_nics_action()
        self.dlg.nics_action.setVisible.assert_called_once_with(True)

    def test_visibility_is_checked_while_building_the_menu(self):
        with patch.object(G, "nics_analyzer_available", return_value=False) as avail:
            G.OrcaResultAnalyzerDialog(
                self.mw, self.parser, self.path, context=self.context
            )
        avail.assert_called()

    def test_a_missing_action_is_not_an_error(self):
        self.dlg.nics_action = None
        with patch.object(G, "nics_analyzer_available") as avail:
            self.dlg._refresh_nics_action()
        avail.assert_not_called()

    def test_the_host_falls_back_to_the_parent_window(self):
        self.dlg.context = None
        self.assertIs(self.dlg._nics_host(), self.mw)

    def test_a_broken_host_lookup_falls_back_to_the_parent_window(self):
        self.context.get_main_window.side_effect = RuntimeError("wrapped C++ deleted")
        self.assertIs(self.dlg._nics_host(), self.mw)


class TestLaunch(_NicsCase):
    def test_the_loaded_file_is_handed_to_the_nics_plugin(self):
        with patch.object(G, "open_nics_analyzer", return_value=(True, "")) as opener:
            with patch.object(G.QMessageBox, "information") as info:
                self.dlg.show_nics_analysis()
        opener.assert_called_once_with(
            self.mw, self.path, exclude_action=self.dlg.nics_action
        )
        info.assert_not_called()

    def test_a_failure_is_explained_to_the_user(self):
        with patch.object(G, "open_nics_analyzer", return_value=(False, "no plugin")):
            with patch.object(G.QMessageBox, "information") as info:
                self.dlg.show_nics_analysis()
        info.assert_called_once()
        self.assertIn("no plugin", info.call_args[0][2])


if __name__ == "__main__":
    unittest.main()

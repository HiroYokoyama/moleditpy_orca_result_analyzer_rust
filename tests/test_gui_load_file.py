"""
tests/test_gui_load_file.py
Covers OrcaResultAnalyzerDialog.load_file: encoding fallback, the NEB
trajectory auto-load that upgrades a path summary to real frames, and the
state reset that has to happen when a new result replaces the previous one.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(__file__))
import gui_harness  # noqa: E402

G = gui_harness.load_isolated("gui")

OUT_TEXT = """
                       * O   R   C   A *

Program Version 5.0.4

CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.000000    0.000000    0.000000
  H      0.960000    0.000000    0.000000

FINAL SINGLE POINT ENERGY      -76.400000000000

                             ****ORCA TERMINATED NORMALLY****
"""


def _xyz(frames):
    out = []
    for energy, x in frames:
        out.append("2")
        out.append(f"Coordinates from ORCA-job E {energy:.8f}")
        out.append("O    0.000000    0.000000    0.000000")
        out.append(f"H    {x:.6f}    0.000000    0.000000")
    return "\n".join(out) + "\n"


class _LoadCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.addCleanup(self._tmp.cleanup)

        self.first = self._write("first.out", OUT_TEXT)

        self.context = MagicMock()
        self.mw = MagicMock()
        self.parser = MagicMock()
        self.parser.filename = self.first
        self.parser.data = {}

        self.dlg = G.OrcaResultAnalyzerDialog(
            self.mw, self.parser, self.first, context=self.context
        )
        # Drawing needs a live 3D viewer; covered separately.
        patcher = patch.object(self.dlg, "load_structure_3d")
        self.draw = patcher.start()
        self.addCleanup(patcher.stop)

    def _write(self, name, text, encoding="utf-8"):
        path = os.path.join(self.tmp, name)
        with open(path, "w", encoding=encoding) as fh:
            fh.write(text)
        return path


class TestLoadFile(_LoadCase):
    def test_the_file_becomes_the_current_result(self):
        path = self._write("second.out", OUT_TEXT)
        self.dlg.load_file(path)
        self.assertEqual(self.dlg.file_path, path)

    def test_the_output_is_parsed(self):
        path = self._write("second.out", OUT_TEXT)
        self.dlg.load_file(path)
        self.assertEqual(self.dlg.parser.data.get("version"), "5.0.4")

    def test_the_geometry_is_parsed(self):
        path = self._write("second.out", OUT_TEXT)
        self.dlg.load_file(path)
        self.assertEqual(len(self.dlg.parser.data.get("coords", [])), 2)

    def test_the_structure_is_reframed_for_a_new_result(self):
        path = self._write("second.out", OUT_TEXT)
        self.dlg.load_file(path)
        self.draw.assert_called_once_with(fit_camera=True)

    def test_open_analysis_windows_are_closed_first(self):
        stale = MagicMock()
        self.dlg.mo_dlg = stale
        self.dlg.load_file(self._write("second.out", OUT_TEXT))
        stale.close.assert_called_once()

    def test_atom_colours_from_the_previous_result_are_cleared(self):
        with patch.object(G, "clear_atom_color_overrides") as clear:
            self.dlg.load_file(self._write("second.out", OUT_TEXT))
        clear.assert_called_once_with(self.dlg.mw)

    def test_the_main_window_title_is_synced(self):
        path = self._write("second.out", OUT_TEXT)
        self.dlg.load_file(path)
        self.assertEqual(self.mw.init_manager.current_file_path, path)

    def test_success_is_reported(self):
        self.dlg.load_file(self._write("second.out", OUT_TEXT))
        self.assertTrue(
            any(
                "Successfully loaded" in str(c.args[0])
                for c in self.context.show_status_message.call_args_list
            )
        )

    def test_a_utf16_output_is_read(self):
        path = self._write("utf16.out", OUT_TEXT, encoding="utf-16")
        self.dlg.load_file(path)
        self.assertEqual(self.dlg.parser.data.get("version"), "5.0.4")

    def test_undecodable_bytes_do_not_stop_the_load(self):
        path = os.path.join(self.tmp, "binary.out")
        with open(path, "wb") as fh:
            fh.write(OUT_TEXT.encode("utf-8") + b"\xff\xfe\x00rubbish")
        self.dlg.load_file(path)
        self.assertEqual(self.dlg.file_path, path)

    def test_a_missing_file_is_reported(self):
        with patch.object(G.QMessageBox, "critical") as crit:
            self.dlg.load_file(os.path.join(self.tmp, "nope.out"))
        crit.assert_called_once()

    def test_a_failed_load_keeps_the_previous_result(self):
        before = self.dlg.file_path
        with patch.object(G.QMessageBox, "critical"):
            self.dlg.load_file(os.path.join(self.tmp, "nope.out"))
        self.assertEqual(self.dlg.file_path, before)


class TestNebTrajectoryAutoload(_LoadCase):
    def test_a_sibling_mep_trajectory_is_picked_up(self):
        path = self._write("neb.out", OUT_TEXT)
        self._write("neb_MEP_trj.xyz", _xyz([(-100.0, 0.96), (-100.5, 0.97)]))
        self.dlg.load_file(path)
        self.assertEqual(len(self.dlg.parser.data.get("scan_steps", [])), 2)

    def test_the_trajectory_replaces_the_summary_geometry(self):
        path = self._write("neb.out", OUT_TEXT)
        self._write("neb_MEP_trj.xyz", _xyz([(-100.0, 0.96), (-100.5, 1.40)]))
        self.dlg.load_file(path)
        # the final frame becomes the displayed structure
        self.assertAlmostEqual(self.dlg.parser.data["coords"][1][0], 1.40, places=5)

    def test_loading_a_trajectory_is_reported(self):
        path = self._write("neb.out", OUT_TEXT)
        self._write("neb_MEP_trj.xyz", _xyz([(-100.0, 0.96)]))
        self.dlg.load_file(path)
        self.assertTrue(
            any(
                "NEB Trajectory" in str(c.args[0])
                for c in self.context.show_status_message.call_args_list
            )
        )

    def test_no_sibling_trajectory_leaves_the_summary_alone(self):
        path = self._write("plain.out", OUT_TEXT)
        self.dlg.load_file(path)
        self.assertEqual(self.dlg.parser.data.get("scan_steps", []), [])

    def test_an_empty_trajectory_is_ignored(self):
        path = self._write("neb.out", OUT_TEXT)
        self._write("neb_MEP_trj.xyz", "")
        self.dlg.load_file(path)
        self.assertEqual(self.dlg.parser.data.get("scan_steps", []), [])

    def test_an_unreadable_trajectory_does_not_stop_the_load(self):
        path = self._write("neb.out", OUT_TEXT)
        self._write("neb_MEP_trj.xyz", "not an xyz file at all\n")
        self.dlg.load_file(path)
        self.assertEqual(self.dlg.file_path, path)


if __name__ == "__main__":
    unittest.main()

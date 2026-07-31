"""
tests/test_mo_selection.py
Covers MODialog's selection handling: auto-loading an already-generated cube
when a row is chosen, the coefficient check that gates the Visualize button
(and the "add these ORCA keywords" hint shown when they are missing), and the
batch generation queue.

The worker that actually evaluates the basis set is stubbed out — these tests
are about the queueing and gating around it.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(__file__))
import gui_harness  # noqa: E402

M = gui_harness.load_isolated("mo_analysis")


def _mos():
    return {
        str(i): {"id": i, "energy": e, "occ": occ, "spin": "restricted"}
        for i, (e, occ) in enumerate([(-1.0, 2.0), (-0.5, 2.0), (0.2, 0.0)])
    }


class _MOCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.addCleanup(self._tmp.cleanup)

        saved = M.__file__
        M.__file__ = os.path.join(self.tmp, "mo_analysis.py")
        self.addCleanup(lambda: setattr(M, "__file__", saved))

        self.out = os.path.join(self.tmp, "job.out")
        self.host = MagicMock()
        self.host.parser.filename = self.out
        self.host.parser.data = {"mo_coeffs": {"1": {"coeffs": [0.1, 0.2]}}}
        self.host.file_path = self.out

        self.dlg = M.MODialog(self.host, _mos(), result_dir=self.tmp)
        self.dlg.parent_dlg = self.host
        self.dlg.mw = MagicMock()

    def _row(self, mo_id):
        """The tree row whose MO column reads `mo_id`."""
        tree = self.dlg.tree
        return next(
            tree.topLevelItem(i)
            for i in range(tree.topLevelItemCount())
            if tree.topLevelItem(i).text(0) == str(mo_id)
        )

    def _select(self, *mo_ids):
        rows = [self._row(i) for i in mo_ids]
        self.dlg.tree.selectedItems = lambda: rows
        return rows

    def _cube_for(self, mo_id):
        cube_dir = os.path.join(self.tmp, "job_cubes")
        os.makedirs(cube_dir, exist_ok=True)
        path = os.path.join(cube_dir, f"job_MO_{mo_id}.cube")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("cube")
        return path


# ---------------------------------------------------------------------------
# Choosing a row
# ---------------------------------------------------------------------------


class TestRowSelected(_MOCase):
    def test_an_already_generated_cube_is_shown(self):
        path = self._cube_for(1)
        with patch.object(self.dlg, "show_cube") as show:
            self.dlg.on_item_changed(self._row(1), None)
        show.assert_called_once_with(path)

    def test_a_row_without_a_cube_shows_nothing(self):
        with patch.object(self.dlg, "show_cube") as show:
            self.dlg.on_item_changed(self._row(1), None)
        show.assert_not_called()

    def test_clearing_the_row_selection_is_ignored(self):
        with patch.object(self.dlg, "show_cube") as show:
            self.dlg.on_item_changed(None, None)
        show.assert_not_called()


class TestVisualiseGating(_MOCase):
    def test_an_orbital_with_coefficients_can_be_visualised(self):
        self._select(1)
        self.dlg.on_selection_changed()
        self.assertTrue(self.dlg.btn_vis.isEnabled())

    def test_an_orbital_without_coefficients_cannot(self):
        self._select(2)
        self.dlg.on_selection_changed()
        self.assertFalse(self.dlg.btn_vis.isEnabled())

    def test_a_missing_orbital_explains_itself(self):
        self._select(2)
        self.dlg.on_selection_changed()
        self.assertIn("No coefficients", self.dlg.btn_vis.toolTip())

    def test_missing_coefficients_surface_the_required_keywords(self):
        self._select(2)
        self.dlg.on_selection_changed()
        self.assertTrue(self.dlg.lbl_warning.isVisible())
        self.assertIn("Print[P_Basis]", self.dlg.lbl_warning.text())
        self.assertTrue(self.dlg.btn_copy_input.isVisible())

    def test_the_hint_is_hidden_once_coefficients_are_present(self):
        self._select(2)
        self.dlg.on_selection_changed()
        self._select(1)
        self.dlg.on_selection_changed()
        self.assertFalse(self.dlg.lbl_warning.isVisible())
        self.assertEqual(self.dlg.btn_vis.toolTip(), "")

    def test_no_selection_disables_visualisation(self):
        self.dlg.tree.selectedItems = lambda: []
        self.dlg.on_selection_changed()
        self.assertFalse(self.dlg.btn_vis.isEnabled())

    def test_no_selection_does_not_show_the_hint(self):
        self.dlg.tree.selectedItems = lambda: []
        self.dlg.on_selection_changed()
        self.assertFalse(self.dlg.lbl_warning.isVisible())

    def test_a_result_without_coefficients_at_all_is_handled(self):
        self.host.parser.data = {}
        self._select(1)
        self.dlg.on_selection_changed()
        self.assertFalse(self.dlg.btn_vis.isEnabled())


# ---------------------------------------------------------------------------
# Batch generation queue
# ---------------------------------------------------------------------------


class TestGenerationQueue(_MOCase):
    def test_visualising_queues_the_selected_orbitals(self):
        self._select(0, 1, 2)
        with patch.object(self.dlg, "_generate_single_mo"):
            self.dlg.visualize_selected_mos()
        # the first is popped and handed to the worker straight away
        self.assertEqual(self.dlg.generation_queue, ["1", "2"])

    def test_the_first_orbital_starts_generating_immediately(self):
        self._select(0, 1)
        with patch.object(self.dlg, "_generate_single_mo") as gen:
            self.dlg.visualize_selected_mos()
        gen.assert_called_once_with("0")

    def test_visualising_nothing_queues_nothing(self):
        self.dlg.tree.selectedItems = lambda: []
        with patch.object(self.dlg, "_generate_single_mo") as gen:
            self.dlg.visualize_selected_mos()
        gen.assert_not_called()
        self.assertEqual(self.dlg.generation_queue, [])

    def test_the_queue_is_replaced_not_appended(self):
        self.dlg.generation_queue = ["stale"]
        self._select(1)
        with patch.object(self.dlg, "_generate_single_mo"):
            self.dlg.visualize_selected_mos()
        self.assertEqual(self.dlg.generation_queue, [])

    def test_the_queue_is_drained_one_at_a_time(self):
        self.dlg.generation_queue = ["0", "1", "2"]
        with patch.object(self.dlg, "_generate_single_mo") as gen:
            self.dlg.process_generation_queue()
        gen.assert_called_once_with("0")
        self.assertEqual(self.dlg.generation_queue, ["1", "2"])

    def test_draining_an_empty_queue_closes_the_progress_dialog(self):
        progress = MagicMock()
        self.dlg.progress_dialog = progress
        self.dlg.generation_queue = []
        self.dlg.process_generation_queue()
        progress.close.assert_called_once()
        self.assertIsNone(self.dlg.progress_dialog)

    def test_draining_an_empty_queue_starts_no_work(self):
        self.dlg.generation_queue = []
        self.dlg.progress_dialog = None
        with patch.object(self.dlg, "_generate_single_mo") as gen:
            self.dlg.process_generation_queue()
        gen.assert_not_called()


if __name__ == "__main__":
    unittest.main()

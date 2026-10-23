"""Scientific data association and empty-view regressions found in the audit."""

import math
from types import SimpleNamespace
from unittest.mock import Mock, patch

import pytest

from ._parser_loader import load_standalone_parser

OrcaParser = load_standalone_parser("audit_parser").OrcaParser


def parse(text, method):
    parser = OrcaParser()
    parser.lines = text.splitlines()
    getattr(parser, method)()
    return parser.data


def test_mo_coefficient_rows_before_a_missing_header_do_not_crash():
    data = parse(
        "MOLECULAR ORBITALS\n------------------\n0C 1s 0.123\n", "parse_mo_coeffs"
    )
    assert data["mo_coeffs"] == {}


@pytest.mark.parametrize("number", ["-2.1e-4", "-2.1D-4", "-.21E-3", "0"])
def test_energy_components_preserve_the_exponent(number):
    data = parse(f"E(CORR)(corrected) ... {number} Eh", "parse_energy_components")
    assert data["energy_components"][0]["value"] == float(number.replace("D", "E"))


def test_a_missing_energy_value_is_not_read_from_the_label():
    data = parse("E(0) ... ***\nT1 diagnostic ... ***", "parse_energy_components")
    assert data["energy_components"] == []


@pytest.mark.parametrize(
    "title,row,expected",
    [
        ("MULLIKEN ATOMIC CHARGES", "0C -0.25", {"charge": -0.25}),
        ("HIRSHFELD ANALYSIS", "0C -0.25 0.75", {"charge": -0.25, "spin": 0.75}),
        (
            "MBIS ANALYSIS",
            "0C -0.25 6.25 0.75",
            {"charge": -0.25, "population": 6.25, "spin": 0.75},
        ),
    ],
)
def test_merged_atom_labels_do_not_shift_charge_columns(title, row, expected):
    data = parse(f"{title}\n-------------------\n\n{row}\n", "parse_charges")
    atoms = next(iter(data["charges"].values()))
    assert atoms[0]["atom_idx"] == 0
    assert atoms[0]["atom_sym"] == "C"
    for key, value in expected.items():
        assert atoms[0][key] == value


def test_optimization_cycle_never_uses_the_next_cycles_geometry():
    data = parse(
        "OPTIMIZATION CYCLE 1\nFINAL SINGLE POINT ENERGY -1.0\n"
        "OPTIMIZATION CYCLE 2\nCARTESIAN COORDINATES (ANGSTROEM)\n---\n"
        "H 1 2 3\n\nFINAL SINGLE POINT ENERGY -2.0\n",
        "parse_trajectory",
    )
    cycles = [step for step in data["scan_steps"] if step["type"] == "opt_cycle"]
    assert [step["step"] for step in cycles] == [2]
    assert cycles[0]["energy"] == -2.0


@pytest.mark.parametrize("first", [0, 1])
def test_gradient_atom_numbering_is_consistent_for_the_whole_block(first):
    rows = "\n".join(f"{first + i} H : {i + 1}.0 0.0 0.0" for i in range(3))
    data = parse("CARTESIAN GRADIENT\n---\n" + rows + "\n", "parse_gradients")
    assert [row["atom_idx"] for row in data["gradients"]] == [0, 1, 2]
    assert [row["vector"][0] for row in data["gradients"]] == [1.0, 2.0, 3.0]


@pytest.mark.parametrize("value", ["nan", "inf", "1e999"])
def test_nonfinite_coordinates_are_not_sent_to_the_renderer(value):
    assert OrcaParser._parse_xyz_row(f"H {value} 0 0") is None
    data = OrcaParser().parse_xyz_content(f"1\nenergy -1\nH {value} 0 0\n")
    assert all(math.isfinite(x) for step in data for xyz in step["coords"] for x in xyz)


def test_intensities_follow_the_printed_mode_id_after_a_missing_frequency():
    data = parse(
        "VIBRATIONAL FREQUENCIES\n---\n0: 0.0 cm**-1\n"
        "1: *** cm**-1\n2: 1200.0 cm**-1\n---\nIR SPECTRUM\n"
        "---\nMode freq Int\nunits\n---\n0: 0.0 0.0\n2: 1200.0 35.0\n",
        "parse_frequencies",
    )
    assert len(data["frequencies"]) == 2
    assert data["frequencies"][1]["freq"] == 1200.0
    assert data["frequencies"][1]["ir"] == 35.0


def test_energy_arrows_do_not_cluster_distant_levels_in_reverse_order():
    from .test_energy_diag import E

    assert E.calculate_arrow_shifts([(0, 100.0), (1, 0.0)], lambda value: value) == {
        0: 0,
        1: 0,
    }


def test_negative_sticks_are_inside_the_automatic_y_range():
    from .test_spectrum_widget import _widget

    widget = _widget([{"energy_nm": 400.0, "osc": -5.0}])
    widget.use_dual_axis = False
    widget.show_gaussian = False
    widget.plot_spectrum()
    low, high = widget.canvas.axes.get_ylim()
    assert low < -5.0 < high


@pytest.mark.parametrize("shifts,axes", [([], Mock()), ([1.0], None)])
def test_clicking_an_empty_nmr_view_or_outside_axes_is_harmless(shifts, axes):
    from . import gui_harness

    module = gui_harness.load_isolated("nmr_plot")
    owner = SimpleNamespace(current_shifts=shifts)
    module._NMRPlotMixin.on_peak_click(
        owner, SimpleNamespace(button=1, xdata=1.0, inaxes=axes)
    )


def test_large_cached_mo_batches_do_not_recurse():
    from . import gui_harness

    module = gui_harness.load_isolated("mo_analysis")
    consumed = []
    owner = SimpleNamespace(
        generation_queue=list(range(1500)),
        progress_dialog=None,
        generation_done_cb=Mock(),
        _is_closing=False,
    )
    owner.process_generation_queue = lambda: module.MODialog.process_generation_queue(
        owner
    )

    def cached(key):
        consumed.append(key)
        owner.process_generation_queue()

    owner._generate_single_mo = cached
    completed = owner.generation_done_cb
    owner.process_generation_queue()
    assert consumed == list(range(1500))
    completed.assert_called_once_with()
    assert owner.generation_done_cb is None


def test_status_message_survives_a_deleted_qt_parent():
    from . import gui_harness

    module = gui_harness.load_isolated("utils")

    class Deleted:
        def parent(self):
            raise RuntimeError("Internal C++ object already deleted")

    assert module.notify(Deleted(), "Analysis completed") is False


def test_detailed_excited_states_get_wavelengths_without_a_summary_table():
    data = parse(
        "TD-DFT EXCITED STATES\nSTATE 1: E= 0.1250 au 3.4000 eV 27422.85 cm**-1\n",
        "parse_tddft",
    )
    assert data["tddft"][0]["energy_ev"] == 3.4
    assert data["tddft"][0]["energy_nm"] == pytest.approx(1239.84193 / 3.4)


def test_comparison_callback_cannot_redraw_after_the_window_closes():
    from . import gui_harness

    module = gui_harness.load_isolated("mo_compare")
    module.MOCompareDialog.render_all(SimpleNamespace(_is_closing=True))


def test_empty_charge_table_cannot_inherit_the_next_population_scheme():
    data = parse(
        "MULLIKEN ATOMIC CHARGES\n---\n\nLOEWDIN ATOMIC CHARGES\n---\n0 C : -0.4\n\n",
        "parse_charges",
    )
    assert "Mulliken" not in data["charges"]
    assert data["charges"]["Loewdin"][0]["charge"] == -0.4


@pytest.mark.parametrize(
    "name,cls",
    [("freq_analysis", "FreqSpectrumWindow"), ("tddft_analysis", "TDDFTDialog")],
)
def test_spectrum_image_export_reports_a_write_failure(name, cls):
    from . import gui_harness

    module = gui_harness.load_isolated(name)
    owner = SimpleNamespace(
        spectrum=Mock(),
        parent=lambda: SimpleNamespace(file_path="job.out"),
        freq_dialog=SimpleNamespace(
            mw=SimpleNamespace(
                init_manager=SimpleNamespace(current_file_path="job.out")
            )
        ),
    )
    owner.spectrum.save_png.side_effect = PermissionError("read only")
    with (
        patch.object(module, "QFileDialog") as chooser,
        patch.object(module, "QMessageBox") as message,
    ):
        chooser.getSaveFileName.return_value = ("graph.png", "")
        getattr(module, cls).save_png(owner)
    message.critical.assert_called_once()


def test_trajectory_export_restores_tooltip_after_a_write_failure():
    from . import gui_harness

    module = gui_harness.load_isolated("traj_analysis")
    owner = SimpleNamespace(annot=Mock(), canvas=Mock(), output_path="job.out")
    owner.annot.get_visible.return_value = True
    owner.canvas.fig.savefig.side_effect = PermissionError("read only")
    with (
        patch.object(module, "QFileDialog") as chooser,
        patch.object(module, "QMessageBox") as message,
    ):
        chooser.getSaveFileName.return_value = ("graph.png", "")
        module.TrajectoryResultDialog.save_graph(owner)
    message.critical.assert_called_once()
    assert owner.annot.set_visible.call_args.args == (True,)


def test_nmr_table_csv_preserves_quotes_and_line_breaks(tmp_path):
    import csv
    from . import gui_harness

    module = gui_harness.load_isolated("nmr_export")
    cell = 'quoted "reference",\nsecond line'
    table = Mock()
    table.columnCount.return_value = 1
    table.rowCount.return_value = 1
    table.horizontalHeaderItem.return_value.text.return_value = "Value"
    table.item.return_value.text.return_value = cell
    owner = SimpleNamespace(current_nucleus="H", file_path="job.out", table=table)
    path = tmp_path / "table.csv"
    with patch.object(module, "QFileDialog") as chooser:
        chooser.getSaveFileName.return_value = (str(path), "")
        module._NMRExportMixin.export_table_csv(owner)
    with path.open(newline="", encoding="utf-8") as stream:
        assert list(csv.reader(stream)) == [["Value"], [cell]]

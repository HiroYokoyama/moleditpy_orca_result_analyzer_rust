"""Regressions at the native parser / MO engine boundary."""

import importlib

import numpy as np
import pytest

from ._parser_loader import load_standalone_parser

_parser = load_standalone_parser("native_result_audit")
native = _parser._rs
mo_native = importlib.import_module(f"{_parser.__package__}.orca_mo_rs")


def shell(**changes):
    return {
        "type": 0,
        "center": [0.0, 0.0, 0.0],
        "exps": [1.0],
        "coeffs": [1.0],
        **changes,
    }


def geometry(x):
    return f"CARTESIAN COORDINATES (ANGSTROEM)\n---\nH {x} 0 0\n\n"


def test_native_empty_input_keeps_singlet_default():
    assert native.parse_all("")["mult"] == 1


@pytest.mark.parametrize("value", ["nan", "inf", "-inf", "1e999"])
def test_native_coordinate_paths_reject_nonfinite_values(value):
    assert native.parse_all(geometry(value))["coords"] == []
    steps = native.parse_xyz_content(f"1\nenergy -1\nH {value} 0 0\n")
    assert steps[0]["coords"] == []


def test_native_xyz_huge_excluded_frame_count_does_not_overflow():
    assert native.parse_xyz_content("18446744073709551615\nTS\n") == []


def test_native_xyz_rejects_zero_atom_counts():
    assert native.parse_xyz_content("0\nenergy -1\n") == []


def test_scan_summary_uses_last_geometry_within_its_own_point():
    data = native.parse_all(
        "RELAXED SURFACE SCAN STEP 1\n"
        + geometry(1)
        + "OPTIMIZATION CYCLE 1\n"
        + geometry(2)
        + "FINAL SINGLE POINT ENERGY -2\nRELAXED SURFACE SCAN STEP 2\n"
        + geometry(3)
        + "FINAL SINGLE POINT ENERGY -3\n"
    )
    summaries = [step for step in data["scan_steps"] if step["type"] == "scan_step"]
    assert [step["coords"][0][0] for step in summaries] == [2.0, 3.0]


def test_geometryless_final_point_does_not_borrow_next_scan_point():
    data = native.parse_all(
        "RELAXED SURFACE SCAN STEP 1\n"
        "FINAL ENERGY EVALUATION AT THE STATIONARY POINT\n"
        "RELAXED SURFACE SCAN STEP 2\n" + geometry(3) + "FINAL SINGLE POINT ENERGY -3\n"
    )
    assert not any(step["type"] == "opt_final" for step in data["scan_steps"])


def test_native_modes_use_printed_ids_and_require_complete_vectors():
    data = native.parse_all(
        geometry(0)
        + "VIBRATIONAL FREQUENCIES\n---\n0: 0 cm**-1\n1: *** cm**-1\n2: 1200 cm**-1\n---\n"
        "NORMAL MODES\n---\n\n\n\n\n\n0 2\n0 0.0 1.0\n1 0.0 2.0\n2 0.0 3.0\n"
    )
    assert data["frequencies"][1]["vector"] == [[1.0, 2.0, 3.0]]
    incomplete = native.parse_all(
        geometry(0) + "VIBRATIONAL FREQUENCIES\n---\n0: 1200 cm**-1\n---\n"
        "NORMAL MODES\n---\n\n\n\n\n\n0\n0 1.0\n1 ***\n2 3.0\n"
    )
    assert incomplete["frequencies"][0]["vector"] == []


@pytest.mark.parametrize("missing", ["type", "center", "exps", "coeffs"])
def test_missing_shell_keys_raise_python_keyerror(missing):
    data = shell()
    del data[missing]
    with pytest.raises(KeyError, match=missing):
        mo_native.BasisSetEngineRust([data])


@pytest.mark.parametrize(
    "changes",
    [
        {"center": []},
        {"center": [0.0, 0.0]},
        {"center": [0.0, float("nan"), 0.0]},
        {"exps": []},
        {"exps": [1.0, 2.0]},
        {"exps": [0.0]},
        {"exps": [-1.0]},
        {"exps": [float("inf")]},
        {"coeffs": [float("nan")]},
    ],
)
def test_invalid_shell_shapes_and_values_raise_python_valueerror(changes):
    with pytest.raises(ValueError):
        mo_native.BasisSetEngineRust([shell(**changes)])


@pytest.mark.parametrize(
    "grid,coefficients",
    [
        ([0.0, 0.0], [1.0]),
        ([0.0, float("nan"), 0.0], [1.0]),
        ([0.0, 0.0, 0.0], []),
        ([0.0, 0.0, 0.0], [1.0, 2.0]),
        ([0.0, 0.0, 0.0], [float("inf")]),
    ],
)
def test_invalid_native_grid_requests_raise_valueerror(grid, coefficients):
    engine = mo_native.BasisSetEngineRust([shell()])
    with pytest.raises(ValueError):
        engine.evaluate_mo_on_grid(grid, coefficients)


def test_python_bridge_refuses_misinterpreting_two_column_grids():
    from .test_mo_engine import _mod

    engine = _mod.BasisSetEngine([shell()])
    with pytest.raises(ValueError, match="shape"):
        engine.evaluate_mo_on_grid(0, np.zeros((3, 2)), [1.0])

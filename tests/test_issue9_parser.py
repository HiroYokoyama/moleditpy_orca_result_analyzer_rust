"""Regression coverage for issue #9 using real ORCA 5/6 outputs."""
from pathlib import Path
import pytest
from tests._parser_loader import load_standalone_parser

OrcaParser = load_standalone_parser().OrcaParser
SAMPLES = Path(__file__).parent / "sample_outputs"


def output(name):
    return (SAMPLES / name).read_text(encoding="utf-8")


@pytest.mark.parametrize("name", ["benzene-opt-ene.out", "benzene-opt-ene_5.out", "benzene-opt-eneQZ.out"])
def test_every_mo_has_all_basis_coefficients(name):
    parser = OrcaParser()
    parser.load_from_memory(output(name))
    size = sum(2 * shell["l"] + 1 for shell in parser.data["basis_set_shells"])
    assert size in (114, 522)
    assert len(parser.data["mo_coeffs"]) == size
    assert all(len(mo["coeffs"]) == size for mo in parser.data["mo_coeffs"].values())


@pytest.mark.parametrize("first,last", [
    ("benzene-opt-ene.out", "benzene-opt-eneQZ.out"),
    ("o2-opt_orb.out", "benzene-opt-ene.out"),
    ("benzene-opt-ene.out", "o2-opt_orb.out"),
])
def test_final_job_matches_loading_it_alone(first, last):
    combined, final = OrcaParser(), OrcaParser()
    combined.load_from_memory(output(first) + "\n" + output(last))
    final.load_from_memory(output(last))
    assert combined.data == final.data


def test_latest_equal_sized_basis_parameters():
    parser = OrcaParser()
    parser.data.update(atoms=["H"], coords=[[0, 0, 0]])
    def basis(exponent):
        return f"BASIS SET IN INPUT FORMAT\n--------------------\nNewGTO H\nS 1\n1 {exponent} 1.0\nend;\n"
    parser.lines = (basis(1.0) + basis(2.0)).splitlines()
    parser.parse_basis_set()
    assert parser.data["basis_set_shells"][0]["exps"] == [2.0]


def test_final_basis_and_coefficients_have_matching_size():
    parser = OrcaParser()
    parser.load_from_memory(output("benzene-opt-ene.out") + "\n" + output("benzene-opt-eneQZ.out"))
    size = sum(2 * shell["l"] + 1 for shell in parser.data["basis_set_shells"])
    assert size == 522
    assert all(len(mo["coeffs"]) == size for mo in parser.data["mo_coeffs"].values())


def test_compound_job_does_not_inherit_previous_results():
    parser = OrcaParser()
    parser.load_from_memory(output("o2-opt_orb.out") + "\nJOB NUMBER 2\n")
    assert parser.data["atoms"] == []
    assert parser.data["mo_coeffs"] == {}
    assert not parser.data.get("spin_s2")
    assert parser.data["version"] == "6.1.1"


def test_repeated_load_resets_all_previous_results():
    parser = OrcaParser()
    parser.load_from_memory(output("acetone-opt.out"))
    assert len(parser.data["scan_steps"]) == 6
    parser.parse_trajectory()
    assert len(parser.data["scan_steps"]) == 6
    parser.load_from_memory("Program Version 6.1.1\n")
    fresh = OrcaParser()
    fresh.load_from_memory("Program Version 6.1.1\n")
    assert parser.data == fresh.data


def test_incomplete_cycle_keeps_geometry_without_inventing_energy():
    text = output("acetone-opt.out")
    start = text.index("GEOMETRY OPTIMIZATION CYCLE   2")
    stop = text.index("CARTESIAN COORDINATES (A.U.)", start)
    parser = OrcaParser()
    parser.load_from_memory(text[:stop])
    step = parser.data["scan_steps"][-1]
    assert step["step"] == 2
    assert len(step["atoms"]) == 10
    assert len(step["coords"]) == 10
    assert step["energy"] is None
    assert parser.data["scan_steps"][0]["energy"] < -192

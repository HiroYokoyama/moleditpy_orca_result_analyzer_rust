"""Check shared Python code against the canonical checkout pinned in SYNC_SOURCE.json."""

import argparse
import ast
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def adapted(text):
    return text.replace("orca_result_analyzer", "orca_result_analyzer_rust")


def tree(text):
    return ast.parse(text)


def check(source):
    manifest = json.loads((ROOT / "SYNC_SOURCE.json").read_text(encoding="utf-8"))
    revision = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    errors = []
    if revision != manifest["commit"]:
        errors.append(f"canonical commit is {revision}, expected {manifest['commit']}")
    for name in manifest["shared_modules"] + manifest["shared_parser_mixins"]:
        expected = adapted(
            (source / "orca_result_analyzer" / name).read_text(encoding="utf-8")
        )
        if name == "__init__.py":
            expected = expected.replace(
                'PLUGIN_SUPPORTED_PYTHON_VERSION = ">=3.9, <3.15"',
                'PLUGIN_SUPPORTED_PYTHON_VERSION = ">=3.9, <3.14"',
            )
        actual = (ROOT / "orca_result_analyzer_rust" / name).read_text(encoding="utf-8")
        if actual != expected:
            errors.append(f"shared module differs: {name}")

    expected_engine = tree(
        adapted(
            (source / "orca_result_analyzer/mo_engine.py").read_text(encoding="utf-8")
        )
    )
    actual_engine = tree(
        (ROOT / "orca_result_analyzer_rust/mo_engine.py").read_text(encoding="utf-8")
    )
    for name in (
        "CubeWriter",
        "CalcWorker",
        "UnsupportedBasisError",
        "read_generation_settings",
    ):
        expected = next(
            n for n in expected_engine.body if getattr(n, "name", None) == name
        )
        actual = next(n for n in actual_engine.body if getattr(n, "name", None) == name)
        if ast.dump(actual) != ast.dump(expected):
            errors.append(f"shared MO API differs: {name}")

    for name in manifest.get("shared_images", []):
        upstream_image = source / "img" / name
        local_image = ROOT / "img" / name
        if (
            not local_image.is_file()
            or upstream_image.read_bytes() != local_image.read_bytes()
        ):
            errors.append(f"shared image differs: {name}")

    # The checker engine is owned by the plugin collection, independently of Python UI parity.
    checker_lock = json.loads(
        (ROOT / "scripts/api-checker.lock.json").read_text(encoding="utf-8")
    )
    checker = (
        (ROOT / "tests/plugin_api_checker.py").read_bytes().replace(b"\r\n", b"\n")
    )
    if hashlib.sha256(checker).hexdigest() != checker_lock["sha256_lf"]:
        errors.append("canonical API checker copy differs from its pinned hash")

    upstream_tests = source / "tests"
    for path in upstream_tests.rglob("*"):
        if not path.is_file() or "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        relative = path.relative_to(upstream_tests)
        dest = ROOT / "tests" / relative
        if str(relative) == "plugin_api_checker.py":
            continue
        if str(relative) == "test_api.py":
            # Retain consumer integration parity without reintroducing shared engine unit tests.
            expected_class = next(
                n
                for n in tree(adapted(path.read_text(encoding="utf-8"))).body
                if isinstance(n, ast.ClassDef) and n.name == "TestAPIChecker"
            )
            actual_class = next(
                n
                for n in tree(dest.read_text(encoding="utf-8")).body
                if isinstance(n, ast.ClassDef) and n.name == "TestAPIChecker"
            )
            if ast.dump(actual_class) != ast.dump(expected_class):
                errors.append("shared API integration check differs")
            continue
        if not dest.is_file():
            errors.append(f"missing test/fixture: {relative}")
        elif path.suffix == ".py":
            expected = adapted(path.read_text(encoding="utf-8"))
            actual = dest.read_text(encoding="utf-8")
            matches = actual == expected
            if not matches:
                errors.append(f"shared test differs: {relative}")
        elif (
            path.suffix not in (".md", ".out", ".cube")
            and path.read_bytes() != dest.read_bytes()
        ):
            errors.append(f"fixture differs: {relative}")
        elif path.suffix in (".out", ".cube"):
            if path.read_text(encoding="utf-8") != dest.read_text(encoding="utf-8"):
                errors.append(f"fixture differs: {relative}")
    return errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=ROOT.parent / "moleditpy_orca_result_analyzer_plugin",
    )
    args = parser.parse_args()
    failures = check(args.source.resolve())
    if failures:
        raise SystemExit("\n".join(failures))
    print("Canonical shared modules, MO APIs, images and all tests/fixtures match.")

"""
set_version.py — Update the plugin version across all version-bearing files.

Usage:
    python set_version.py 2.5.2
    python set_version.py 2.5.2.1

Files updated:
    orca_result_analyzer_rust/__init__.py   PLUGIN_VERSION = "..."
    pyproject.toml                          version = "..."
    crates/orca_parser_rs/Cargo.toml        version = "..."  (SemVer only — 3-part)
    crates/orca_mo_rs/Cargo.toml            version = "..."

Notes:
    Cargo requires strict SemVer (MAJOR.MINOR.PATCH).  If you supply a 4-part
    version (e.g. 2.5.1.1) the Cargo files are set to the 3-part prefix
    (e.g. 2.5.1) while all Python files receive the full version string.
"""

import re
import sys
import os

ROOT = os.path.dirname(os.path.abspath(__file__))

FILES_PYTHON = [
    os.path.join(ROOT, "orca_result_analyzer_rust", "__init__.py"),
    os.path.join(ROOT, "pyproject.toml"),
]

FILES_CARGO = [
    os.path.join(ROOT, "crates", "orca_parser_rs", "Cargo.toml"),
    os.path.join(ROOT, "crates", "orca_mo_rs",    "Cargo.toml"),
]


def semver_prefix(version: str) -> str:
    """Return the 3-part SemVer prefix of a version string."""
    parts = version.split(".")
    return ".".join(parts[:3])


def replace_in_file(path: str, pattern: str, replacement: str):
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()
    new_content, count = re.subn(pattern, replacement, content, flags=re.MULTILINE)
    if count == 0:
        print(f"  WARNING: pattern not found in {os.path.relpath(path, ROOT)}")
        return
    with open(path, "w", encoding="utf-8") as f:
        f.write(new_content)
    print(f"  {os.path.relpath(path, ROOT)}")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    version = sys.argv[1].strip()
    cargo_version = semver_prefix(version)

    if not re.fullmatch(r"\d+\.\d+\.\d+(\.\d+)?", version):
        print(f"ERROR: version must be X.Y.Z or X.Y.Z.N, got: {version!r}")
        sys.exit(1)

    print(f"Setting plugin version: {version}")
    if cargo_version != version:
        print(f"  (Cargo files will use 3-part prefix: {cargo_version})")

    print("\nPython / TOML files:")
    for path in FILES_PYTHON:
        if path.endswith("__init__.py"):
            replace_in_file(path,
                r'(PLUGIN_VERSION\s*=\s*")[^"]*(")',
                rf'\g<1>{version}\g<2>')
        else:  # pyproject.toml
            replace_in_file(path,
                r'(^version\s*=\s*")[^"]*(")',
                rf'\g<1>{version}\g<2>',)

    print("\nCargo files (3-part SemVer):")
    for path in FILES_CARGO:
        replace_in_file(path,
            r'(^version\s*=\s*")[^"]*(")',
            rf'\g<1>{cargo_version}\g<2>')

    print("\nDone.")


if __name__ == "__main__":
    main()

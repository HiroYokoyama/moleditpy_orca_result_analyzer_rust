"""
build.py — Build both Rust extensions for the ORCA Result Analyzer.

  crates/orca_parser_rs/  ->  orca_parser_rs.*.pyd   (ORCA output parser)
  crates/orca_mo_rs/      ->  orca_mo_rs.*.pyd        (GTO/MO grid evaluation)

Usage:
    python build.py            # development build + copy in-place + create plugin zip
    python build.py --release  # optimised release build
    python build.py --wheel    # produce wheels only (no in-place copy or zip)

Requires:
    pip install maturin
"""

import subprocess
import sys
import argparse
import os
import glob
import shutil
import zipfile

ROOT     = os.path.dirname(os.path.abspath(__file__))
PKG_DIR  = os.path.join(ROOT, "orca_result_analyzer_rust")
BUILD_DIR = os.path.join(ROOT, "build")

CRATES = [
    {
        "manifest": os.path.join(ROOT, "crates", "orca_parser_rs", "Cargo.toml"),
        "pattern":  "orca_parser_rs",
    },
    {
        "manifest": os.path.join(ROOT, "crates", "orca_mo_rs", "Cargo.toml"),
        "pattern":  "orca_mo_rs",
    },
]


def run(cmd, **kwargs):
    print(f">> {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=ROOT, **kwargs)
    if result.returncode != 0:
        sys.exit(result.returncode)


def build_crate(maturin, manifest, release):
    cmd = maturin + [
        "build",
        "--manifest-path", manifest,
        "--out", BUILD_DIR,
        "--interpreter", sys.executable,
    ]
    if release:
        cmd += ["--release"]
    run(cmd)


def extract_pyd(wheel, pattern):
    """Extract only the .pyd/.so matching *pattern* from *wheel* into PKG_DIR."""
    with zipfile.ZipFile(wheel) as zf:
        ext_files = [n for n in zf.namelist()
                     if pattern in n and n.endswith((".pyd", ".so"))]
        if not ext_files:
            print(f"ERROR: Could not find '{pattern}' .pyd/.so in {os.path.basename(wheel)}")
            print("Contents:", zf.namelist())
            sys.exit(1)
        for name in ext_files:
            dest = os.path.join(PKG_DIR, os.path.basename(name))
            with zf.open(name) as src, open(dest, "wb") as dst:
                shutil.copyfileobj(src, dst)
            print(f"Copied  {os.path.basename(name)}  ->  orca_result_analyzer_rust/")


def newest_wheel(pattern):
    wheels = glob.glob(os.path.join(BUILD_DIR, f"*{pattern}*.whl"))
    if not wheels:
        print(f"ERROR: No wheel matching '*{pattern}*.whl' found in build/")
        sys.exit(1)
    return max(wheels, key=os.path.getmtime)


def main():
    ap = argparse.ArgumentParser(description="Build Rust extensions for ORCA Result Analyzer.")
    ap.add_argument("--release", action="store_true",
                    help="Optimised release build (slower compile, faster runtime).")
    ap.add_argument("--wheel",   action="store_true",
                    help="Build wheels only; skip in-place copy and plugin zip.")
    args = ap.parse_args()

    # Check maturin
    check = subprocess.run(
        [sys.executable, "-m", "maturin", "--version"],
        capture_output=True, cwd=ROOT,
    )
    if check.returncode != 0:
        print("ERROR: maturin not found.  Install with:  pip install maturin")
        sys.exit(1)

    maturin = [sys.executable, "-m", "maturin"]
    os.makedirs(BUILD_DIR, exist_ok=True)

    # Build both crates
    for crate in CRATES:
        build_crate(maturin, crate["manifest"], args.release)

    if args.wheel:
        print(f"\nWheels written to: {BUILD_DIR}/")
        return

    # Extract .pyd from each wheel
    for crate in CRATES:
        wheel = newest_wheel(crate["pattern"])
        print(f"\nExtracting from: {os.path.basename(wheel)}")
        extract_pyd(wheel, crate["pattern"])

    # Derive version + platform tag from the parser wheel name
    parser_wheel = newest_wheel("orca_parser_rs")
    stem = os.path.splitext(os.path.basename(parser_wheel))[0]
    parts = stem.split("-")
    if len(parts) >= 5:
        platform_tag = "-".join(parts[2:])   # e.g. cp313-cp313-win_amd64
    else:
        platform_tag = "unknown"

    # Read plugin version from __init__.py by parsing the text, NOT by importing
    # it — __init__.py imports PyQt6, which is absent in the build environment
    # (that made the zip name fall back to "unknown").
    import re
    version_tag = "unknown"
    try:
        with open(os.path.join(PKG_DIR, "__init__.py"), encoding="utf-8") as fh:
            m = re.search(
                r'^\s*PLUGIN_VERSION\s*=\s*["\'](.+?)["\']', fh.read(), re.MULTILINE
            )
        if m:
            version_tag = m.group(1)
    except Exception:
        version_tag = "unknown"

    # Package the plugin folder as a zip for the MoleditPy plugin system
    zip_name = f"orca_result_analyzer_rust-{version_tag}-{platform_tag}.zip"
    zip_path = os.path.join(BUILD_DIR, zip_name)
    pkg_name = os.path.basename(PKG_DIR)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, dirs, files in os.walk(PKG_DIR):
            dirs[:] = [d for d in dirs if d != "__pycache__"]
            for fname in files:
                if fname.endswith(".pyc"):
                    continue
                full    = os.path.join(root, fname)
                arcname = os.path.join(pkg_name, os.path.relpath(full, PKG_DIR))
                zf.write(full, arcname)
    print(f"Plugin zip:  {zip_name}  ->  build/")

    print("\nDone.")


if __name__ == "__main__":
    main()

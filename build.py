"""
build.py — Build the Rust extension (orca_parser_rs) for the ORCA Result Analyzer.

Usage:
    python build.py            # build and copy .pyd/.so in-place (no virtualenv needed)
    python build.py --release  # optimised release build (slower compile, faster runtime)
    python build.py --wheel    # build a redistributable wheel into build/ only

Requires:
    pip install maturin

The compiled extension is copied directly into orca_result_analyzer/ so Python
can import it as:
    from orca_result_analyzer import orca_parser_rs
"""

import subprocess
import sys
import argparse
import os
import glob
import shutil

ROOT = os.path.dirname(os.path.abspath(__file__))
PKG_DIR = os.path.join(ROOT, "orca_result_analyzer_rust")
BUILD_DIR = os.path.join(ROOT, "build")


def run(cmd, **kwargs):
    print(f">> {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=ROOT, **kwargs)
    if result.returncode != 0:
        sys.exit(result.returncode)


def main():
    parser = argparse.ArgumentParser(description="Build the Rust parser extension.")
    parser.add_argument(
        "--release",
        action="store_true",
        help="Build in release mode (optimised, slower compile).",
    )
    parser.add_argument(
        "--wheel",
        action="store_true",
        help="Build a wheel into build/ only, do not copy in-place.",
    )
    args = parser.parse_args()

    # Ensure maturin is available
    check = subprocess.run(
        [sys.executable, "-m", "maturin", "--version"],
        capture_output=True,
        cwd=ROOT,
    )
    if check.returncode != 0:
        print("ERROR: maturin not found.  Install it with:  pip install maturin")
        sys.exit(1)

    maturin = [sys.executable, "-m", "maturin"]
    os.makedirs(BUILD_DIR, exist_ok=True)

    # Always build a wheel into build/ first (works without a virtualenv)
    cmd = maturin + ["build", "--out", BUILD_DIR, "--interpreter", sys.executable]
    if args.release:
        cmd += ["--release"]
    run(cmd)

    if args.wheel:
        print(f"\nWheel written to: {BUILD_DIR}/")
        return

    # Find the built wheel and extract the .pyd / .so from it
    wheels = glob.glob(os.path.join(BUILD_DIR, "*.whl"))
    if not wheels:
        print("ERROR: No wheel found in build/ after maturin build.")
        sys.exit(1)

    wheel = max(wheels, key=os.path.getmtime)  # newest
    print(f"\nExtracting extension from: {os.path.basename(wheel)}")

    import zipfile
    with zipfile.ZipFile(wheel) as zf:
        # The extension file is named orca_parser_rs*.pyd or orca_parser_rs*.so
        ext_files = [n for n in zf.namelist()
                     if "orca_parser_rs" in n and not n.endswith("/")]
        if not ext_files:
            print("ERROR: Could not find orca_parser_rs extension inside wheel.")
            print("Contents:", zf.namelist())
            sys.exit(1)

        for ext_name in ext_files:
            dest = os.path.join(PKG_DIR, os.path.basename(ext_name))
            with zf.open(ext_name) as src, open(dest, "wb") as dst:
                shutil.copyfileobj(src, dst)
            print(f"Copied  {os.path.basename(ext_name)}  →  orca_result_analyzer/")

    # Derive platform tag from the wheel name, e.g. cp313-cp313-win_amd64
    wheel_stem = os.path.splitext(os.path.basename(wheel))[0]  # strip .whl
    # wheel name: <dist>-<ver>-<pytag>-<abitag>-<platform>.whl
    parts = wheel_stem.split("-")
    if len(parts) >= 5:
        version_tag = parts[1]               # e.g. 2.4.1
        platform_tag = "-".join(parts[2:])   # e.g. cp313-cp313-win_amd64
    else:
        version_tag = "unknown"
        platform_tag = "unknown"

    # Package the plugin folder as a zip for the MoleditPy plugin system
    import zipfile as _zf
    zip_name = f"orca_result_analyzer_rust-{version_tag}-{platform_tag}.zip"
    zip_path = os.path.join(BUILD_DIR, zip_name)
    pkg_name = os.path.basename(PKG_DIR)
    with _zf.ZipFile(zip_path, "w", _zf.ZIP_DEFLATED) as zf:
        for root, dirs, files in os.walk(PKG_DIR):
            # Skip __pycache__
            dirs[:] = [d for d in dirs if d != "__pycache__"]
            for fname in files:
                if fname.endswith(".pyc"):
                    continue
                full = os.path.join(root, fname)
                arcname = os.path.join(pkg_name, os.path.relpath(full, PKG_DIR))
                zf.write(full, arcname)
    print(f"Plugin zip:  {zip_name}  →  build/")

    print("\nDone. The Rust parser is ready.")


if __name__ == "__main__":
    main()

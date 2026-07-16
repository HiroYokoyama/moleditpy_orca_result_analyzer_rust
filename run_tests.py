#!/usr/bin/env python3
"""CI-faithful test runner.

Runs the suite the way CI does: with third-party pytest plugin autoload
disabled. The tests stub PyQt6/PySide6 at module level; a locally-installed
pytest-qt otherwise imports a real binding at collection and defeats those
stubs, causing spurious failures/segfaults. Disabling autoload keeps local
runs matching CI regardless of what is installed.

Requires the Rust extensions to be built first (see build.py).

Usage:
    python run_tests.py                # run the whole suite
    python run_tests.py -k parser      # forward any args to pytest
"""
import os
import subprocess
import sys
from pathlib import Path

os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"

ROOT = Path(__file__).resolve().parent
args = sys.argv[1:] or ["tests/"]
raise SystemExit(subprocess.call([sys.executable, "-m", "pytest", *args], cwd=ROOT))

PLUGIN_NAME = "ORCA Result Analyzer"
PLUGIN_VERSION = "3.9.2"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Comprehensive analyzer for ORCA output files (.out). Includes Vibrational, MO, TDDFT, and NMR analysis."
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"

from PyQt6.QtWidgets import QApplication, QMessageBox  # noqa: E402
import logging  # noqa: E402

_context = None  # Stored from initialize() so run() can use the registry API


def _read_orca_file(path, parent_widget):
    """Read an ORCA output file trying several encodings. Returns content string or None."""
    encodings = ["utf-8", "utf-16", "latin-1", "cp1252"]
    for enc in encodings:
        try:
            with open(path, "r", encoding=enc) as f:
                return f.read()
        except UnicodeError:
            continue
        except Exception as e:
            QMessageBox.critical(
                parent_widget, "Error Reading File", f"Could not read file:\n{e}"
            )
            return None
    # Fallback with error replace
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            return f.read()
    except Exception as e:
        QMessageBox.critical(
            parent_widget, "Error Reading File", f"Could not read file:\n{e}"
        )
        return None


def _open_orca_file(path, context):
    """Parse and display an ORCA output file using the context registry for window management."""
    QApplication.processEvents()
    mw = context.get_main_window()

    content = _read_orca_file(path, mw)
    if content is None:
        return

    from .parser import OrcaParser

    parser = OrcaParser()
    parser.load_from_memory(content, path)

    # Close existing window if open
    existing = context.get_window("analyzer")
    if existing is not None:
        try:
            existing.close()
        except Exception as _e:
            logging.warning("silenced: %s", _e)

    # New result loaded — any atom colors applied to the previous molecule's
    # indices must not bleed onto this (possibly differently-indexed) one.
    from .utils import clear_atom_color_overrides

    clear_atom_color_overrides(mw)

    from .gui import OrcaResultAnalyzerDialog

    win = OrcaResultAnalyzerDialog(mw, parser, path, context)
    context.register_window("analyzer", win)

    win.show()
    win.raise_()
    win.activateWindow()

    QApplication.processEvents()
    # fit_camera=True: this is a fresh file load (drop or file-opener), so
    # reframe the camera after drawing. Without it the newly loaded molecule
    # keeps the previous view's zoom/center and can fall outside the sphere.
    win.load_structure_3d(fit_camera=True)
    QApplication.processEvents()


def _open_orca_analyzer_empty(context):
    """Open the analyzer dialog with no file pre-loaded."""
    QApplication.processEvents()
    mw = context.get_main_window()

    # Close existing window if open
    existing = context.get_window("analyzer")
    if existing is not None:
        try:
            # If already open, just raise it
            existing.show()
            existing.raise_()
            existing.activateWindow()
            return
        except Exception as _e:
            logging.warning("silenced: %s", _e)

    from .parser import OrcaParser
    from .gui import OrcaResultAnalyzerDialog
    from .utils import clear_atom_color_overrides

    parser = OrcaParser()  # empty — no file loaded
    clear_atom_color_overrides(mw)

    win = OrcaResultAnalyzerDialog(mw, parser, "", context)
    context.register_window("analyzer", win)

    win.show()
    win.raise_()
    win.activateWindow()
    QApplication.processEvents()


def _on_document_reset(context):
    """Close the analyzer window and drop stale 3D overlay state on File->New.

    File->New wipes the plotter (mw.clear_3d_view()) but does not know about
    plugin-owned state: the analyzer window would otherwise keep showing the
    previous file's data, and any atom-color overrides applied by the Atomic
    Charges view would survive and could bleed onto a differently-indexed
    molecule loaded afterwards.
    """
    win = context.get_window("analyzer")
    if win is not None:
        try:
            win.close()
        except Exception as e:
            logging.warning("silenced: %s", e)

    from .utils import clear_atom_color_overrides

    clear_atom_color_overrides(context.get_main_window())


def initialize(context):
    """Initialize the ORCA Result Analyzer plugin.

    Registers file openers for .out with HIGH PRIORITY (100).
    """
    global _context
    _context = context

    def open_orca_file(path):
        _open_orca_file(path, context)

    def handle_drop(path):
        # Check for standard ORCA output
        if path.lower().endswith(".out"):
            try:
                with open(path, "r", encoding="utf-8", errors="ignore") as f:
                    header = f.read(2048)
                if "* O   R   C   A *" in header or "Program Version" in header:
                    _open_orca_file(path, context)
                    return True
            except Exception as _e:
                logging.warning("silenced: %s", _e)
        return False

    context.register_file_opener(".out", open_orca_file, priority=100)
    context.register_drop_handler(handle_drop, priority=100)
    context.register_document_reset_handler(lambda: _on_document_reset(context))

    # Extensions menu entry — opens analyzer without requiring a file first
    context.add_menu_action(
        "Extensions/ORCA Result Analyzer",
        lambda: _open_orca_analyzer_empty(context),
    )


def run(mw):
    """Legacy run() entry: called from Plugins menu by the host.

    Opens the analyzer directly without requiring a file — the user can
    use 'Select File' or 'Select from Directory' inside the dialog.
    """
    context = _context
    if context is None:
        return

    _open_orca_analyzer_empty(context)

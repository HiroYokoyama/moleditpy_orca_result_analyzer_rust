"""Plugin entry points and metadata for the ORCA Result Analyzer plugin."""

PLUGIN_NAME = "ORCA Result Analyzer"
PLUGIN_VERSION = "5.0.2"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Comprehensive analyzer for ORCA output files (.out). Includes Vibrational, MO, TDDFT, and NMR analysis."
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=5.0.0, <6.0.0"
PLUGIN_SUPPORTED_PYTHON_VERSION = ">=3.9, <3.14"
PLUGIN_DEPENDENCIES = ["rdkit", "PySide6", "matplotlib", "Pillow"]
PLUGIN_OPTIONAL_DEPENDENCIES = ["nmrsim"]

from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402
import logging  # noqa: E402

_context = None  # Stored from initialize() so run() can use the registry API


def _read_orca_file(path, parent_widget):
    """Read an ORCA output file trying several encodings. Returns content string or None."""
    from .loading import read_orca_text

    try:
        return read_orca_text(path)
    except (OSError, ValueError) as e:
        QMessageBox.critical(
            parent_widget, "Error Reading File", f"Could not read file:\n{e}"
        )
        return None


def _open_orca_file(path, context):
    """Parse and display an ORCA output file using the context registry for window management."""
    QApplication.processEvents()
    mw = context.get_main_window()

    from .loading import load_orca_parser

    try:
        parser = load_orca_parser(path, mw)
    except (OSError, ValueError) as e:
        QMessageBox.critical(mw, "Error Reading File", f"Could not read file:\n{e}")
        return
    if parser is None:  # cancelled by the user
        return

    # Close existing window if open
    existing = context.get_window("analyzer")
    if existing is not None:
        try:
            if existing.close() is False:
                return
        except (RuntimeError, AttributeError) as e:
            logging.warning("Could not close the previous analyzer window: %s", e)

    # New result loaded — any atom colors applied to the previous molecule's
    # indices must not bleed onto this (possibly differently-indexed) one.
    from .utils import clear_atom_color_overrides, sync_main_window_file

    clear_atom_color_overrides(mw)
    sync_main_window_file(mw, path, context)

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
        except (RuntimeError, AttributeError) as e:
            logging.warning("Could not re-raise the existing analyzer window: %s", e)

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
        except (RuntimeError, AttributeError) as e:
            logging.warning(
                "Could not close the analyzer window on document reset: %s", e
            )

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
            except (OSError, ValueError) as e:
                logging.warning(
                    "Could not read dropped file %s to check for an ORCA header: %s",
                    path,
                    e,
                )
        return False

    context.register_file_opener(".out", open_orca_file, priority=100)
    context.register_drop_handler(handle_drop, priority=100)
    context.register_document_reset_handler(lambda: _on_document_reset(context))

    # Extensions menu entry — opens analyzer without requiring a file first
    context.add_menu_action(
        "Extensions/ORCA Result Analyzer",
        lambda: _open_orca_analyzer_empty(context),
    )


def run(mw):  # pylint: disable=unused-argument  # host calls run(mw); reopens via module context
    """Legacy run() entry: called from Plugins menu by the host.

    Opens the analyzer directly without requiring a file — the user can
    use 'Select File' or 'Select from Directory' inside the dialog.
    """
    context = _context
    if context is None:
        return

    _open_orca_analyzer_empty(context)

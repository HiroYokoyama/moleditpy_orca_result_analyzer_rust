PLUGIN_NAME = "ORCA Result Analyzer (Rust)"
PLUGIN_VERSION = "2.5.1.2"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Comprehensive analyzer for ORCA output files (.out) with Rust-powered parser. Includes Vibrational, MO, TDDFT, and NMR analysis."

from PyQt6.QtWidgets import QMessageBox, QFileDialog
import logging

# Global reference to keep window alive
_analyzer_window = None


def initialize(context):
    """
    Initialize the ORCA Result Analyzer plugin.
    Registers file openers for .out with HIGH PRIORITY (100).
    """

    def open_orca_file(path):
        global _analyzer_window
        from PyQt6.QtWidgets import QApplication

        # Ensure main window is initialized and processed
        # This helps when starting from CLI
        QApplication.processEvents()

        mw = context.get_main_window()

        # Read file to memory
        content = ""
        encodings = ["utf-8", "utf-16", "latin-1", "cp1252"]
        for enc in encodings:
            try:
                with open(path, "r", encoding=enc) as f:
                    content = f.read()
                break
            except UnicodeError:
                continue
            except Exception as e:
                # If it's not an encoding error, standard error
                QMessageBox.critical(
                    mw, "Error Reading File", f"Could not read file:\n{e}"
                )
                return

        if not content:
            # Fallback with error replace
            try:
                with open(path, "r", encoding="utf-8", errors="replace") as f:
                    content = f.read()
            except Exception as e:
                QMessageBox.critical(
                    mw, "Error Reading File", f"Could not read file:\n{e}"
                )
                return

        # Initialize Parser
        from .parser import OrcaParser

        parser = OrcaParser()
        parser.load_from_memory(content, path)

        # Close existing if open
        if _analyzer_window is not None:
            try:
                _analyzer_window.close()
            except Exception as _e:
                logging.warning("silenced: %s", _e)
            _analyzer_window = None

        # Open Dialog (Modeless)
        from .gui import OrcaResultAnalyzerDialog

        _analyzer_window = OrcaResultAnalyzerDialog(mw, parser, path, context)

        _analyzer_window.show()
        _analyzer_window.raise_()
        _analyzer_window.activateWindow()

        # Ensure UI shows up before logic
        QApplication.processEvents()

        # Auto-load 3D structure
        _analyzer_window.load_structure_3d()

        # Final UI flush
        QApplication.processEvents()

    # Register Drop Handler
    def handle_drop(path):
        # 1. Check for standard ORCA output
        if path.lower().endswith(".out"):
            # content check for ORCA
            try:
                with open(path, "r", encoding="utf-8", errors="ignore") as f:
                    header = f.read(2048)  # Read first 2KB
                    if "* O   R   C   A *" in header or "Program Version" in header:
                        open_orca_file(path)
                        return True
            except ImportError:
                pass
            except Exception as _e:
                logging.warning("silenced: %s", _e)
        return False

    # Register Opener
    # Priority 100
    context.register_file_opener(".out", open_orca_file, priority=100)

    context.register_drop_handler(handle_drop, priority=100)

    # Add to Main Menu
    def menu_action():
        mw = context.get_main_window()
        path, _ = QFileDialog.getOpenFileName(
            mw, "Open ORCA Output", "", "ORCA Output (*.out)"
        )
        if path:
            if not handle_drop(path):
                open_orca_file(
                    path
                )  # Fallback to standard opener logic if handle_drop returns false (e.g. extension check)

    # context.add_menu_action("Analysis/ORCA Result Analyzer", menu_action)


def run(mw):
    from PyQt6.QtWidgets import QFileDialog, QApplication

    path, _ = QFileDialog.getOpenFileName(
        mw, "Open ORCA Output", "", "ORCA Output (*.out);;All Files (*)"
    )
    if not path:
        return

    from moleditpy.plugins.plugin_interface import PluginContext

    context = PluginContext(mw.plugin_manager, PLUGIN_NAME)

    global _analyzer_window
    QApplication.processEvents()

    content = ""
    encodings = ["utf-8", "utf-16", "latin-1", "cp1252"]
    for enc in encodings:
        try:
            with open(path, "r", encoding=enc) as f:
                content = f.read()
            break
        except UnicodeError:
            continue
        except Exception as e:
            QMessageBox.critical(mw, "Error Reading File", f"Could not read file:\n{e}")
            return

    if not content:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
        except Exception as e:
            QMessageBox.critical(mw, "Error Reading File", f"Could not read file:\n{e}")
            return

    from .parser import OrcaParser

    parser = OrcaParser()
    parser.load_from_memory(content, path)

    if _analyzer_window is not None:
        try:
            _analyzer_window.close()
        except Exception as _e:
            logging.warning("silenced: %s", _e)
        _analyzer_window = None

    from .gui import OrcaResultAnalyzerDialog

    _analyzer_window = OrcaResultAnalyzerDialog(mw, parser, path, context)
    _analyzer_window.show()
    _analyzer_window.raise_()
    _analyzer_window.activateWindow()
    QApplication.processEvents()
    _analyzer_window.load_structure_3d()
    QApplication.processEvents()

"""Optional hand-off to the ORCA NICS Analyzer plugin.

NICS analysis lives in a separate plugin, so this module only ever probes for
it: when that plugin is not installed the menu entry stays hidden.
"""

import logging
import sys

NICS_PLUGIN_NAME = "ORCA NICS Analyzer"
NICS_MODULE_HINT = "orca_nics_analyzer"


def _looks_like_nics(module) -> bool:
    if module is None:
        return False
    if getattr(module, "PLUGIN_NAME", "") == NICS_PLUGIN_NAME:
        return True
    return NICS_MODULE_HINT in str(getattr(module, "__name__", "")).lower()


def find_nics_module(main_window):
    """Return the loaded NICS Analyzer plugin module, or None."""
    manager = getattr(main_window, "plugin_manager", None)
    try:
        for entry in getattr(manager, "plugins", None) or []:
            if not isinstance(entry, dict):
                continue
            module = entry.get("module")
            name = str(entry.get("name", ""))
            if _looks_like_nics(module) or NICS_MODULE_HINT in name.lower():
                return module
    except (AttributeError, TypeError, RuntimeError) as exc:
        logging.debug("NICS plugin lookup via the plugin manager failed: %s", exc)

    # Second chance: the module is imported but the host registry is unusable.
    try:
        for module in list(sys.modules.values()):
            if (
                _looks_like_nics(module)
                and getattr(module, "_context", None) is not None
            ):
                return module
    except (AttributeError, TypeError, RuntimeError) as exc:
        logging.debug("NICS plugin lookup via sys.modules failed: %s", exc)
    return None


def nics_analyzer_available(main_window) -> bool:
    """Return whether the NICS Analyzer plugin is loaded in the host."""
    return find_nics_module(main_window) is not None


def _trigger_menu_action(main_window, exclude=None) -> bool:
    """Fire the NICS plugin's own menu entry (used when we have no file).

    *exclude* is the analyzer's own "NICS Analysis..." action: findChildren is
    recursive and the analyzer window is a child of the main window, so without
    it the search could find and re-trigger the entry that called us.
    """
    try:
        from PySide6.QtGui import QAction

        for action in main_window.findChildren(QAction):
            if exclude is not None and action is exclude:
                continue
            if "nics" in action.text().lower():
                action.trigger()
                return True
    except (ImportError, AttributeError, TypeError, RuntimeError) as exc:
        logging.debug("Triggering the NICS menu action failed: %s", exc)
    return False


def open_nics_analyzer(main_window, file_path=None, exclude_action=None):
    """Open the NICS Analyzer, on *file_path* when one is loaded.

    Returns ``(ok, message)``; *message* explains the failure when not ok.
    *exclude_action* is never triggered by the menu fallback.
    """
    module = find_nics_module(main_window)
    if module is None:
        return False, (
            "The ORCA NICS Analyzer plugin is not installed.\n\n"
            "Install it from the Plugin Manager to analyze NICS data."
        )

    opener = getattr(module, "open_file", None) or getattr(module, "_open_file", None)
    context = getattr(module, "_context", None)
    if file_path and callable(opener) and context is not None:
        try:
            opener(file_path, context)
            return True, ""
        except Exception as exc:  # plugin boundary: never take the analyzer down  # pylint: disable=broad-exception-caught
            logging.warning("NICS Analyzer hand-off failed", exc_info=True)
            return False, f"The ORCA NICS Analyzer could not read this file:\n{exc}"

    if _trigger_menu_action(main_window, exclude=exclude_action):
        return True, ""
    return False, "The ORCA NICS Analyzer plugin could not be opened."

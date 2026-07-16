import logging
import os

# Known dummy / pseudo-atom labels used in ORCA and other QC output files.
_DUMMY_SYMBOLS: frozenset[str] = frozenset(
    {"*", "-", "X", "DA", "DU", "DUM", "DUMMY", "Q", "BQ", "LP"}
)


def get_default_export_path(base_path, suffix="_analyzed", extension=""):
    """
    Generates a default export path based on the input file path.
    Example: 'job.out' -> 'job_analyzed.csv'
    """
    if not base_path:
        return ""
    dirname = os.path.dirname(base_path)
    filename_base = os.path.splitext(os.path.basename(base_path))[0]
    new_filename = f"{filename_base}{suffix}{extension}"
    return os.path.join(dirname, new_filename)


def normalize_atom_symbol(raw: str) -> str:
    """Return a valid RDKit atom symbol for *raw*, mapping dummy labels to '*'.

    Handles:
    - Labels with a colon suffix  e.g. ``X:1`` → ``*``
    - Known dummy labels          e.g. ``DA``, ``DU``, ``BQ`` → ``*``
    - Unknown / non-periodic      e.g. ``Xx`` → ``*``
    - Normal elements             e.g. ``Fe``, ``C`` → returned as-is (capitalised)
    """
    sym = raw.strip()
    # Strip ORCA-style colon-suffixed labels like "X:1" or "C:2"
    if ":" in sym:
        sym = sym.split(":")[0]
    if sym.upper() in _DUMMY_SYMBOLS:
        return "*"
    sym = sym.capitalize()
    try:
        from rdkit import Chem

        if Chem.GetPeriodicTable().GetAtomicNumber(sym) <= 0:
            return "*"
    except Exception:  # noqa: BLE001
        return "*"
    return sym


def determine_bonds_without_dummies(mol, charge: int = 0, bond_orders: bool = True):
    """Run RDKit bond determination on *mol*, skipping dummy ('*') atoms.

    Builds a sub-molecule containing only real (non-dummy) atoms, calls
    ``DetermineConnectivity`` (and optionally ``DetermineBondOrders``) on it,
    then copies the resulting bonds back to the original *mol* (which must be
    an ``RWMol``).  Any failure is caught and logged — the function is
    intentionally non-fatal.

    Parameters
    ----------
    mol:
        An RDKit ``RWMol`` with a conformer already attached.
    charge:
        Formal charge to pass to ``DetermineBondOrders``.
    bond_orders:
        If *True* (default) also determine bond orders.  Pass *False*
        during animation playback to avoid per-frame latency.
    """
    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D
        from rdkit.Chem import rdDetermineBonds

        conf = mol.GetConformer()

        # Build index map: sub_idx -> orig_idx (only real atoms)
        real_indices = [
            i
            for i in range(mol.GetNumAtoms())
            if mol.GetAtomWithIdx(i).GetSymbol() != "*"
        ]

        if not real_indices:
            return  # nothing to do

        # Build sub-molecule with only real atoms
        sub = Chem.RWMol()
        sub_conf = Chem.Conformer(len(real_indices))
        for sub_i, orig_i in enumerate(real_indices):
            sym = mol.GetAtomWithIdx(orig_i).GetSymbol()
            sub.AddAtom(Chem.Atom(sym))
            pos = conf.GetAtomPosition(orig_i)
            sub_conf.SetAtomPosition(sub_i, Point3D(pos.x, pos.y, pos.z))
        sub.AddConformer(sub_conf)

        # Bond determination on the sub-molecule
        rdDetermineBonds.DetermineConnectivity(sub)
        if bond_orders:
            try:
                rdDetermineBonds.DetermineBondOrders(sub, charge=charge)
            except Exception as e:
                logging.debug(
                    "DetermineBondOrders failed, falling back to connectivity: %s", e
                )

        # Copy bonds back to the original molecule
        for bond in sub.GetBonds():
            orig_i = real_indices[bond.GetBeginAtomIdx()]
            orig_j = real_indices[bond.GetEndAtomIdx()]
            mol.AddBond(orig_i, orig_j, bond.GetBondType())

    except Exception as exc:  # noqa: BLE001
        logging.debug("determine_bonds_without_dummies: non-fatal — %s", exc)


def clear_atom_color_overrides(mw) -> None:
    """Clear any 3D atom-color overrides left by the Atomic Charges view.

    Loading a *different* result (new file, reload, or a fresh analyzer
    window) must not let colors keyed to the previous molecule's atom
    indices bleed onto a differently-indexed molecule. Called on every
    "new file" entry point (drag-drop / Select File / Reload / opening a
    fresh analyzer from the main GUI) as well as on document reset.
    """
    v3d = getattr(mw, "view_3d_manager", None) if mw is not None else None
    if v3d is not None and hasattr(v3d, "_plugin_color_overrides"):
        try:
            v3d._plugin_color_overrides.clear()
        except Exception as exc:  # noqa: BLE001
            logging.warning("clear_atom_color_overrides: %s", exc)


def list_orca_output_files(directory: str) -> list[str]:
    """Return a sorted list of ``*.out`` filenames found in *directory*.

    Only the bare filenames (not full paths) are returned.  An empty list is
    returned if *directory* does not exist or cannot be listed.

    Parameters
    ----------
    directory:
        Path to the directory to scan.
    """
    try:
        return sorted(f for f in os.listdir(directory) if f.lower().endswith(".out"))
    except Exception as exc:  # noqa: BLE001
        logging.debug("list_orca_output_files: cannot list '%s' — %s", directory, exc)
        return []

"""Parsing mixin for MO coefficients, orbital energies, basis set and SCF trace."""

import re
import logging


class _ElectronicParsingMixin:
    """MO coefficient, orbital energy, basis set and SCF trace parsing for OrcaParser."""

    def parse_mo_coeffs(self) -> None:
        """Extract per-orbital coefficients, energy, occupancy and spin into self.data."""
        self.data[
            "mo_coeffs"
        ] = {}  # mo_idx -> { 'coeffs': list, 'energy': float, 'occ': float, 'spin': 'alpha'/'beta' }

        # Blocks to look for:
        # "MOLECULAR ORBITALS" (RHF or UHF if (UHF) present)
        # "SPIN UP ORBITALS" (UHF Alpha)
        # "SPIN DOWN ORBITALS" (UHF Beta)

        start_indices = []
        for i, line in enumerate(self.lines):
            uu = line.upper()
            if (
                "MOLECULAR ORBITALS" in uu
                and i + 1 < len(self.lines)
                and "---" in self.lines[i + 1]
            ):
                start_indices.append((i, "restricted"))
            elif (
                "SPIN UP ORBITALS" in uu
                and i + 1 < len(self.lines)
                and "---" in self.lines[i + 1]
            ):
                start_indices.append((i, "alpha"))
            elif (
                "SPIN DOWN ORBITALS" in uu
                and i + 1 < len(self.lines)
                and "---" in self.lines[i + 1]
            ):
                start_indices.append((i, "beta"))

        if not start_indices:
            return

        generic_starts = [idx for idx, spin in start_indices if spin == "restricted"]
        if generic_starts:
            start_indices = [
                (idx, spin) for idx, spin in start_indices if idx >= generic_starts[-1]
            ]

        # Only keep the LAST occurrence for each spin type
        final_indices = {}
        for idx, spin in start_indices:
            final_indices[spin] = idx

        filtered_indices = []
        for spin, idx in final_indices.items():
            filtered_indices.append((idx, spin))
        filtered_indices.sort()

        for start_idx, spin in filtered_indices:
            curr = start_idx + 2

            # Check if this "restricted" block is actually UHF
            header_line = self.lines[start_idx].upper()
            current_spin = spin
            if spin == "restricted" and "(UHF)" in header_line:
                current_spin = "alpha"

            last_first_mo_idx = -1
            current_mos = []

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    curr += 1
                    continue
                if "TIMINGS" in line:
                    break
                if "--------" in line:
                    curr += 1
                    continue
                if (
                    "ORBITALS" in line
                    and curr + 1 < len(self.lines)
                    and "--------" in self.lines[curr + 1]
                ):
                    break  # Next block

                parts = line.split()
                if not parts:
                    curr += 1
                    continue

                # Header check: Integers? "0   1   2..."
                is_header = all(p.isdigit() for p in parts)

                if is_header:
                    current_mos = [int(p) for p in parts]

                    # Detect spin switch (Index Reset) logic for implicit UHF
                    if current_mos[0] <= last_first_mo_idx and last_first_mo_idx != -1:
                        # If we were in alpha/restricted and index dropped, assume beta
                        if current_spin in ("alpha", "restricted"):
                            current_spin = "beta"

                    last_first_mo_idx = current_mos[0]

                    # Init storage dicts
                    for idx in current_mos:
                        key = f"{idx}_{current_spin}"  # Use current_spin
                        # If duplicate key (e.g. from previous block override), reset it
                        # BUT be careful not to wipe if we are appending chunks (not the case here, cols are complete MOs)
                        if key not in self.data["mo_coeffs"]:
                            self.data["mo_coeffs"][key] = {
                                "coeffs": [],
                                "spin": current_spin,
                                "energy": 0.0,
                                "occ": 0.0,
                                "id": idx,
                            }

                    # Try to parse Energy / Occ lines immediately following
                    # Usually:
                    #                  -10.000    -5.000
                    #                   2.0000     1.000
                    try:
                        next_lines = [
                            self.lines[curr + 1].strip(),
                            self.lines[curr + 2].strip(),
                        ]
                        # Check if numbers
                        vals1 = next_lines[0].split()
                        vals2 = next_lines[1].split()

                        if len(vals1) == len(current_mos) and len(vals2) == len(
                            current_mos
                        ):
                            # Assume line 1 is Energy (Eh), Line 2 is Occ
                            # Verify they look like floats
                            try:
                                energies = [float(v) for v in vals1]
                                occs = [float(v) for v in vals2]
                                for k, idx in enumerate(current_mos):
                                    key = f"{idx}_{current_spin}"
                                    if key in self.data["mo_coeffs"]:
                                        self.data["mo_coeffs"][key]["energy"] = (
                                            energies[k]
                                        )
                                        self.data["mo_coeffs"][key]["occ"] = occs[k]
                                curr += 2  # Skip these 2 lines
                            except (KeyError, IndexError, TypeError, ValueError) as e:
                                logging.warning(
                                    "MO coefficients: could not parse energy/occupation values for MOs %s: %s",
                                    current_mos,
                                    e,
                                )
                    except (
                        AttributeError,
                        KeyError,
                        IndexError,
                        TypeError,
                        ValueError,
                    ) as e:
                        logging.warning(
                            "MO coefficients: could not read the energy/occupation lines after MOs %s: %s",
                            current_mos,
                            e,
                        )

                    curr += 1
                    continue

                # Coefficient line: "0   C  1s   0.000 ..." or "0C  1s ..."
                if len(parts) >= 2:
                    try:
                        atom_idx = -1
                        sym = ""
                        orb = ""
                        val_strs = []

                        # Check for merged format "0C"
                        match_merged = re.match(r"^(\d+)([A-Za-z]+)$", parts[0])
                        if match_merged:
                            atom_idx = int(match_merged.group(1))
                            sym = match_merged.group(2)
                            orb = parts[1]
                            val_strs = parts[2:]
                        elif len(parts) >= 3 and parts[0].isdigit():
                            atom_idx = int(parts[0])
                            sym = parts[1]
                            orb = parts[2]
                            val_strs = parts[3:]
                        else:
                            curr += 1
                            continue

                        val_strs = re.findall(
                            r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?",
                            " ".join(val_strs),
                        )
                        if len(val_strs) == len(current_mos):
                            for k, v_str in enumerate(val_strs):
                                mo_idx = current_mos[k]
                                key = f"{mo_idx}_{current_spin}"
                                try:
                                    val = float(v_str)
                                    if key in self.data["mo_coeffs"]:
                                        self.data["mo_coeffs"][key]["coeffs"].append(
                                            {
                                                "atom_idx": atom_idx,
                                                "sym": sym,
                                                "orb": orb,
                                                "coeff": val,
                                            }
                                        )
                                except (
                                    KeyError,
                                    IndexError,
                                    TypeError,
                                    ValueError,
                                ) as e:
                                    logging.warning(
                                        "MO coefficients: could not parse coefficient %r for atom %d orbital %s: %s",
                                        v_str,
                                        atom_idx,
                                        orb,
                                        e,
                                    )
                    except (
                        AttributeError,
                        KeyError,
                        IndexError,
                        TypeError,
                        ValueError,
                    ) as e:
                        logging.warning(
                            "MO coefficients: could not parse the coefficient line %r: %s",
                            line,
                            e,
                        )
                curr += 1

    def parse_orbital_energies(self):
        """Parse orbital energies from ORBITAL ENERGIES section"""
        self.data["orbital_energies"] = []
        self.data["mos"] = []  # Backward compatibility

        # Look for "ORBITAL ENERGIES" section
        start_indices = []
        for i, line in enumerate(self.lines):
            uu = line.upper()
            if (
                "ORBITAL ENERGIES" in uu
                and i + 1 < len(self.lines)
                and "---" in self.lines[i + 1]
            ):
                start_indices.append((i, "restricted"))  # Default
            elif "SPIN UP ORBITALS" in uu:
                # Real ORCA prints the "NO OCC E(Eh) E(eV)" column header on the
                # next line (not a "---" separator), so do not require one here.
                start_indices.append((i, "alpha"))
            elif "SPIN DOWN ORBITALS" in uu:
                start_indices.append((i, "beta"))

        if not start_indices:
            return

        generic_starts = [idx for idx, spin in start_indices if spin == "restricted"]
        if generic_starts:
            start_indices = [
                (idx, spin) for idx, spin in start_indices if idx >= generic_starts[-1]
            ]

        # Only keep the LAST occurrence for each spin type
        final_indices = {}
        for idx, spin in start_indices:
            final_indices[spin] = idx

        # Open-shell/UHF: explicit spin sections exist. The generic
        # "ORBITAL ENERGIES" header points at the alpha block and would
        # otherwise be read as a single "restricted" block (swallowing the
        # alpha orbitals and dropping the entire beta manifold), so discard it
        # in favor of the dedicated alpha/beta blocks.
        if "alpha" in final_indices or "beta" in final_indices:
            final_indices.pop("restricted", None)

        filtered_indices = []
        for spin, idx in final_indices.items():
            filtered_indices.append((idx, spin))
        filtered_indices.sort()

        for start_idx, spin in filtered_indices:
            # Skip the section title, then find the column header. The header
            # sits at start_idx+1 for spin blocks (no separator) and at
            # start_idx+2 for the restricted block (after the "---" rule), so
            # scan from start_idx+1 to cover both.
            curr = start_idx + 1
            while curr < len(self.lines) and curr < start_idx + 10:
                line = self.lines[curr].strip()
                # Check for "NO", "OCC", and "Eh" or "eV"
                if "NO" in line and (
                    "OCC" in line or "E(Eh)" in line or "E(eV)" in line
                ):
                    curr += 1
                    break
                curr += 1

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if (
                    not line
                    or "---" in line
                    or "****" in line
                    or "MULLIKEN" in line
                    or "ORBITALS" in line.upper()
                ):
                    break
                if line.startswith("*"):
                    curr += 1
                    continue

                parts = line.split()
                # Format: NO   OCC          E(Eh)            E(eV)
                #           0   2.0000     -10.186408      -277.1862
                if len(parts) >= 4 and parts[0].lstrip("-").isdigit():
                    try:
                        orbital_idx = int(parts[0])
                        occupation = float(parts[1])
                        energy_eh = float(parts[2])
                        energy_ev = float(parts[3])

                        orb_data = {
                            "index": orbital_idx,
                            "id": orbital_idx,  # For mo_analysis
                            "occupation": occupation,
                            "occ": occupation,  # For mo_analysis
                            "energy_eh": energy_eh,
                            "energy_ev": energy_ev,
                            "energy": energy_eh,  # For backward compatibility
                            "spin": spin,
                            "type": "occupied" if occupation > 0.1 else "virtual",
                        }
                        self.data["orbital_energies"].append(orb_data)
                        self.data["mos"].append(orb_data)
                    except (
                        AttributeError,
                        KeyError,
                        IndexError,
                        TypeError,
                        ValueError,
                    ) as e:
                        logging.warning(
                            "Orbital energies: could not parse line %r: %s", line, e
                        )

                curr += 1

    def parse_basis_set(self):
        """Parse Basis Set information needed for MO visualization"""
        self.data["basis_set_shells"] = []

        # Look for "BASIS SET IN INPUT FORMAT"
        start_idx = -1
        for i, line in enumerate(self.lines):
            if line.strip() == "BASIS SET IN INPUT FORMAT":
                start_idx = i

        if start_idx == -1:
            return

        curr = start_idx + 2

        basis_defs = {}  # Sym -> List of shells
        current_sym = None
        current_shells = []

        while curr < len(self.lines):
            line = self.lines[curr].strip()

            # Stop conditions
            if "--------" in line and curr > start_idx + 10:
                break
            if "AUXILIARY BASIS" in line:
                break

            # Start of Atom block: "NewGTO H"
            if line.startswith("NewGTO"):
                parts = line.split()
                if len(parts) >= 2:
                    current_sym = parts[1]
                    current_shells = []
                curr += 1
                continue

            # End of Atom block: "end" or "end;"
            if line.startswith("end"):
                if current_sym:
                    basis_defs[current_sym] = current_shells
                curr += 1
                continue

            # Shell definition header: "S   3" or "P   2"
            parts = line.split()
            if len(parts) >= 2 and parts[0].upper() in ["S", "P", "D", "F", "G"]:
                sh_type = parts[0].upper()
                try:
                    n_prim = int(parts[1])
                    if n_prim > 50:  # Sanity check
                        curr += 1
                        continue

                    curr += 1

                    exps = []
                    coeffs = []

                    # Read primitives
                    for _ in range(n_prim):
                        if curr >= len(self.lines):
                            break
                        pl = self.lines[curr].strip()
                        pp = pl.split()
                        if len(pp) >= 3:
                            exps.append(float(pp[1]))
                            coeffs.append(float(pp[2]))
                        curr += 1

                    l_map = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}
                    l_val = l_map.get(sh_type, 0)

                    if exps:
                        current_shells.append(
                            {"l": l_val, "exps": exps, "coeffs": coeffs}
                        )

                    continue
                except (
                    AttributeError,
                    KeyError,
                    IndexError,
                    TypeError,
                    ValueError,
                ) as e:
                    logging.warning(
                        "Basis set: could not parse the %s shell for atom %s: %s",
                        sh_type,
                        current_sym,
                        e,
                    )

            curr += 1
            if curr > start_idx + 5000:
                break  # Safety break

        # Expand to actual atoms
        atoms = self.data.get("atoms", [])
        coords = self.data.get("coords", [])

        full_shells = []

        if not atoms:
            # Try to recover atoms from parser data if parse_basic hasn't run or failed?
            # No, parse_basis_set relies on atoms being parsed.
            return

        for idx, (sym, coord) in enumerate(zip(atoms, coords)):
            defs = basis_defs.get(sym, [])
            for d in defs:
                full_shells.append(
                    {
                        "atom_idx": idx,
                        "origin": coord,
                        "l": d["l"],
                        "exps": d["exps"],
                        "coeffs": d["coeffs"],
                    }
                )

        self.data["basis_set_shells"] = full_shells

    def parse_scf_trace(self):
        """Extract SCF iteration energies for each block found."""
        self.data["scf_traces"] = []

        # We search for blocks starting with D-I-I-S or S-O-S-C-F or SCF ITERATIONS
        current_step_label = "Initial"

        i = 0
        while i < len(self.lines):
            line = self.lines[i]
            uu = line.upper()

            if "OPTIMIZATION CYCLE" in uu:
                try:
                    parts = line.split()
                    # Find the index of CYCLE and take the next part
                    cycle_part = (
                        parts[parts.index("CYCLE") + 1]
                        if "CYCLE" in parts
                        else parts[-2]
                    )
                    current_step_label = f"Cycle {cycle_part}"
                except IndexError:
                    current_step_label = "Opt Cycle"
            elif "SCAN STEP" in uu:
                # The banner is "*   RELAXED SURFACE SCAN STEP   1   *"; the
                # last token is the closing asterisk, not the step number.
                m_step = re.search(r"SCAN STEP\s+(\d+)", uu)
                current_step_label = (
                    f"Scan Step {m_step.group(1)}" if m_step else "Scan Step"
                )
            elif "ORCA PROPERTIES" in uu or "ORCA PROPERTY" in uu:
                current_step_label = "Property/Final"
            elif "OPTIMIZATION HAS CONVERGED" in uu:
                current_step_label = "Post-Opt/Final"

            if (
                "SCF ITERATIONS" in uu
                or "ORCA LEAN-SCF" in uu
                or "INCREMENTAL FOCK MATRIX" in uu
            ):
                header_idx = -1
                for k in range(1, 15):
                    if i + k >= len(self.lines):
                        break
                    uu_k = self.lines[i + k].upper()
                    if "ITER" in uu_k and "ENERGY" in uu_k:
                        header_idx = i + k
                        break

                if header_idx != -1:
                    trace = []
                    idx = header_idx + 1
                    # Check for separator line right after header
                    if idx < len(self.lines) and "---" in self.lines[idx]:
                        idx += 1

                    while idx < len(self.lines):
                        l_scf = self.lines[idx].strip()
                        if not l_scf or "SUCCESS" in l_scf or "Energy Check" in l_scf:
                            if trace:
                                break
                            idx += 1
                            continue

                        parts = l_scf.split()
                        if len(parts) >= 2:
                            try:
                                it_no = int(parts[0])
                                it_en = float(parts[1])
                                trace.append({"iter": it_no, "energy": it_en})
                            except (IndexError, TypeError, ValueError):
                                # ORCA prints '***' for overflow values; skip unparseable lines
                                logging.debug(
                                    "Skipping unparseable SCF iteration line: %r",
                                    l_scf,
                                    exc_info=True,
                                )
                        idx += 1

                    if trace:
                        # Check if we should append or start new
                        # If it's a new "SCF ITERATIONS" block, we want a new entry in the traces list
                        # unless it's genuinely part of the same convergence process (rare in ORCA output stream)

                        # Heuristic: if last trace in self.data["scf_traces"] has same label,
                        # suffix the new one if they are distinct blocks.
                        # Exact match or an earlier "(n)" suffix of it: a bare
                        # startswith() counted "Cycle 10".."Cycle 19" as
                        # repeats of "Cycle 1".
                        same_label_count = sum(
                            1
                            for t in self.data["scf_traces"]
                            if t["step"] == current_step_label
                            or t["step"].startswith(f"{current_step_label} (")
                        )

                        label = current_step_label
                        if same_label_count > 0:
                            label = f"{current_step_label} ({same_label_count + 1})"

                        self.data["scf_traces"].append(
                            {"step": label, "iterations": trace}
                        )
                        i = idx
                        continue
            i += 1

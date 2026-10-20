"""Parsing mixin for NMR, TD-DFT, thermochemistry and vibrational frequencies."""

import re
import logging


#: Modes below this magnitude are translations/rotations or numerical noise,
#: not genuine imaginary frequencies (cm^-1).
IMAGINARY_FREQ_THRESHOLD = 10.0


class _SpectraParsingMixin:
    """NMR, TD-DFT, thermochemistry and frequency parsing for OrcaParser."""

    def parse_nmr(self) -> None:
        """Extract NMR shielding and spin-spin coupling values into self.data."""
        self.data[
            "nmr_shielding"
        ] = []  # List of {atom_idx, atom_sym, shielding, shift=None}
        self.data["nmr_couplings"] = []

        # Look for "CHEMICAL SHIELDING SUMMARY (PPM)" or individual nucleus blocks
        # Last block wins, like every other section: an optimization + NMR
        # job, or a compound job, can print more than one summary.
        summary_start = -1
        for i, line in enumerate(self.lines):
            if "CHEMICAL SHIELDING SUMMARY (PPM)" in line.upper():
                summary_start = i

        if summary_start != -1:
            curr = summary_start + 1
            # Skip until we hit the header "N Nucleus Shielding" or similar
            while curr < len(self.lines) and curr < summary_start + 10:
                l_up = self.lines[curr].upper()
                if "N" in l_up and "SHIELDING" in l_up:
                    curr += 1
                    break
                if "NUCLEUS" in l_up and ("ISOTROPIC" in l_up or "ELEMENT" in l_up):
                    curr += 1
                    break
                curr += 1

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    break
                if "---" in line:
                    if len(self.data["nmr_shielding"]) > 0:
                        break
                    curr += 1
                    continue

                parts = line.split()
                # Format: Index Nucleus Shielding
                if len(parts) >= 3:
                    try:
                        idx = int(parts[0])
                        sym = parts[1]
                        val = float(parts[2])
                        self.data["nmr_shielding"].append(
                            {"atom_idx": idx, "atom_sym": sym, "shielding": val}
                        )
                    except (KeyError, IndexError, TypeError, ValueError) as e:
                        logging.warning(
                            "NMR: could not parse the shielding row %r: %s", line, e
                        )
                curr += 1

        # Parse Couplings
        # "SUMMARY OF ISOTROPIC COUPLING CONSTANTS (Hz)"
        header_found = False
        start_idx = -1
        for i, line in enumerate(self.lines):
            # Make (Hz) optional and case-insensitive
            upper_line = line.upper()
            if "SUMMARY OF ISOTROPIC COUPLING CONSTANTS" in upper_line:
                start_idx = i
                header_found = True

        if header_found:
            curr = start_idx + 1
            # Skip until first separator or data
            while curr < len(self.lines):
                if "----------------" in self.lines[curr]:
                    curr += 1  # Skip first separator
                    break
                curr += 1

            # Now we iterate through blocks
            current_col_indices = []  # List of atom indices for current block of columns

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    curr += 1
                    continue

                if (
                    "Maximum memory used" in line
                    or "Timings" in line
                    or "ORCA TERMINATED" in line
                ):
                    break

                parts = line.split()
                if len(parts) == 0:
                    curr += 1
                    continue

                # Is it a header line?
                is_header = False
                # Header format: "0 C   1 C ..."
                # Row format: "0 C   0.000 ..."

                if len(parts) >= 3:
                    # Check 3rd token (index 2).
                    token3 = parts[2]
                    # In header, token3 is an index (int). In data, it is a value (float).
                    # simple check: isdigit handles positive integers correctly.
                    # 0.000 is not digit. -1.23 is not digit.
                    if token3.isdigit():
                        is_header = True
                    else:
                        is_header = False
                elif len(parts) == 2:
                    # "0 C" only -> could be header of last column block
                    is_header = True

                if is_header:
                    # Parse column indices
                    # ["0", "C", "1", "C", ...]
                    current_col_indices = []
                    p_idx = 0
                    while p_idx < len(parts) - 1:
                        if parts[p_idx].isdigit():
                            current_col_indices.append(int(parts[p_idx]))
                            p_idx += 2
                        else:
                            p_idx += 1
                else:
                    # Parse Row
                    if len(parts) >= 2 and parts[0].isdigit():
                        try:
                            row_atom_idx = int(parts[0])

                            values = parts[2:]
                            for c_i, val_str in enumerate(values):
                                if c_i < len(current_col_indices):
                                    col_atom_idx = current_col_indices[c_i]
                                    val = float(val_str)

                                    # Store unique couplings (J_AB)
                                    if row_atom_idx < col_atom_idx:
                                        self.data["nmr_couplings"].append(
                                            {
                                                "atom_idx1": row_atom_idx,
                                                "atom_idx2": col_atom_idx,
                                                "coupling": val,
                                            }
                                        )
                        except (
                            AttributeError,
                            KeyError,
                            IndexError,
                            TypeError,
                            ValueError,
                        ) as e:
                            logging.warning(
                                "NMR: could not parse a spin-spin coupling value from %r: %s",
                                line,
                                e,
                            )

                curr += 1

    def parse_tddft(self):
        """
        Parse TD-DFT/TDA excited states and spectra (Absorption/CD).
        Updated Version:
        - Correctly parses both Length (Electric Dipole) and Velocity gauges.
        - Stores data in 'osc_len', 'osc_vel', 'rot_len', 'rot_vel'.
        - Handles '0-1A -> 1-1A' transition lines and auto-calculates eV from nm.
        """
        self.data["tddft"] = []
        states_dict = {}  # state_id (int) -> dict

        # Helper to ensure state dict exists with all gauge keys initialized
        def get_state(idx):
            if idx not in states_dict:
                states_dict[idx] = {
                    "state": idx,
                    "energy_ev": 0.0,
                    "energy_nm": 0.0,
                    "energy_cm": 0.0,
                    "osc": 0.0,  # Default (usually Length)
                    "osc_len": 0.0,  # Length Gauge
                    "osc_vel": 0.0,  # Velocity Gauge
                    "rotatory_strength": 0.0,  # Default (usually Length)
                    "rot_len": 0.0,  # Length Gauge
                    "rot_vel": 0.0,  # Velocity Gauge
                    "transitions": [],
                }
            return states_dict[idx]

        # -------------------------------------------------------------------------
        # PASS 1: Parse Detailed Excited States (Energy + Transitions)
        # -------------------------------------------------------------------------
        current_state_id = -1

        # Only the last TD-DFT run counts (an excited-state optimization
        # prints one per cycle); states from earlier runs used to linger
        # whenever a later run reported fewer roots.
        header_re = re.compile(
            r"^[A-Z0-9/\-\s]*EXCITED STATES(\s*\((SINGLETS|TRIPLETS)\))?$"
        )
        singlet_heads = [
            i
            for i, line in enumerate(self.lines)
            if header_re.match(line.strip()) and "TRIPLETS" not in line
        ]
        run_start = singlet_heads[-1] if singlet_heads else 0

        # Triplet states restart at STATE 1 and used to overwrite the singlet
        # of the same number; they carry no oscillator strength, so skip them.
        in_triplets = False
        for i in range(run_start, len(self.lines)):
            line = self.lines[i]
            line_strip = line.strip()
            line_upper = line.upper()

            if header_re.match(line_strip):
                in_triplets = "TRIPLETS" in line_strip
                current_state_id = -1
                continue
            if in_triplets:
                continue

            # Detect State Header
            match_state = re.search(r"STATE\s+(\d+)\s*:", line_upper)
            if match_state:
                try:
                    current_state_id = int(match_state.group(1))
                    state_entry = get_state(current_state_id)

                    # Extract Energies
                    match_ev = re.search(r"([-+]?\d*\.\d+)\s*eV", line, re.IGNORECASE)
                    if match_ev:
                        state_entry["energy_ev"] = float(match_ev.group(1))

                    match_cm = re.search(r"([-+]?\d*\.\d+)\s*cm", line, re.IGNORECASE)
                    if match_cm:
                        state_entry["energy_cm"] = float(match_cm.group(1))

                    match_nm = re.search(r"([-+]?\d*\.\d+)\s*nm", line, re.IGNORECASE)
                    if match_nm:
                        state_entry["energy_nm"] = float(match_nm.group(1))
                except (AttributeError, KeyError, IndexError, TypeError, ValueError):
                    current_state_id = -1

            # Detect Transitions
            elif current_state_id != -1:
                if "-------" in line or "SPECTRUM" in line_upper:
                    pass
                elif "->" in line and ":" in line:
                    parts = line_strip.split(":")
                    if len(parts) >= 2:
                        trans_desc = parts[0].strip()
                        coeff = parts[1].strip()
                        t_str = f"{trans_desc} (coeff: {coeff})"
                        if t_str not in get_state(current_state_id)["transitions"]:
                            get_state(current_state_id)["transitions"].append(t_str)

        # -------------------------------------------------------------------------
        # PASS 2: Parse Summary Tables (All Gauges)
        # -------------------------------------------------------------------------

        def parse_summary_table(start_idx, data_key):
            curr = start_idx + 1
            header_found = False

            # 1. Header detection (only confirms the table has started)
            while curr < len(self.lines) and curr < start_idx + 30:
                line = self.lines[curr].strip()
                if not line or "--------" in line:
                    curr += 1
                    continue

                u_line = line.upper()
                # A table needs TRANSITION/STATE plus ENERGY/WAVELENGTH
                if ("TRANSITION" in u_line or "STATE" in u_line) and (
                    "ENERGY" in u_line or "WAVELENGTH" in u_line
                ):
                    header_found = True
                    break

                # Fallback for older output formats
                parts = line.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    header_found = True
                    break

                curr += 1

            if not header_found:
                return

            # 2. Parse the data rows
            data_parsing_started = False
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    curr += 1
                    continue

                # Check for "TOTAL" explicitly (always stop)
                if "TOTAL" in line.upper():
                    break

                # Skip the first separator line encountered (header separator)
                if "--------" in line:
                    if not data_parsing_started:
                        curr += 1
                        # If the very next line is also a separator, skip it too?
                        # Or if we just skipped one and haven't started data, we are now ready.
                        # Let's set a flag to verify we found the header separator.
                        data_parsing_started = True
                        continue
                    else:
                        break  # End of table

                parts = line.split()
                # If we encounter a valid data line (starts with digit), mark as started
                if (len(parts) > 0 and parts[0].isdigit()) or "->" in parts:
                    data_parsing_started = True

                # --- Pattern A: standard format containing the "->" arrow ---
                # e.g. 0-1A  ->  1-1A   2.780   22422   446.0   0.000 ...
                # idx:  0    1     2      3       4       5       6
                if "->" in parts:
                    try:
                        arrow_idx = parts.index("->")
                        # Columns are read relative to the arrow position
                        # arrow+1: Target State (1-1A)
                        # arrow+2: Energy (eV)
                        # arrow+3: Energy (cm-1)
                        # arrow+4: Wavelength (nm)
                        # arrow+5: Value (Strength)

                        if len(parts) > arrow_idx + 5:
                            # 1. State ID. "0-1A -> 1-3A" is a
                            # spin-forbidden triplet whose number clashes
                            # with singlet 1; skip multiplicity changes.
                            target_state_str = parts[arrow_idx + 1]
                            src = re.match(r"^\d+-(\d+)", parts[arrow_idx - 1])
                            dst = re.match(r"^\d+-(\d+)", target_state_str)
                            if src and dst and src.group(1) != dst.group(1):
                                match = None
                            else:
                                match = re.search(r"^(\d+)", target_state_str)
                            if match:
                                s_id = int(match.group(1))
                                entry = get_state(s_id)

                                # 2. Values (intensity)
                                entry[data_key] = float(parts[arrow_idx + 5])

                                # 3. Energies (fill in when zero / overwrite)
                                # Reading eV here matters when the detailed block (PASS 1) was absent
                                try:
                                    entry["energy_ev"] = float(parts[arrow_idx + 2])
                                except (IndexError, TypeError, ValueError) as e:
                                    logging.debug(
                                        "TD-DFT: could not parse the eV energy for state %d: %s",
                                        s_id,
                                        e,
                                    )

                                try:
                                    entry["energy_cm"] = float(parts[arrow_idx + 3])
                                except (IndexError, TypeError, ValueError) as e:
                                    logging.debug(
                                        "TD-DFT: could not parse the cm**-1 energy for state %d: %s",
                                        s_id,
                                        e,
                                    )

                                try:
                                    entry["energy_nm"] = float(parts[arrow_idx + 4])
                                except (IndexError, TypeError, ValueError) as e:
                                    logging.debug(
                                        "TD-DFT: could not parse the nm wavelength for state %d: %s",
                                        s_id,
                                        e,
                                    )
                    except (
                        AttributeError,
                        KeyError,
                        IndexError,
                        TypeError,
                        ValueError,
                    ):
                        # Skip rows that fail to parse
                        logging.debug("Skipping unparseable TDDFT line", exc_info=True)

                # --- Pattern B: short format, no arrow ---
                #   State  Energy  Wavelength  Value ...
                #   idx 0     1         2        3
                # Column 1 is cm-1 in ORCA <= 5 but eV in some layouts, and
                # both label it only in a second header line. Assuming eV put
                # a ~22000 eV excitation in the table whenever the detailed
                # STATE block was absent, so infer the unit from column 2 (nm),
                # which both layouts agree on.
                elif len(parts) >= 4 and parts[0].isdigit():
                    try:
                        s_id = int(parts[0])
                        entry = get_state(s_id)

                        entry[data_key] = float(parts[3])

                        col1 = float(parts[1])
                        nm = float(parts[2])
                        if entry["energy_nm"] == 0:
                            entry["energy_nm"] = nm

                        if nm > 0.1 and col1 > 0:
                            as_cm = 1e7 / nm
                            as_ev = 1239.84193 / nm
                            col1_is_cm = abs(col1 - as_cm) < abs(col1 - as_ev)
                            if col1_is_cm:
                                if entry["energy_cm"] == 0:
                                    entry["energy_cm"] = col1
                                if entry["energy_ev"] == 0:
                                    entry["energy_ev"] = as_ev
                            elif entry["energy_ev"] == 0:
                                entry["energy_ev"] = col1
                        elif entry["energy_ev"] == 0:
                            entry["energy_ev"] = col1
                    except (
                        AttributeError,
                        KeyError,
                        IndexError,
                        TypeError,
                        ValueError,
                    ) as e:
                        logging.warning(
                            "TD-DFT: could not parse the short-format excitation row %r: %s",
                            line,
                            e,
                        )

                curr += 1

        # Locate tables and parse SPECIFIC gauges, from the last run only.
        # SOC-corrected spectra number spin-orbit states differently and used
        # to overwrite the spin-free strengths of the same index.
        for i in range(run_start, len(self.lines)):
            line_upper = self.lines[i].upper()
            if "SOC" in line_upper or "SPIN ORBIT" in line_upper:
                continue

            # Absorption Spectrum
            if "ABSORPTION SPECTRUM" in line_upper:
                if "ELECTRIC DIPOLE" in line_upper:
                    # Parse as 'osc_len' AND default 'osc'
                    parse_summary_table(i, "osc_len")
                    parse_summary_table(i, "osc")
                elif "VELOCITY DIPOLE" in line_upper:
                    # Parse as 'osc_vel'
                    parse_summary_table(i, "osc_vel")

            # CD Spectrum
            if "CD SPECTRUM" in line_upper:
                if "ELECTRIC DIPOLE" in line_upper:
                    # Parse as 'rot_len' AND default 'rotatory_strength'
                    parse_summary_table(i, "rot_len")
                    parse_summary_table(i, "rotatory_strength")
                elif "VELOCITY DIPOLE" in line_upper:
                    # Parse as 'rot_vel'
                    parse_summary_table(i, "rot_vel")

        # -------------------------------------------------------------------------
        # FINALIZATION: Fix Missing Units & Sort
        # -------------------------------------------------------------------------
        all_items = list(states_dict.values())

        # 1. Auto-Calculate Missing eV from nm
        for item in all_items:
            if item["energy_ev"] == 0.0:
                if item["energy_nm"] > 0.1:
                    item["energy_ev"] = 1239.84193 / item["energy_nm"]
                elif item["energy_cm"] > 0.1:
                    item["energy_ev"] = item["energy_cm"] / 8065.54425

        # 2. Filter & Sort
        valid_items = [item for item in all_items if item["energy_ev"] > 0]
        valid_items.sort(key=lambda x: x["energy_ev"])

        self.data["tddft"] = valid_items

    def parse_thermal(self):
        """Extract the thermochemistry block (enthalpy, entropy, free energy) into self.data."""
        self.data["thermal"] = {}
        # ORCA Thermochem block
        # Look for "THERMOCHEMISTRY AT 298.15 K"

        start_line = -1
        # Scan for the start of thermochemistry section
        for i, line in enumerate(self.lines):
            line_upper = line.upper()
            if "THERMOCHEMISTRY AT" in line_upper:
                start_line = i

        if start_line != -1:
            curr = start_line

            # Mapping of ORCA labels to JSON keys
            # Using uppercase for case-insensitive matching
            thermal_keys = {
                "ELECTRONIC ENERGY": "electronic_energy",
                "ZERO POINT ENERGY": "zpe",
                "THERMAL VIBRATIONAL CORRECTION": "corr_vib",
                "THERMAL ROTATIONAL CORRECTION": "corr_rot",
                "THERMAL TRANSLATIONAL CORRECTION": "corr_trans",
                "TOTAL THERMAL ENERGY": "thermal_energy",
                "TOTAL THERMAL CORRECTION": "corr_thermal_total",
                "NON-THERMAL (ZPE) CORRECTION": "corr_zpe",
                "TOTAL CORRECTION": "corr_total",
                "TOTAL ENTHALPY": "enthalpy",
                "THERMAL ENTHALPY CORRECTION": "thermal_enthalpy_corr",  # Just RT
                "ELECTRONIC ENTROPY": "s_el",
                "VIBRATIONAL ENTROPY": "s_vib",
                "ROTATIONAL ENTROPY": "s_rot",
                "TRANSLATIONAL ENTROPY": "s_trans",
                "FINAL ENTROPY TERM": "entropy",
                "FINAL GIBBS FREE ENERGY": "gibbs",
                "G-E(EL)": "gibbs_corr",
            }

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                line_upper = line.upper()

                if "TIMINGS FOR INDIVIDUAL MODULES" in line_upper:
                    break

                # Temperature check
                if (
                    "TEMPERATURE" in line_upper
                    and "K" in line_upper
                    and "..." in line_upper
                ):
                    match = re.search(r"(\d+\.\d+)\s*K", line)
                    if match:
                        self.data["thermal"]["temperature"] = float(match.group(1))

                # Iterate through expected keys
                for key_upper, val_key in thermal_keys.items():
                    if key_upper in line_upper:
                        # Extract the numeric value, prioritizing heartrees (Eh)
                        if "Eh" in line:
                            # Split by "Eh" and take the part before it
                            pre_eh = line.split("Eh")[0]
                            # Find the last float-like number in the pre-Eh part
                            val_matches = re.findall(
                                r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?", pre_eh
                            )
                            if val_matches:
                                try:
                                    val = float(val_matches[-1])
                                    self.data["thermal"][val_key] = val
                                except ValueError as e:
                                    logging.warning(
                                        "Thermal: could not parse %s value from %r: %s",
                                        val_key,
                                        pre_eh,
                                        e,
                                    )
                        else:
                            # Fallback: take the last numeric match
                            matches = re.findall(
                                r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?", line
                            )
                            if matches:
                                try:
                                    val = float(matches[-1])
                                    self.data["thermal"][val_key] = val
                                except ValueError as e:
                                    logging.warning(
                                        "Thermal: could not parse %s value from %r: %s",
                                        val_key,
                                        line,
                                        e,
                                    )
                curr += 1

            # Post-processing: Calculate H correction (H - E_el) more accurately
            # ORCA's "Thermal Enthalpy correction" is often just RT.
            # We want the total correction including ZPE and Thermal effects.
            t_data = self.data["thermal"]
            if "enthalpy" in t_data and "electronic_energy" in t_data:
                t_data["enthalpy_corr"] = (
                    t_data["enthalpy"] - t_data["electronic_energy"]
                )
            elif "corr_total" in t_data and "thermal_enthalpy_corr" in t_data:
                # Fallback: H_corr = Total_corr (U-E_el) + RT_corr
                t_data["enthalpy_corr"] = (
                    t_data["corr_total"] + t_data["thermal_enthalpy_corr"]
                )

            # Similarly for Gibbs if not directly found
            if "gibbs" in t_data and "electronic_energy" in t_data:
                t_data["gibbs_corr"] = t_data["gibbs"] - t_data["electronic_energy"]

            # Count imaginary frequencies. Unprojected translations/rotations
            # come through as small negative values; counting those reported a
            # transition state for a perfectly good minimum.
            freqs = self.data.get("frequencies", [])
            imaginary_count = sum(
                1 for f in freqs if f.get("freq", 0) < -IMAGINARY_FREQ_THRESHOLD
            )
            t_data["imaginary_freq_count"] = imaginary_count

            return

        # 2. Final Energy
        for line in reversed(self.lines):
            target = "FINAL SINGLE POINT ENERGY"
            if target in line:
                try:
                    parts = line.split()
                    val = float(parts[-1])
                    self.data["scf_energy"] = val
                    break
                except (KeyError, IndexError, TypeError, ValueError) as e:
                    logging.warning(
                        "Thermal: could not parse the final single-point energy from %r: %s",
                        line,
                        e,
                    )

    def parse_frequencies(self):
        """Extract vibrational frequencies with their IR/Raman intensities and vectors."""
        self.data["frequencies"] = []

        # 1. Frequencies
        freq_start = -1
        for i, line in enumerate(self.lines):
            if "VIBRATIONAL FREQUENCIES" in line.upper():
                freq_start = i
                # Don't break immediately, could be multiple? usually last one matters or first?
                # In optimization + freq, it's at end.

        if freq_start != -1:
            curr = freq_start + 1
            # Skip until we find first data line (Index: Value)
            while curr < len(self.lines) and curr < freq_start + 10:
                line = self.lines[curr].strip()
                if ":" in line and ("cm**-1" in line or "cm-1" in line):
                    break  # Found data match
                curr += 1

            # Now parse data
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "NORMAL MODES" in line or "-------" in line:
                    # If we already have freqs, a dashed line means end
                    if len(self.data["frequencies"]) > 0:
                        break

                if not line:
                    curr += 1
                    continue

                # Format:   0:         0.00 cm**-1
                if ":" in line and ("cm**-1" in line or "cm-1" in line):
                    parts = line.split()
                    # Find value after colon?
                    # usually parts: "0:", "0.00", "cm**-1"
                    try:
                        # Value is typically index 1 if index 0 ends with colon
                        val_str = parts[1]
                        val = float(val_str)
                        self.data["frequencies"].append(
                            {"freq": val, "ir": 0.0, "raman": 0.0, "vector": []}
                        )
                    except (KeyError, IndexError, TypeError, ValueError) as e:
                        logging.warning(
                            "Frequencies: could not parse the vibrational frequency row %r: %s",
                            line,
                            e,
                        )
                elif len(self.data["frequencies"]) > 0 and ":" not in line:
                    # Maybe end of block
                    break

                curr += 1

        # IR, Raman and normal modes are read only after the last frequency
        # block, so a later run without (say) an IR table cannot inherit an
        # earlier run's intensities onto its own modes.
        if freq_start == -1:
            return
        after_freqs = range(freq_start + 1, len(self.lines))

        # 2. IR Spectrum
        ir_start = -1
        for i in after_freqs:
            if "IR SPECTRUM" in self.lines[i].upper():
                ir_start = i

        if ir_start != -1:
            # Locate the km/mol column from the header rather than assuming
            # position 3: the layout differs between ORCA versions, and a
            # mismatch used to be swallowed into ir=0.0, which is
            # indistinguishable from a genuinely IR-inactive mode.
            ir_col = self._find_intensity_column(ir_start, ("int", "t**2"), default=3)

            curr = ir_start + 5
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "The first frequency" in line or "-----" in line:
                    if "The first frequency" in line:
                        break
                    if "-----" in line and curr > ir_start + 10:
                        break  # End of block

                parts = line.split()
                if len(parts) > ir_col and ":" in parts[0]:
                    try:
                        idx = int(parts[0].replace(":", ""))
                        inten = float(parts[ir_col])
                        if 0 <= idx < len(self.data["frequencies"]):
                            self.data["frequencies"][idx]["ir"] = inten
                            deriv = self._parse_dipole_derivative(line)
                            if deriv is not None:
                                self.data["frequencies"][idx]["dipole_deriv"] = deriv
                    except ValueError:
                        logging.debug(
                            "IR row not parsed at column %d: %r", ir_col, line
                        )
                curr += 1

        # 3. Raman
        raman_start = -1
        for i in after_freqs:
            if "RAMAN SPECTRUM" in self.lines[i].upper():
                raman_start = i

        if raman_start != -1:
            raman_col = self._find_intensity_column(
                raman_start, ("activity",), default=2
            )

            curr = raman_start + 5
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "The first frequency" in line or "-----" in line:
                    if "The first frequency" in line:
                        break
                    if "-----" in line and curr > raman_start + 10:
                        break

                parts = line.split()
                if len(parts) > raman_col and ":" in parts[0]:
                    try:
                        idx = int(parts[0].replace(":", ""))
                        act = float(parts[raman_col])
                        if 0 <= idx < len(self.data["frequencies"]):
                            self.data["frequencies"][idx]["raman"] = act
                    except ValueError:
                        logging.debug(
                            "Raman row not parsed at column %d: %r", raman_col, line
                        )
                curr += 1

        # 4. Normal Modes
        modes_start = -1
        for i in after_freqs:
            if "NORMAL MODES" in self.lines[i].upper():
                modes_start = i

        if modes_start != -1 and self.data["atoms"]:
            n_atoms = len(self.data["atoms"])
            n_coords = n_atoms * 3
            curr = modes_start + 7

            mode_buffer = {}  # m_idx -> [vals]

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    curr += 1
                    continue
                if "IR SPECTRUM" in line or "--------" in line:
                    break

                try:
                    headers = [int(x) for x in line.split()]
                    start_data = curr + 1
                    for r in range(n_coords):
                        if start_data + r >= len(self.lines):
                            break
                        dline = self.lines[start_data + r]
                        dparts = dline.split()
                        values = [float(x) for x in dparts[1:]]

                        for c, m_idx in enumerate(headers):
                            if m_idx not in mode_buffer:
                                mode_buffer[m_idx] = []
                            if c < len(values):
                                mode_buffer[m_idx].append(values[c])

                    curr = start_data + n_coords
                except ValueError:
                    curr += 1
                    continue

            # Process collected mode buffer AFTER parsing all mode blocks
            for m_idx, vec_flat in mode_buffer.items():
                if 0 <= m_idx < len(self.data["frequencies"]):
                    vecs = []
                    for k in range(0, len(vec_flat), 3):
                        if k + 2 < len(vec_flat):
                            vecs.append((vec_flat[k], vec_flat[k + 1], vec_flat[k + 2]))
                    self.data["frequencies"][m_idx]["vector"] = vecs

    def _find_intensity_column(self, section_start, wanted, default):
        """Index of the first column in `wanted` within a spectrum header.

        ORCA's IR/Raman tables carry a "Mode freq ..." header whose columns
        move between versions. The header names the mode and frequency
        columns too, so its token positions line up with the data rows —
        but only once inline unit tokens are dropped. The Raman header reads
        "Mode freq (cm**-1) Activity Depolarization" while its rows carry no
        unit column, so counting "(cm**-1)" put every lookup one place to the
        right: Raman activity was read from the depolarization column, and
        ORCA 4's "Mode freq (cm**-1) T**2 ..." IR header resolved to TX.

        Falls back to `default` when no header is recognized.
        """
        for i in range(section_start + 1, min(section_start + 6, len(self.lines))):
            tokens = self.lines[i].split()
            if not tokens or tokens[0].lower() != "mode":
                continue
            lowered = [t.lower().strip(":") for t in tokens if not t.startswith("(")]
            for name in wanted:
                if name in lowered:
                    return lowered.index(name)
        return default

    @staticmethod
    def _parse_dipole_derivative(line):
        """The (TX TY TZ) triple an IR row carries, or None.

        ORCA prints the mode's transition dipole derivative in parentheses
        after the intensity columns -- "( 0.040181 -0.009108  0.002223)" --
        in atomic units. Older versions omit it entirely, so a row without
        the group is normal rather than an error. Matched on the parenthesis
        group instead of fixed columns because the leading columns move
        between ORCA versions (the same reason ir_col is looked up).
        """
        m = re.search(r"\(([^)]*)\)", line)
        if not m:
            return None
        parts = m.group(1).split()
        if len(parts) != 3:
            return None
        try:
            return tuple(float(p) for p in parts)
        except ValueError:
            logging.debug("Unparseable dipole derivative group: %r", line)
            return None

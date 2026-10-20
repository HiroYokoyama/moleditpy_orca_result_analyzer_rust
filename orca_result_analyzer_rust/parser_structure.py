"""Parsing mixin for geometry, trajectory, gradients and relaxed-scan steps."""

from __future__ import annotations

import re
import logging

_XYZ_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


class _StructureParsingMixin:
    """Geometry, trajectory, gradient and scan parsing for OrcaParser."""

    @staticmethod
    def _parse_xyz_row(line: str) -> tuple[str, list[float]] | None:
        """(symbol, [x, y, z]) from a "SYM X Y Z" row, or None if it is not one.

        A malformed row used to raise ValueError straight out of parse_all and
        fail the whole file load.
        """
        parts = line.split()
        if len(parts) < 4:
            return None
        try:
            return parts[0], [float(parts[1]), float(parts[2]), float(parts[3])]
        except ValueError:
            logging.debug("Not a coordinate row: %r", line)
            return None

    def parse_xyz_content(self, content: str) -> list[dict]:
        """Parse multi-frame XYZ content."""
        lines = content.splitlines()
        steps = []
        i = 0
        n_lines = len(lines)

        while i < n_lines:
            line = lines[i].strip()
            if not line:
                i += 1
                continue

            # Number of atoms
            try:
                natoms = int(line)
            except ValueError:
                i += 1
                continue

            if natoms <= 0:
                logging.warning("XYZ: ignoring invalid atom count %d", natoms)
                i += 1
                continue

            i += 1
            if i >= n_lines:
                break

            # Comment line (extract energy if possible)
            comment = lines[i].strip()

            is_excluded = False
            upper_comment = comment.upper()

            # Check TS
            if re.search(r"\bTS\b", upper_comment):
                is_excluded = True

            # Check CI, but carefully
            # Only exclude if it's NOT "CI-NEB" or similar method string
            if re.search(r"\bCI\b", upper_comment):
                # It contains CI. Check if it is part of CI-NEB
                if "CI-NEB" not in upper_comment:
                    is_excluded = True

            if is_excluded:
                # Skip this step (atoms + coords)
                # We need to advance i by natoms
                i += 1 + natoms
                continue

            # Try robust extraction of Energy and Distance/Coord from comment
            energy = 0.0
            dist_val = None

            # 1. Look for Energy Label: "Energy: -123.4" or "Energy -123.4" or "E -123.4"
            # The label must be a whole word: without \b the trailing "e" of a
            # preceding word (e.g. "Coordinate 1.2") matches the bare-"E"
            # alternative case-insensitively and steals the wrong number.
            e_match = re.search(
                rf"\b(?:Energy|E)\b[=:\s]+({_XYZ_NUMBER})",
                comment,
                re.IGNORECASE,
            )
            if e_match:
                try:
                    energy = float(e_match.group(1))
                except (IndexError, TypeError, ValueError) as e:
                    logging.debug(
                        "XYZ: could not parse the energy label from comment %r: %s",
                        comment,
                        e,
                    )
            else:
                # Fallback: Just take the last float (usually energy)
                floats = re.findall(_XYZ_NUMBER, comment)
                if floats:
                    try:
                        energy = float(floats[-1])
                    except (IndexError, TypeError, ValueError) as e:
                        logging.debug(
                            "XYZ: could not parse a fallback energy value from comment %r: %s",
                            comment,
                            e,
                        )

            # 2. Look for Distance/Coordinate Label: "Dist 1.2" or "Coord 1.2"
            d_match = re.search(
                rf"\b(?:Dist(?:ance)?|Coord(?:inate)?|Scan)[=:\s]+({_XYZ_NUMBER})",
                comment,
                re.IGNORECASE,
            )
            if d_match:
                try:
                    dist_val = float(d_match.group(1))
                except (IndexError, TypeError, ValueError) as e:
                    logging.debug(
                        "XYZ: could not parse the scan coordinate from comment %r: %s",
                        comment,
                        e,
                    )

            i += 1

            atoms = []
            coords = []

            for _ in range(natoms):
                if i >= n_lines:
                    break
                parts = lines[i].split()
                if len(parts) >= 4:
                    # A single malformed coordinate line must not abort the
                    # whole multi-frame parse — skip just that atom.
                    try:
                        xyz = [float(parts[1]), float(parts[2]), float(parts[3])]
                    except ValueError:
                        i += 1
                        continue
                    atoms.append(parts[0])
                    coords.append(xyz)
                i += 1

            steps.append(
                {
                    "type": "neb_step",
                    "energy": energy,
                    "dist": dist_val,
                    "scan_coord": dist_val,
                    "atoms": atoms,
                    "coords": coords,
                }
            )

        return steps

    def parse_basic(self):
        """Parse basic info: SCF Energy, Convergence, Geometry."""
        for i, line in enumerate(self.lines):
            if "Program Version" in line:
                try:
                    self.data["version"] = line.split("Version")[-1].strip().split()[0]
                except (KeyError, IndexError) as e:
                    logging.warning(
                        "Could not parse the ORCA version from %r: %s", line, e
                    )

            line = line.strip()
            uu = line.upper()
            if "FINAL SINGLE POINT ENERGY" in uu:
                try:
                    self.data["scf_energy"] = float(line.split()[-1])
                except (KeyError, IndexError, TypeError, ValueError) as e:
                    logging.warning(
                        "Could not parse the final single-point energy from %r: %s",
                        line,
                        e,
                    )
            if "TOTAL CHARGE" in uu:
                # Could be "Total Charge 0" or "Total Charge ... 0".
                # The phrase also heads a population-analysis column whose
                # value is formatted "0.000000"; int() rejecting that is the
                # intended way to skip it, so it is not worth a warning --
                # the real line elsewhere in the file still sets the value.
                try:
                    parts = line.split()
                    val = int(parts[-1])
                    self.data["charge"] = val
                except (KeyError, IndexError, TypeError, ValueError) as e:
                    logging.debug("Not the charge line: %r (%s)", line.strip(), e)

            if "MULTIPLICITY" in uu:
                try:
                    parts = line.split()
                    val = int(parts[-1])
                    self.data["mult"] = val
                except (KeyError, IndexError, TypeError, ValueError) as e:
                    logging.debug("Not the multiplicity line: %r (%s)", line.strip(), e)
            if (
                "SCF CONVERGED" in uu
                or "OPTIMIZATION CONVERGED" in uu
                or "HURRAY" in uu
            ):
                self.data["converged"] = True

            if "RELAXED SURFACE SCAN" in uu:
                self.data["is_scan"] = True
            if "NUDGED ELASTIC BAND" in uu or " NEB " in uu:
                self.data["is_neb"] = True

            if "CURRENT TRAJECTORY WILL BE WRITTEN TO" in uu:
                # Robust regex extraction
                match = re.search(
                    r"Current trajectory will be written to\s*\.+\s*(.+)",
                    line,
                    re.IGNORECASE,
                )
                if match:
                    self.data["neb_trj_file"] = match.group(1).strip()
                else:
                    # Fallback
                    try:
                        self.data["neb_trj_file"] = line.split()[-1].strip()
                    except (KeyError, IndexError) as e:
                        logging.warning(
                            "Could not parse the NEB trajectory filename from %r: %s",
                            line,
                            e,
                        )

            if "CARTESIAN COORDINATES (ANGSTROEM)" in uu:
                # Read geometry
                self.data["atoms"] = []
                self.data["coords"] = []
                curr = i + 2
                while curr < len(self.lines):
                    l_geo = self.lines[curr].strip()
                    if not l_geo or "---" in l_geo:
                        break
                    xyz = self._parse_xyz_row(l_geo)
                    if xyz is not None:
                        self.data["atoms"].append(xyz[0])
                        self.data["coords"].append(xyz[1])
                    curr += 1
        self.parse_termination_status()

    def parse_termination_status(self):
        """Parse the termination status from the bottom of the file."""
        self.data["termination_status"] = "Running"
        if not self.lines:
            return

        # Check last 100 lines
        last_lines = self.lines[-100:]
        content_block = "\n".join(last_lines)
        uu_block = content_block.upper()

        if "ORCA TERMINATED NORMALLY" in uu_block:
            self.data["termination_status"] = "Terminated normally"
            return

        # Check if there is any error signature
        if (
            "ORCA FINISHED BY ERROR TERMINATION" in uu_block
            or "INPUT ERROR" in uu_block
            or "ERROR !!!" in uu_block
            or "ORCA FINISHED WITH ERROR RETURN" in uu_block
            or re.search(
                r"\[file\s+[^,\]]+,\s*line\s*\d+\]", content_block, re.IGNORECASE
            )
        ):
            self.data["termination_status"] = "ERROR"

    def parse_trajectory(self):
        """Parse Optimization, Scan, and NEB Trajectories."""
        self.data["scan_steps"] = []

        # Look for "RELAXED SURFACE SCAN STEP" or "GEOMETRY OPTIMIZATION CYCLE" or "NEB"
        # And capture Energy + Geometry
        # Actually usually Step header -> Energy -> ... -> Coordinates

        current_scan_step = None

        # Helper to find coords after a header
        def read_block_at(header_idx):
            atoms, coords = [], []
            c_idx = header_idx + 2
            while c_idx < len(self.lines):
                cl = self.lines[c_idx].strip()
                if not cl or "-------" in cl:
                    break
                xyz = self._parse_xyz_row(cl)
                if xyz is not None:
                    atoms.append(xyz[0])
                    coords.append(xyz[1])
                c_idx += 1
            return atoms, coords

        def read_coords_from(idx):
            # First "CARTESIAN COORDINATES (ANGSTROEM)" block within 1000 lines.
            for k in range(idx, min(idx + 1000, len(self.lines))):
                if "CARTESIAN COORDINATES (ANGSTROEM)" in self.lines[k].upper():
                    atoms, coords = read_block_at(k)
                    return atoms, coords, True
            return [], [], False

        def read_last_coords_between(start, end):
            # A relaxed-scan step optimizes its geometry, so the step's
            # structure is the last block before the next step, not the first
            # (which is the unrelaxed starting point).
            for k in range(end - 1, start - 1, -1):
                if "CARTESIAN COORDINATES (ANGSTROEM)" in self.lines[k].upper():
                    atoms, coords = read_block_at(k)
                    return atoms, coords, True
            return read_coords_from(start)

        for i, line in enumerate(self.lines):
            uu_line = line.upper()

            # --- NEB Parsing ---
            if "PATH SUMMARY" in uu_line and i > 0 and "----" in self.lines[i - 1]:
                # Found the summary table
                # Skip header lines
                # Line i: ---------------------------
                # Line i+1:          PATH SUMMARY
                # Line i+2: ---------------------------
                # Line i+3: All forces in Eh/Bohr.
                # Line i+4: Image Dist.(Ang.)    E(Eh)   dE(kcal/mol)  max(|Fp|)  RMS(Fp)

                curr = i + 1
                header_found = False
                while curr < len(self.lines) and curr < i + 10:
                    if "Image" in self.lines[curr] and "E(Eh)" in self.lines[curr]:
                        header_found = True
                        curr += 1
                        break
                    curr += 1

                if header_found:
                    while curr < len(self.lines):
                        l_row = self.lines[curr].strip()
                        if not l_row:
                            break

                        parts = l_row.split()
                        if len(parts) >= 3 and parts[0].isdigit():
                            try:
                                img_idx = int(parts[0])
                                dist = float(parts[1])
                                en = float(parts[2])

                                # We have no geometry, but user said "NO STRUCTURE IS OK"
                                # We provide empty atoms/coords
                                self.data["scan_steps"].append(
                                    {
                                        "type": "neb_image",
                                        "scan_step_id": img_idx,
                                        "step": len(self.data["scan_steps"]) + 1,
                                        "id": img_idx,
                                        "dist": dist,
                                        "energy": en,
                                        "atoms": [],
                                        "coords": [],
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
                                    "Trajectory: could not parse the NEB path-summary row %r: %s",
                                    l_row,
                                    e,
                                )
                        curr += 1

            # Scan Step Header
            if "RELAXED SURFACE SCAN STEP" in uu_line:
                step_idx = 0
                match = re.search(r"STEP\s+(\d+)", uu_line)
                if match:
                    step_idx = int(match.group(1))
                current_scan_step = step_idx

                # Find next step to bound search
                next_marker = len(self.lines)
                for m in range(i + 1, len(self.lines)):
                    if "RELAXED SURFACE SCAN STEP" in self.lines[m].upper():
                        next_marker = m
                        break

                en = None
                conv_info = {}
                coord_val = None
                for k in range(i, next_marker):
                    uu = self.lines[k].strip().upper()
                    if "ACTUAL SCAN COORDINATE" in uu:
                        try:
                            # Format: Actual scan coordinate      ...   1.500000
                            coord_val = float(self.lines[k].split()[-1])
                        except (IndexError, TypeError, ValueError) as e:
                            logging.warning(
                                "Trajectory: could not parse the scan coordinate at step %d from %r: %s",
                                step_idx,
                                self.lines[k],
                                e,
                            )
                    if "FINAL SINGLE POINT ENERGY" in uu:
                        try:
                            en = float(self.lines[k].split()[-1])
                        except (IndexError, TypeError, ValueError) as e:
                            logging.warning(
                                "Trajectory: could not parse the final single-point energy at scan step %d: %s",
                                step_idx,
                                e,
                            )
                    elif "TOTAL ENERGY" in uu and ":" in uu and "EH" in uu:
                        # For ORCA 6: Total Energy       :        -79.79102291629319 Eh
                        try:
                            parts = self.lines[k].split(":")
                            en = float(parts[1].split()[0])
                        except (IndexError, TypeError, ValueError) as e:
                            logging.warning(
                                "Trajectory: could not parse the total energy at scan step %d: %s",
                                step_idx,
                                e,
                            )
                    elif "CURRENT ENERGY" in uu and "...." in uu:
                        # For ORCA relaxation blocks: Current Energy                          ....   -79.800115921 Eh
                        try:
                            parts = self.lines[k].split("....")
                            en = float(parts[1].split()[0])
                        except (IndexError, TypeError, ValueError) as e:
                            logging.warning(
                                "Trajectory: could not parse the current energy at scan step %d: %s",
                                step_idx,
                                e,
                            )
                    elif "GEOMETRY CONVERGENCE" in uu or "CONVERGENCE CRITERIA" in uu:
                        c_idx = k + 1
                        while c_idx < next_marker and c_idx < k + 30:
                            cl = self.lines[c_idx].strip()
                            # Only break on rule if we've already found some data lines
                            # ORCA 6.1.1 has intermediate separators, so we shouldn't break immediately.
                            if "---" in cl:
                                c_idx += 1
                                continue

                            p = cl.split()
                            if len(p) >= 4:
                                s = p[-1]
                                t = p[-2]
                                v = p[-3]
                                n = " ".join(p[:-3]).strip().lower()

                                # Check if it's a standard Yes/No criterion
                                if s.upper() in ["YES", "NO"]:
                                    if n and n != "item":
                                        conv_info[n] = {
                                            "value": v,
                                            "tolerance": t,
                                            "converged": s,
                                        }
                                elif "max(" in cl.lower():
                                    # Parse Max(...) stats
                                    # e.g. Max(Bonds) 0.123 Max(Angles) 0.0
                                    matches = re.findall(
                                        r"(Max\([^)]+\))\s+([-\d\.]+)",
                                        cl,
                                        re.IGNORECASE,
                                    )
                                    for label, val in matches:
                                        conv_info[label] = {
                                            "value": val,
                                            "tolerance": "",
                                            "converged": "INFO",
                                        }
                            c_idx += 1

                # Find matching gradients for this step
                step_grads = []
                candidates = []
                for g_block in self.data.get("all_gradients", []):
                    if g_block["line"] >= i and g_block["line"] < next_marker:
                        candidates.append(g_block["grads"])

                if candidates:
                    step_grads = candidates[-1]
                # Fallback: if not found between markers, maybe it's slightly before the marker?
                # Or just use the one closest to the coordinate block.

                atoms, coords, found = read_last_coords_between(i, next_marker)
                if found:
                    self.data["scan_steps"].append(
                        {
                            "type": "scan_step",
                            "scan_step_id": current_scan_step,
                            "step": step_idx,
                            "energy": en,
                            "scan_coord": coord_val,
                            "atoms": atoms,
                            "coords": coords,
                            "convergence": conv_info,
                            "gradients": step_grads,
                        }
                    )

            elif "OPTIMIZATION CYCLE" in uu_line:
                cycle_idx = 0
                match = re.search(r"CYCLE\s+(\d+)", uu_line)
                if match:
                    cycle_idx = int(match.group(1))

                # Find next cycle to bound search
                next_marker = len(self.lines)
                for m in range(i + 1, len(self.lines)):
                    uu_m = self.lines[m].upper()
                    # Termination markers
                    if "OPTIMIZATION CYCLE" in uu_m:
                        next_marker = m
                        break
                    if (
                        "OPTIMIZATION HAS CONVERGED" in uu_m
                        or "OPTIMIZATION HAS RUN OUT OF CYCLES" in uu_m
                    ):
                        next_marker = m
                        break
                    if "ORCA TERMINATED NORMALLY" in uu_m:
                        next_marker = m
                        break

                en = None
                conv_info = {}

                for k in range(i, next_marker):
                    uu = self.lines[k].strip().upper()
                    if "FINAL SINGLE POINT ENERGY" in uu:
                        try:
                            en = float(self.lines[k].split()[-1])
                        except (IndexError, TypeError, ValueError) as e:
                            logging.warning(
                                "Trajectory: could not parse the final single-point energy at optimization cycle %d: %s",
                                cycle_idx,
                                e,
                            )
                    elif "TOTAL ENERGY" in uu and ":" in uu and "EH" in uu:
                        try:
                            parts = self.lines[k].split(":")
                            en = float(parts[1].split()[0])
                        except (IndexError, TypeError, ValueError) as e:
                            logging.warning(
                                "Trajectory: could not parse the total energy at optimization cycle %d: %s",
                                cycle_idx,
                                e,
                            )
                    elif "CURRENT ENERGY" in uu and "...." in uu:
                        try:
                            parts = self.lines[k].split("....")
                            en = float(parts[1].split()[0])
                        except (IndexError, TypeError, ValueError) as e:
                            logging.warning(
                                "Trajectory: could not parse the current energy at optimization cycle %d: %s",
                                cycle_idx,
                                e,
                            )
                    elif "GEOMETRY CONVERGENCE" in uu or "CONVERGENCE CRITERIA" in uu:
                        c_idx = k + 1
                        while c_idx < next_marker and c_idx < k + 30:
                            cl = self.lines[c_idx].strip()
                            if not cl:
                                c_idx += 1
                                continue
                            if "---" in cl:
                                c_idx += 1
                                continue

                            p = cl.split()
                            if len(p) >= 4:
                                s = p[-1]
                                t = p[-2]
                                v = p[-3]
                                n = " ".join(p[:-3]).strip().lower()

                                # Check if it's a standard Yes/No criterion
                                if s.upper() in ["YES", "NO"]:
                                    if n and n != "item":
                                        conv_info[n] = {
                                            "value": v,
                                            "tolerance": t,
                                            "converged": s,
                                        }
                                elif "max(" in cl.lower():
                                    # Parse Max(...) stats
                                    matches = re.findall(
                                        r"(Max\([^)]+\))\s+([-\d\.]+)",
                                        cl,
                                        re.IGNORECASE,
                                    )
                                    for label, val in matches:
                                        conv_info[label] = {
                                            "value": val,
                                            "tolerance": "",
                                            "converged": "INFO",
                                        }
                            c_idx += 1

                # Find matching gradients for this cycle
                # Gradient block for cycle N is usually printed AFTER the convergence checks of cycle N
                step_grads = []
                # Strategy:
                # 1. Take the LAST gradient block that appeared before next_marker
                # 2. But it MUST be at or after the current cycle index (i)
                candidates = []
                for g_block in self.data.get("all_gradients", []):
                    if g_block["line"] >= i and g_block["line"] < next_marker:
                        candidates.append(g_block["grads"])

                if candidates:
                    step_grads = candidates[
                        -1
                    ]  # Usually only one, but take the last if multiple

                # Special case: if we are at cycle N, and the gradient was printed just BEFORE the header?
                # This doesn't usually happen in ORCA, but for robustness we could check the previous few lines.

                # If we still don't have gradients, look slightly further back?
                # sometimes printed just before? No, usually after.

                atoms, coords, found = read_coords_from(i)
                if found:
                    self.data["scan_steps"].append(
                        {
                            "type": "opt_cycle",
                            "scan_step_id": current_scan_step,
                            "step": cycle_idx,
                            "energy": en,
                            "atoms": atoms,
                            "coords": coords,
                            "convergence": conv_info,
                            "gradients": step_grads,
                        }
                    )

            elif "FINAL ENERGY EVALUATION AT THE STATIONARY POINT" in uu_line:
                # Read energy from this section
                final_en = None
                for k in range(i, min(i + 1500, len(self.lines))):
                    uu_k = self.lines[k].strip().upper()
                    if "FINAL SINGLE POINT ENERGY" in uu_k:
                        try:
                            final_en = float(self.lines[k].split()[-1])
                        except (IndexError, TypeError, ValueError) as e:
                            logging.warning(
                                "Trajectory: could not parse the final stationary-point energy: %s",
                                e,
                            )
                        break
                f_atoms, f_coords, f_found = read_coords_from(i)
                if f_found:
                    last_cycle = max(
                        (
                            s["step"]
                            for s in self.data["scan_steps"]
                            if s.get("type") == "opt_cycle"
                            and s.get("scan_step_id") == current_scan_step
                        ),
                        default=0,
                    )
                    self.data["scan_steps"].append(
                        {
                            "type": "opt_final",
                            "scan_step_id": current_scan_step,
                            "step": last_cycle + 1,
                            "energy": final_en,
                            "atoms": f_atoms,
                            "coords": f_coords,
                            "convergence": {},
                            "gradients": [],
                        }
                    )

    def parse_scan(self):
        """Alias for parse_trajectory."""
        self.parse_trajectory()

    def parse_gradient(self):
        """Alias for parse_gradients."""
        self.parse_gradients()

    def parse_gradients(self):
        """Parse all Cartesian Gradient blocks found in the file."""
        self.data["gradients"] = []  # The last one (default)
        self.data["all_gradients"] = []  # List of {line: int, grads: []}

        gradient_starts = []
        for i, line in enumerate(self.lines):
            stripped = line.strip().upper()
            if "CARTESIAN GRADIENT" in stripped and "NORM" not in stripped:
                gradient_starts.append(i)

        if not gradient_starts:
            return

        for start_idx in gradient_starts:
            block_grads = []
            curr = start_idx + 1
            # Skip header separators or empty lines until data matches format
            found_data = False
            while curr < len(self.lines) and curr < start_idx + 15:
                line = self.lines[curr].strip()
                parts = line.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    found_data = True
                    break
                curr += 1

            if not found_data:
                continue

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "-------" in line or "Difference to" in line:
                    break
                if not line:
                    break

                parts = line.split()
                # Format: 0 C : X Y Z (6 parts) or 1 C X Y Z (5 parts)
                if len(parts) >= 5:
                    try:
                        if parts[2] == ":":
                            if len(parts) >= 6:
                                idx_raw = int(parts[0])
                                sym = parts[1]
                                vx, vy, vz = (
                                    float(parts[3]),
                                    float(parts[4]),
                                    float(parts[5]),
                                )
                                # Usually 1-indexed
                                idx = idx_raw - 1 if idx_raw > 0 else 0
                            else:
                                curr += 1
                                continue
                        else:
                            idx_raw = int(parts[0])
                            sym = parts[1]
                            vx, vy, vz = (
                                float(parts[2]),
                                float(parts[3]),
                                float(parts[4]),
                            )
                            idx = idx_raw - 1 if idx_raw > 0 else 0

                        block_grads.append(
                            {"atom_idx": idx, "atom_sym": sym, "vector": [vx, vy, vz]}
                        )
                    except (
                        AttributeError,
                        KeyError,
                        IndexError,
                        TypeError,
                        ValueError,
                    ) as e:
                        logging.warning(
                            "Gradients: could not parse the gradient row %r: %s",
                            line,
                            e,
                        )
                curr += 1

            if block_grads:
                self.data["all_gradients"].append(
                    {"line": start_idx, "grads": block_grads}
                )

        if self.data["all_gradients"]:
            # Set the last one as the default "gradients"
            self.data["gradients"] = self.data["all_gradients"][-1]["grads"]

    def parse_scan_results_table(self):
        """
        Parses the specific 1D scan summary table from the ORCA output.
        """
        data_start = -1
        # Search from reverse for the summary table
        for i, line in enumerate(reversed(self.lines)):
            if "Actual Energy" in line:
                data_start = len(self.lines) - 1 - i + 1
                break

        if data_start == -1:
            return

        table_vals = []
        for i in range(data_start, len(self.lines)):
            line = self.lines[i].strip()
            # Stop on blank line only after we've collected some data
            if not line:
                if table_vals:
                    break
                continue

            parts = line.split()
            if len(parts) >= 2:
                try:
                    coord = float(parts[0])
                    en = float(parts[-1])
                    table_vals.append({"coord": coord, "energy": en})
                except (IndexError, TypeError, ValueError):
                    # Non-numeric line after data = end of table
                    if table_vals:
                        break

        if not table_vals:
            return

        if "scan_steps" not in self.data:
            self.data["scan_steps"] = []

        steps = self.data["scan_steps"]
        if not steps:
            # Reconstruct entire trajectory from the summary table
            for idx, v in enumerate(table_vals):
                steps.append(
                    {
                        "type": "scan_step_summary",
                        "scan_step_id": idx,
                        "step": idx,
                        "energy": v["energy"],
                        "scan_coord": v["coord"],
                        "atoms": [],
                        "coords": [],
                    }
                )
        else:
            # Map coordinates to existing optimization steps
            sids = [
                s.get("scan_step_id")
                for s in steps
                if s.get("scan_step_id") is not None
            ]
            if sids:
                offset = min(sids)
                for s in steps:
                    sid = s.get("scan_step_id", None)
                    if sid is not None:
                        idx = sid - offset
                        if 0 <= idx < len(table_vals):
                            s["scan_coord"] = table_vals[idx]["coord"]
            else:
                for idx, s in enumerate(steps[: len(table_vals)]):
                    s["scan_coord"] = table_vals[idx]["coord"]

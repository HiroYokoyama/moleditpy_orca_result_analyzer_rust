"""Parsing mixin for dipole, atomic charges, NBO, Mayer bonds and energy components."""

import math
import re
import logging


#: 1 atomic unit of electric dipole moment (e*a0) in Debye.
AU_TO_DEBYE = 2.541746473


class _PropertyParsingMixin:
    """Dipole, charge, NBO, Mayer and energy-component parsing for OrcaParser."""

    def parse_dipole(self) -> None:
        """Extract the total dipole moment (and per-origin dipoles) into self.data."""
        # Look for "Total Dipole Moment"
        self.data["dipoles"] = None
        self.data["dipole"] = None

        candidates = []
        for i, line in enumerate(self.lines):
            uu = line.upper()
            if "TOTAL DIPOLE MOMENT" in uu and ":" in line:
                candidates.append(i)

        if not candidates:
            return

        # Take the last occurrence
        idx = candidates[-1]
        line = self.lines[idx]
        parts = line.split(":")
        if len(parts) > 1:
            try:
                vec_str = parts[1].strip().split()
                if len(vec_str) >= 3:
                    x, y, z = float(vec_str[0]), float(vec_str[1]), float(vec_str[2])

                    mag_au = None
                    mag_debye = None
                    for k in range(idx + 1, min(idx + 7, len(self.lines))):
                        lk = self.lines[k]
                        if "Magnitude" in lk and ":" in lk:
                            try:
                                val = float(lk.split(":")[1].strip())
                                if "Debye" in lk:
                                    mag_debye = val
                                else:
                                    mag_au = val
                            except (IndexError, TypeError, ValueError):
                                logging.debug(
                                    "Unparseable dipole magnitude line: %r",
                                    lk,
                                    exc_info=True,
                                )
                    # ORCA prints the components in atomic units (e*a0) and
                    # the magnitudes separately in both a.u. and Debye. Keep
                    # them apart: reporting the a.u. vector as Debye put the
                    # components and the magnitude on different scales.
                    if mag_au is None:
                        mag_au = math.sqrt(x * x + y * y + z * z)
                    if mag_debye is None:
                        mag_debye = mag_au * AU_TO_DEBYE

                    # Prefer ORCA's own pair of magnitudes for the scale factor
                    # so the vector stays consistent with what it printed.
                    if mag_au > 1e-12:
                        to_debye = mag_debye / mag_au
                    else:
                        to_debye = AU_TO_DEBYE

                    self.data["dipoles"] = {
                        "vector_au": (x, y, z),
                        "vector_debye": (x * to_debye, y * to_debye, z * to_debye),
                        "magnitude_au": mag_au,
                        "magnitude_debye": mag_debye,
                        # Back-compat aliases (Debye, matching the UI labels).
                        "vector": (x * to_debye, y * to_debye, z * to_debye),
                        "magnitude": mag_debye,
                    }
                    self.data["dipole"] = self.data["dipoles"]
            except (AttributeError, KeyError, IndexError, TypeError, ValueError) as e:
                logging.warning(
                    "Dipole moment: could not parse the dipole vector block: %s", e
                )

    def parse_spin_contamination(self):
        """Extract the UHF/UKS spin expectation value <S**2>.

        ORCA prints, for open-shell calculations:
            Expectation value of <S**2>     :     2.007028
            Ideal value S*(S+1) for S=1.0   :     2.000000

        Stores {actual, ideal, contamination} (contamination = actual - ideal),
        or None for restricted/closed-shell jobs that print no such block. The
        last occurrence wins (final geometry in an optimization).
        """
        self.data["spin_s2"] = None
        actual = None
        ideal = None
        for line in self.lines:
            if "Expectation value of <S**2>" in line and ":" in line:
                try:
                    actual = float(line.split(":")[1].strip().split()[0])
                except (IndexError, TypeError, ValueError) as e:
                    logging.warning(
                        "Spin contamination: could not parse the actual <S**2> value from %r: %s",
                        line,
                        e,
                    )
            elif "Ideal value" in line and "S*(S+1)" in line and ":" in line:
                try:
                    ideal = float(line.split(":")[1].strip().split()[0])
                except (IndexError, TypeError, ValueError) as e:
                    logging.warning(
                        "Spin contamination: could not parse the ideal S*(S+1) value from %r: %s",
                        line,
                        e,
                    )
        if actual is not None:
            self.data["spin_s2"] = {
                "actual": actual,
                "ideal": ideal,
                "contamination": (actual - ideal) if ideal is not None else None,
            }

    def parse_dispersion(self):
        """Extract the London dispersion correction energy (DFT-D3/D4), if present.

        ORCA prints e.g.:  "Dispersion correction           -0.000882087"
        Stored as a float (Eh) or None. The case-sensitive match avoids the
        prose line "... London dispersion correction". Last occurrence wins.
        """
        self.data["dispersion"] = None
        for line in self.lines:
            m = re.search(
                r"Dispersion correction\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)", line
            )
            if m:
                try:
                    self.data["dispersion"] = float(m.group(1))
                except (KeyError, IndexError, TypeError, ValueError) as e:
                    logging.warning(
                        "Dispersion correction: could not parse value from %r: %s",
                        m.group(1),
                        e,
                    )

    def parse_energy_components(self):
        """Parse post-HF correlation energy components (MP2 / CCSD(T) / ...).

        Captures whichever labels are present, e.g.:
            MP2 CORRELATION ENERGY   :   -0.201139520 Eh
            MP2 TOTAL ENERGY:            -76.162123498 Eh
            E(0)                  ...    -75.960983979
            E(CORR)(corrected)    ...     -0.210803526
            Triples Correction (T)...     -0.002885334
            Final correlation energy ...  -0.213688859
            E(CCSD)               ...    -76.171787504
            E(CCSD(T))            ...    -76.174672838
            T1 diagnostic         ...      0.005944123
        Stored (in this order) as a list of {label, value, dimensionless}.
        Matching is case-sensitive so DLPNO's "SL-MP2 correlation energy"
        intermediate lines do not collide with canonical "MP2 CORRELATION
        ENERGY". The last occurrence of each label wins.
        """
        self.data["energy_components"] = []
        specs = [
            ("E(0)", "Reference energy E(0)", False),
            ("MP2 CORRELATION ENERGY", "MP2 correlation energy", False),
            ("MP2 TOTAL ENERGY", "MP2 total energy", False),
            ("E(CORR)(corrected)", "CCSD correlation energy", False),
            ("Triples Correction (T)", "(T) triples correction", False),
            ("Final correlation energy", "Total correlation energy", False),
            ("E(CCSD)", "E(CCSD)", False),
            ("E(CCSD(T))", "E(CCSD(T))", False),
            ("T1 diagnostic", "T1 diagnostic", True),
        ]
        float_re = re.compile(r"[-+]?\d+\.\d+")
        found = {}
        for line in self.lines:
            for sub, label, dimensionless in specs:
                if sub in line:
                    nums = float_re.findall(line)
                    if nums:
                        try:
                            found[label] = (float(nums[-1]), dimensionless)
                        except ValueError:
                            logging.debug(
                                "Unparseable energy component value in line: %r",
                                line,
                                exc_info=True,
                            )
        self.data["energy_components"] = [
            {"label": label, "value": found[label][0], "dimensionless": found[label][1]}
            for _sub, label, _dim in specs
            if label in found
        ]

    def parse_mayer_bond_orders(self):
        """Parse the Mayer bond-order matrix.

        ORCA prints, after the Mayer population analysis:
            Mayer bond orders larger than 0.100000
            B(  0-C ,  1-C ) :   1.3925 B(  0-C ,  5-C ) :   1.3924 ...
            B(  1-C ,  2-C ) :   1.3927 ...

        Stored as a list of {atom_idx1, atom_sym1, atom_idx2, atom_sym2, order}
        with atom_idx1 < atom_idx2 (0-based), or [] if absent. The last block
        wins (final geometry in an optimization).
        """
        self.data["mayer_bond_orders"] = []

        start = -1
        for i, line in enumerate(self.lines):
            if "Mayer bond orders" in line:
                start = i

        if start == -1:
            return

        pattern = re.compile(
            r"B\(\s*(\d+)-(\w+)\s*,\s*(\d+)-(\w+)\s*\)\s*:\s*(-?\d+\.\d+)"
        )

        bonds = []
        for j in range(start + 1, len(self.lines)):
            line = self.lines[j]
            if not line.strip():
                break
            matches = pattern.findall(line)
            if not matches:
                if bonds:
                    break
                continue
            for idx1, sym1, idx2, sym2, order in matches:
                a1, s1, a2, s2 = int(idx1), sym1, int(idx2), sym2
                if a1 > a2:
                    a1, s1, a2, s2 = a2, s2, a1, s1
                bonds.append(
                    {
                        "atom_idx1": a1,
                        "atom_sym1": s1,
                        "atom_idx2": a2,
                        "atom_sym2": s2,
                        "order": float(order),
                    }
                )

        self.data["mayer_bond_orders"] = bonds

    def parse_nbo_orbitals(self):
        """Parse the "NATURAL BOND ORBITALS (Summary)" table (NBO analysis).

        One NBO per line, e.g.:
            1. CR ( 1) O  1             1.99996   -18.75648
            4. BD ( 1) O  1- H  2       1.99864    -0.60737  23(v)
            6. BD*( 1) O  1- H  2       0.00013     0.38562
        Stored as a list of {index, type, atoms, occupancy, energy}, or [].
        Types: CR (core), LP (lone pair), BD (bond), BD* (antibond),
        RY/RY* (Rydberg), LV (lone vacancy).

        If the section appears more than once (e.g. multiple NBO analyses
        across an optimization/scan job), the last occurrence wins.
        """
        self.data["nbo_orbitals"] = []
        start = -1
        for i, line in enumerate(self.lines):
            if "NATURAL BOND ORBITALS (Summary)" in line:
                start = i
        if start == -1:
            return
        pattern = re.compile(
            r"^\s*(\d+)\.\s+(\S+?)\s*\(\s*\d+\)\s+(.+?)\s{2,}(\d+\.\d+)\s+(-?\d+\.\d+)"
        )
        orbitals = []
        for j in range(start + 1, min(start + 2000, len(self.lines))):
            line = self.lines[j]
            if "NBO analysis completed" in line:
                break
            m = pattern.match(line)
            if m:
                atoms_str = m.group(3).strip()
                atom_indices = [int(n) - 1 for n in re.findall(r"\d+", atoms_str)]
                orbitals.append(
                    {
                        "index": int(m.group(1)),
                        "type": m.group(2),
                        "atoms": atoms_str,
                        "atom_indices": atom_indices,
                        "occupancy": float(m.group(4)),
                        "energy": float(m.group(5)),
                    }
                )
        self.data["nbo_orbitals"] = orbitals

    def parse_nbo_perturbation(self):
        """Parse NBO second-order perturbation theory donor->acceptor analysis.

        Lines under "SECOND ORDER PERTURBATION THEORY ANALYSIS ...":
            2. LP ( 1) O  1            18. RY ( 2) H  2            1.95    1.87   0.054
        giving donor NBO, acceptor NBO, E(2) [kcal/mol], E(NL)-E(L) [a.u.],
        and F(L,NL) [a.u.]. Stored as a list of
        {donor, acceptor, e2_kcal, e_diff, fock}, or [].

        If the section appears more than once (e.g. multiple NBO analyses
        across an optimization/scan job), the last occurrence wins.
        """
        self.data["nbo_perturbation"] = []
        start = -1
        for i, line in enumerate(self.lines):
            if "SECOND ORDER PERTURBATION THEORY ANALYSIS" in line.upper():
                start = i
        if start == -1:
            return
        pattern = re.compile(
            r"^\s*(\d+)\.\s+(\S.*?\S)\s+(\d+)\.\s+(\S.*?\S)\s+"
            r"(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
        )
        inter = []
        for j in range(start + 1, len(self.lines)):
            line = self.lines[j]
            if (
                "NATURAL BOND ORBITALS (Summary)" in line
                or "NBO analysis completed" in line
            ):
                break
            m = pattern.match(line)
            if m:
                inter.append(
                    {
                        "donor": f"{m.group(1)}. {m.group(2)}",
                        "acceptor": f"{m.group(3)}. {m.group(4)}",
                        "e2_kcal": float(m.group(5)),
                        "e_diff": float(m.group(6)),
                        "fock": float(m.group(7)),
                    }
                )
        self.data["nbo_perturbation"] = inter

    def parse_nbo_hybrids(self):
        """Augment nbo_orbitals with per-atom hybridization (%s / %p / %d).

        From the detailed "(Occupancy) Bond orbital / Coefficients / Hybrids"
        block, e.g.:
            1. (1.99996) CR ( 1) O  1            s(100.00%)
            2. (1.99698) LP ( 1) O  1            s(  0.00%)p 1.00( 99.92%) ...
            4. (1.99864) BD ( 1) O  1- H  2
                       ( 72.35%)   0.8506* O  1 s( 22.56%)p 3.43( 77.29%) ...
                       ( 27.65%)   0.5258* H  2 s( 99.76%)p 0.00(  0.24%)
        Adds a "hybrids" list to each matching NBO: one component per atom with
        {atom_sym, atom_idx, s_pct, p_pct, d_pct, weight_pct?, label}.

        If the section appears more than once, the last occurrence wins (it
        corresponds to the same NBO analysis run as nbo_orbitals/
        nbo_perturbation, which both already pick their last occurrence).
        """
        orbitals = self.data.get("nbo_orbitals", [])
        if not orbitals:
            return

        start = -1
        for i, line in enumerate(self.lines):
            if "Bond orbital / Coefficients / Hybrids" in line:
                start = i
        if start == -1:
            return

        s_re = re.compile(r"s\(\s*([\d.]+)%\)")
        p_re = re.compile(r"p\s*[\d.]*\(\s*([\d.]+)%\)")
        d_re = re.compile(r"d\s*[\d.]*\(\s*([\d.]+)%\)")
        header_re = re.compile(
            r"^\s*(\d+)\.\s*\(\s*[\d.]+\)\s*\S+?\s*\(\s*\d+\)\s*(.*)$"
        )
        pol_re = re.compile(
            r"^\s*\(\s*([\d.]+)%\)\s+[-\d.]+\*?\s*([A-Za-z]+)\s+(\d+)\s+(.*)$"
        )

        def _label(s_pct, p_pct, d_pct):
            # Below ~5% s the sp^n ratio explodes and is meaningless (typical of
            # diffuse Rydberg orbitals); fall back to the dominant character.
            if s_pct < 5.0:
                return "d" if d_pct > p_pct else "p"
            if p_pct <= 1.0 and d_pct <= 1.0:
                return "s"
            base = f"sp{p_pct / s_pct:.2f}" if p_pct > 1.0 else "s"
            # Annotate d only when chemically meaningful (> ~1%).
            if d_pct > 1.0:
                base += f"d{d_pct / s_pct:.2f}"
            return base

        def _extract(text, atom_sym=None, atom_num=None, weight=None):
            sm = s_re.search(text)
            if not sm:
                return None
            pm = p_re.search(text)
            dm = d_re.search(text)
            s_pct = float(sm.group(1))
            p_pct = float(pm.group(1)) if pm else 0.0
            d_pct = float(dm.group(1)) if dm else 0.0
            comp = {
                "atom_sym": atom_sym,
                "atom_idx": (atom_num - 1) if atom_num is not None else None,
                "s_pct": s_pct,
                "p_pct": p_pct,
                "d_pct": d_pct,
                "label": _label(s_pct, p_pct, d_pct),
                "raw": text[sm.start() :].rstrip(),
            }
            if weight is not None:
                comp["weight_pct"] = weight
            return comp

        hybrids = {}  # nbo index -> list of components
        current_idx = None
        end_markers = (
            "SECOND ORDER PERTURBATION",
            "NATURAL BOND ORBITALS (Summary)",
            "NBO analysis completed",
        )
        for j in range(start + 1, len(self.lines)):
            line = self.lines[j]
            if any(mk in line for mk in end_markers):
                break

            hm = header_re.match(line)
            if hm:
                current_idx = int(hm.group(1))
                rest = hm.group(2)
                comp = _extract(rest)
                if comp is not None:
                    am = re.match(r"\s*([A-Za-z]+)\s+(\d+)", rest)
                    if am:
                        comp["atom_sym"] = am.group(1)
                        comp["atom_idx"] = int(am.group(2)) - 1
                    hybrids.setdefault(current_idx, []).append(comp)
                continue

            pm = pol_re.match(line)
            if pm and current_idx is not None:
                comp = _extract(
                    pm.group(4), pm.group(2), int(pm.group(3)), float(pm.group(1))
                )
                if comp is not None:
                    hybrids.setdefault(current_idx, []).append(comp)

        for orb in orbitals:
            orb["hybrids"] = hybrids.get(orb["index"], [])

    def parse_charges(self):
        """Extract every available atomic-charge scheme into self.data["charges"]."""
        self.data["charges"] = {}  # type -> list of {atom_idx, atom_sym, charge}

        # Section Markers
        mulliken_start = -1
        loewdin_start = -1
        hirshfeld_start = -1
        mayer_start = -1
        nbo_start = -1
        chelpg_start = -1
        mbis_start = -1
        resp_start = -1
        fmo_start = -1

        for i, line in enumerate(self.lines):
            uu = line.upper()
            if "MULLIKEN ATOMIC CHARGES" in uu:
                mulliken_start = i
            elif "LOEWDIN ATOMIC CHARGES" in uu:
                loewdin_start = i
            elif "HIRSHFELD ANALYSIS" in uu:
                hirshfeld_start = i
            elif "MAYER POPULATION ANALYSIS" in uu:
                mayer_start = i
            elif "NATURAL POPULATIONS" in uu:
                nbo_start = i
            elif (
                "CHELPG CHARGES" in uu
                and "GENERATION" not in uu
                and "CALCULATED" not in uu
            ):
                # Real ORCA prints the table header as "CHELPG Charges"
                # (not "CHELPG ATOMIC CHARGES"); exclude the surrounding
                # "CHELPG CHARGES GENERATION" / "CHELPG charges calculated"
                # status lines that also contain the phrase.
                chelpg_start = i
            elif "MBIS ANALYSIS" in uu:
                mbis_start = i
            elif (
                "RESP CHARGES" in uu
                and "GENERATION" not in uu
                and "CALCULATED" not in uu
            ):
                # Real ORCA 6 prints the table header as "RESP Charges"
                # (not "RESP ATOMIC CHARGES"); exclude the surrounding
                # "RESP CHARGES GENERATION" / "RESP charges calculated" lines.
                resp_start = i
            elif "FRONTIER MOLECULAR ORBITAL POPULATION ANALYSIS" in uu:
                fmo_start = i

        def parse_standard_block(
            start_idx, header_lines=2, hirshfeld=False, mbis=False
        ):
            res = []
            if start_idx == -1:
                return res
            curr = start_idx + header_lines
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line or "---" in line or "Sum of" in line:
                    if res:
                        break
                    curr += 1
                    continue

                parts = line.split()
                if len(parts) >= 3:
                    try:
                        # Handle merged "0C" case
                        first_part = parts[0]
                        match = re.match(r"^(\d+)([a-zA-Z]+)$", first_part)
                        if match:
                            idx_str = match.group(1)
                            sym = match.group(2)
                        else:
                            idx_str = first_part.strip(":")
                            # If "0 C", sym is parts[1]
                            sym = parts[1].strip(":")

                        if hirshfeld or mbis:
                            val = float(parts[2])  # 0 C Charge Spin
                        else:
                            # 0 C : -1.23 or 0 C -1.23
                            val = (
                                float(parts[3]) if parts[2] == ":" else float(parts[2])
                            )

                        atom_data = {
                            "atom_idx": int(idx_str),
                            "atom_sym": sym,
                            "charge": val,
                        }

                        # Capture spin if available (Hirshfeld/MBIS)
                        if hirshfeld:
                            if len(parts) >= 4:
                                atom_data["spin"] = float(parts[3])
                        elif mbis:
                            if len(parts) >= 4:
                                # ATOM CHARGE POPULATION [SPIN]
                                atom_data["population"] = float(parts[3])
                            if len(parts) >= 5:
                                atom_data["spin"] = float(parts[4])

                        res.append(atom_data)
                    except (
                        AttributeError,
                        KeyError,
                        IndexError,
                        TypeError,
                        ValueError,
                    ) as e:
                        # The block is scanned line by line and the numeric
                        # conversion is what separates data rows from the
                        # headers and blank lines around them, so failing
                        # here is the normal case, not a problem to report.
                        logging.debug("Not a charge row: %r (%s)", line.strip(), e)
                curr += 1
            return res

        self.data["charges"]["Mulliken"] = parse_standard_block(mulliken_start)
        self.data["charges"]["Loewdin"] = parse_standard_block(loewdin_start)
        self.data["charges"]["Hirshfeld"] = parse_standard_block(
            hirshfeld_start, hirshfeld=True
        )
        self.data["charges"]["CHELPG"] = parse_standard_block(chelpg_start)
        self.data["charges"]["MBIS"] = parse_standard_block(
            mbis_start, header_lines=3, mbis=True
        )
        self.data["charges"]["RESP"] = parse_standard_block(resp_start)

        # Clean up empty entries
        self.data["charges"] = {k: v for k, v in self.data["charges"].items() if v}

        # Mayer Parsing (QA is Mulliken charge)
        if mayer_start != -1:
            mayer_res = []
            curr = mayer_start + 1
            while curr < len(self.lines) and curr < mayer_start + 15:
                if "ATOM" in self.lines[curr] and "QA" in self.lines[curr]:
                    curr += 1
                    break
                curr += 1

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line or "---" in line or "Mayer bond" in line:
                    break
                parts = line.split()
                # QA lives at parts[4], so a >= 4 guard let a short row through
                # to an IndexError that the blanket handler below swallowed --
                # the row just vanished from the table with no trace.
                if len(parts) >= 5:
                    try:
                        idx = int(parts[0])
                        sym = parts[1]
                        # NA, ZA, QA, VA, BVA, FA
                        # parts[2] = NA
                        # parts[3] = ZA (Atomic Number)
                        # parts[4] = QA (Charge)
                        qa = float(parts[4])

                        # Extra Mayer Indices
                        extra = {}
                        if len(parts) >= 8:
                            extra["valency"] = float(parts[5])  # VA
                            extra["bonded_valency"] = float(parts[6])  # BVA
                            extra["free_valency"] = float(parts[7])  # FA

                        mayer_res.append(
                            {"atom_idx": idx, "atom_sym": sym, "charge": qa, **extra}
                        )
                    except (
                        AttributeError,
                        KeyError,
                        IndexError,
                        TypeError,
                        ValueError,
                    ) as e:
                        # Log the row, not `sym`: if int(parts[0]) fails on the
                        # first row, `sym` is unbound and the handler itself
                        # raised UnboundLocalError out of parse_all.
                        logging.warning(
                            "Mayer charges: could not parse the valency row %r: %s",
                            line,
                            e,
                        )
                curr += 1
            if mayer_res:
                self.data["charges"]["Mayer"] = mayer_res
                if not self.data["charges"].get("Mulliken", None):
                    self.data["charges"]["Mulliken"] = mayer_res

        # NBO Parsing
        if nbo_start != -1:
            nbo_charges = []

            # 1. Try to find Detailed Summary Table first
            # "Summary of Natural Population Analysis"
            summary_start = -1
            # Search a bit deeper than nbo_start, as "NATURAL POPULATIONS" might be just one header
            # But usually the summary is a distinct block.
            # Let's search the whole file? No, usually near the end or after NBO analysis
            # Optimization: Look forward from nbo_start

            for i in range(nbo_start, min(nbo_start + 2000, len(self.lines))):
                if "Summary of Natural Population Analysis" in self.lines[i]:
                    summary_start = i
                    break

            if summary_start != -1:
                # Parse the detailed table
                curr = summary_start + 1
                while curr < len(self.lines):
                    line = self.lines[curr].strip()
                    if "Atom No" in line and "Charge" in line and "Core" in line:
                        curr += 1  # Skip header line
                        if curr < len(self.lines) and "----" in self.lines[curr]:
                            curr += 1  # Skip separator
                        break
                    curr += 1

                while curr < len(self.lines):
                    line = self.lines[curr].strip()
                    if "====" in line or "Total" in line:
                        break
                    if not line:
                        curr += 1
                        continue

                    parts = line.split()
                    # C  1   -0.49495      1.99999     4.48487    0.01009     6.49495
                    if len(parts) >= 7:
                        try:
                            sym = parts[0]
                            idx = int(parts[1])
                            chg = float(parts[2])
                            core = float(parts[3])
                            val = float(parts[4])
                            ryd = float(parts[5])
                            tot = float(parts[6])

                            nbo_charges.append(
                                {
                                    "atom_idx": idx - 1,
                                    "atom_sym": sym,
                                    "charge": chg,
                                    "core": core,
                                    "valence": val,
                                    "rydberg": ryd,
                                    "total": tot,
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
                                "NBO charges: could not parse the summary row %r: %s",
                                line,
                                e,
                            )
                    curr += 1

            # 2. Fallback if no summary found (or parsing failed), try simple block near nbo_start
            if not nbo_charges:
                curr = nbo_start + 1
                while curr < len(self.lines):
                    line = self.lines[curr].strip()
                    if "---" in line:
                        curr += 1
                        continue
                    if (
                        "================" in line
                        or "Natural Electron Configuration" in line
                    ):
                        if nbo_charges:
                            break
                        curr += 1
                        continue
                    if not line:
                        if nbo_charges:
                            break  # Stop on empty line if we have data
                        curr += 1
                        continue

                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            # Simple format: Atom No Charge ...
                            # or format: C 1 -0.123
                            sym = parts[0]
                            idx = int(parts[1])
                            chg = float(parts[2])
                            nbo_charges.append(
                                {"atom_idx": idx - 1, "atom_sym": sym, "charge": chg}
                            )
                        except (IndexError, TypeError, ValueError) as e:
                            logging.warning(
                                "NBO charges: could not parse the fallback-format row %r: %s",
                                line,
                                e,
                            )
                    curr += 1

            if nbo_charges:
                self.data["charges"]["NBO"] = nbo_charges

        # FMO Parsing
        if fmo_start != -1:
            fmo_data = []
            curr = fmo_start + 1

            table_start = False
            # Look for table header
            while curr < len(self.lines) and curr < fmo_start + 40:
                if "--------" in self.lines[curr]:
                    # Check previous lines for "Atom" "HOMO" etc.
                    # Usually two or three lines of header, then barrier
                    # Or "Atom   Q(Mulliken) ..."
                    if curr > 0 and "Atom" in self.lines[curr - 1]:
                        table_start = True
                        curr += 1
                        break
                    if curr > 1 and "HOMO" in self.lines[curr - 2]:
                        table_start = True
                        curr += 1
                        break
                curr += 1

            if table_start:
                while curr < len(self.lines):
                    line = self.lines[curr].strip()
                    if not line or "--------" in line:
                        if fmo_data:
                            break
                        curr += 1
                        continue

                    parts = line.split()
                    # 0-C 0.937 0.906 0.804 0.755
                    if len(parts) >= 5:
                        try:
                            atom_lbl = parts[0]  # 0-C
                            if "-" in atom_lbl:
                                # Split only on first hyphen to handle negative indices or strange names?
                                # Usually ORCA format is 0-C, 1-H etc.
                                p_lbl = atom_lbl.split("-")
                                idx_str = p_lbl[0]
                                sym = p_lbl[1]
                                idx = int(idx_str)
                            else:
                                # Fallback if just C or just 0
                                idx = len(fmo_data)
                                sym = atom_lbl

                            homo_m = float(parts[1])
                            homo_l = float(parts[2])
                            lumo_m = float(parts[3])
                            lumo_l = float(parts[4])

                            fmo_data.append(
                                {
                                    "atom_idx": idx,
                                    "atom_sym": sym,
                                    # Use Mulliken HOMO as primary visual if asked for "charge"
                                    "charge": homo_m,
                                    "homo_mulliken": homo_m,
                                    "homo_loewdin": homo_l,
                                    "lumo_mulliken": lumo_m,
                                    "lumo_loewdin": lumo_l,
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
                                "FMO charges: could not parse the HOMO/LUMO row %r: %s",
                                line,
                                e,
                            )
                    curr += 1

                if fmo_data:
                    self.data["charges"]["FMO"] = fmo_data

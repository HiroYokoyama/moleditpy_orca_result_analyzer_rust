// orca_parser_rs — Rust-powered parser for ORCA quantum chemistry output files.
// Exposed to Python via PyO3.  Drop-in replacement for parser.py's OrcaParser.
//
// Public surface:
//   parse_all(content: str) -> dict          — parse a complete ORCA output string
//   parse_xyz_content(content: str) -> list  — parse multi-frame XYZ content

use once_cell::sync::Lazy;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use regex::Regex;
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Compiled regex cache (compiled once, reused for every call)
// ---------------------------------------------------------------------------

static RE_FLOAT: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+\.?").unwrap());
static RE_ENERGY_LABEL: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r"(?i)(?:Energy|E)[=:\s]+([-+]?\d*\.\d+|[-+]?\d+\.?)").unwrap()
});
static RE_DIST_LABEL: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r"(?i)(?:Dist(?:ance)?|Coord(?:inate)?|Scan)[=:\s]+([-+]?\d*\.\d+|[-+]?\d+\.?)")
        .unwrap()
});
static RE_STATE_HEADER: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"(?i)STATE\s+(\d+)\s*:").unwrap());
static RE_ENERGY_EV: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"([-+]?\d*\.\d+)\s*eV").unwrap());
static RE_ENERGY_CM: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"([-+]?\d*\.\d+)\s*cm").unwrap());
static RE_ENERGY_NM: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"([-+]?\d*\.\d+)\s*nm").unwrap());
static RE_TARGET_STATE: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^(\d+)").unwrap());
static RE_CYCLE: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)CYCLE\s+(\d+)").unwrap());
static RE_SCAN_STEP: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)STEP\s+(\d+)").unwrap());
static RE_TEMP: Lazy<Regex> = Lazy::new(|| Regex::new(r"(\d+\.\d+)\s*K").unwrap());
static RE_MERGED_ATOM: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^(\d+)([A-Za-z]+)$").unwrap());
static RE_NEB_TRJ: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r"(?i)Current trajectory will be written to\s*\.+\s*(.+)").unwrap()
});
static RE_MAX_STATS: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"(?i)(Max\([^)]+\))\s+([-\d\.]+)").unwrap());
static RE_FLOATS_IN_LINE: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?").unwrap()
});

// ---------------------------------------------------------------------------
// Tiny helpers
// ---------------------------------------------------------------------------

fn parse_f64(s: &str) -> Option<f64> {
    s.trim().parse::<f64>().ok()
}

fn parse_i32(s: &str) -> Option<i32> {
    s.trim().parse::<i32>().ok()
}

fn last_float(line: &str) -> Option<f64> {
    RE_FLOATS_IN_LINE
        .find_iter(line)
        .last()
        .and_then(|m| parse_f64(m.as_str()))
}

fn all_floats(line: &str) -> Vec<f64> {
    RE_FLOATS_IN_LINE
        .find_iter(line)
        .filter_map(|m| parse_f64(m.as_str()))
        .collect()
}

// ---------------------------------------------------------------------------
// Internal data structures
// ---------------------------------------------------------------------------

#[derive(Default)]
struct ParsedData {
    scf_energy: Option<f64>,
    converged: bool,
    atoms: Vec<String>,
    coords: Vec<[f64; 3]>,
    charge: i32,
    mult: i32,
    version: Option<String>,
    is_scan: bool,
    is_neb: bool,
    neb_trj_file: Option<String>,

    scf_traces: Vec<ScfTrace>,
    frequencies: Vec<FreqEntry>,
    orbital_energies: Vec<OrbitalEnergy>,
    charges: HashMap<String, Vec<ChargeEntry>>,
    dipole: Option<Dipole>,
    tddft: Vec<TddftState>,
    nmr_shielding: Vec<NmrShielding>,
    nmr_couplings: Vec<NmrCoupling>,
    thermal: HashMap<String, f64>,
    scan_steps: Vec<ScanStep>,
    all_gradients: Vec<GradientBlock>,
    gradients: Vec<Gradient>,
    basis_set_shells: Vec<BasisShell>,
    mo_coeffs: HashMap<String, MoCoeff>,
}

#[derive(Default, Clone)]
struct ScfTrace {
    step: String,
    iterations: Vec<ScfIter>,
}

#[derive(Default, Clone)]
struct ScfIter {
    iter: i32,
    energy: f64,
}

#[derive(Default, Clone)]
struct FreqEntry {
    freq: f64,
    ir: f64,
    raman: f64,
    vector: Vec<[f64; 3]>,
}

#[derive(Default, Clone)]
struct OrbitalEnergy {
    index: i32,
    occupation: f64,
    energy_eh: f64,
    energy_ev: f64,
    spin: String,
    orb_type: String,
}

#[derive(Default, Clone)]
struct ChargeEntry {
    atom_idx: i32,
    atom_sym: String,
    charge: f64,
    spin: Option<f64>,
    population: Option<f64>,
    valency: Option<f64>,
    bonded_valency: Option<f64>,
    free_valency: Option<f64>,
    core: Option<f64>,
    valence: Option<f64>,
    rydberg: Option<f64>,
    total: Option<f64>,
    homo_mulliken: Option<f64>,
    homo_loewdin: Option<f64>,
    lumo_mulliken: Option<f64>,
    lumo_loewdin: Option<f64>,
}

#[derive(Default, Clone)]
struct Dipole {
    x: f64,
    y: f64,
    z: f64,
    magnitude: f64,
}

#[derive(Default, Clone)]
struct TddftState {
    state: i32,
    energy_ev: f64,
    energy_nm: f64,
    energy_cm: f64,
    osc: f64,
    osc_len: f64,
    osc_vel: f64,
    rotatory_strength: f64,
    rot_len: f64,
    rot_vel: f64,
    transitions: Vec<String>,
}

#[derive(Default, Clone)]
struct NmrShielding {
    atom_idx: i32,
    atom_sym: String,
    shielding: f64,
}

#[derive(Default, Clone)]
struct NmrCoupling {
    atom_idx1: i32,
    atom_idx2: i32,
    coupling: f64,
}

#[derive(Default, Clone)]
struct ScanStep {
    step_type: String,
    scan_step_id: Option<i32>,
    step: i32,
    energy: f64,
    scan_coord: Option<f64>,
    dist: Option<f64>,
    atoms: Vec<String>,
    coords: Vec<[f64; 3]>,
    convergence: HashMap<String, ConvEntry>,
    gradients: Vec<Gradient>,
}

#[derive(Default, Clone)]
struct ConvEntry {
    value: String,
    tolerance: String,
    converged: String,
}

#[derive(Default, Clone)]
struct GradientBlock {
    line: usize,
    grads: Vec<Gradient>,
}

#[derive(Default, Clone)]
struct Gradient {
    atom_idx: usize,
    atom_sym: String,
    vector: [f64; 3],
}

#[derive(Default, Clone)]
struct BasisShell {
    atom_idx: usize,
    origin: [f64; 3],
    l: i32,
    exps: Vec<f64>,
    coeffs: Vec<f64>,
}

#[derive(Default, Clone)]
struct MoCoeff {
    coeffs: Vec<MoCoeffEntry>,
    spin: String,
    energy: f64,
    occ: f64,
    id: i32,
}

#[derive(Default, Clone)]
struct MoCoeffEntry {
    atom_idx: i32,
    sym: String,
    orb: String,
    coeff: f64,
}

// ---------------------------------------------------------------------------
// Main parser struct
// ---------------------------------------------------------------------------

struct OrcaRustParser<'a> {
    lines: Vec<&'a str>,
}

impl<'a> OrcaRustParser<'a> {
    fn new(content: &'a str) -> Self {
        OrcaRustParser {
            lines: content.lines().collect(),
        }
    }

    fn get(&self, i: usize) -> &str {
        self.lines.get(i).copied().unwrap_or("")
    }

    fn upper(&self, i: usize) -> String {
        self.get(i).to_uppercase()
    }

    fn len(&self) -> usize {
        self.lines.len()
    }

    // -----------------------------------------------------------------------
    // parse_all — orchestrates all parsing methods
    // -----------------------------------------------------------------------
    fn parse(&self) -> ParsedData {
        let mut data = ParsedData::default();
        self.parse_basic(&mut data);
        self.parse_gradients(&mut data);
        self.parse_trajectory(&mut data);
        self.parse_frequencies(&mut data);
        self.parse_thermal(&mut data);
        self.parse_orbital_energies(&mut data);
        self.parse_mo_coeffs(&mut data);
        self.parse_charges(&mut data);
        self.parse_dipole(&mut data);
        self.parse_tddft(&mut data);
        self.parse_nmr(&mut data);
        self.parse_basis_set(&mut data);
        self.parse_scf_trace(&mut data);
        self.parse_scan_results_table(&mut data);
        data
    }

    // -----------------------------------------------------------------------
    // parse_basic
    // -----------------------------------------------------------------------
    fn parse_basic(&self, data: &mut ParsedData) {
        for (i, &line) in self.lines.iter().enumerate() {
            if line.contains("Program Version") {
                if let Some(v) = line.split("Version").nth(1) {
                    if let Some(tok) = v.trim().split_whitespace().next() {
                        data.version = Some(tok.to_string());
                    }
                }
            }

            let stripped = line.trim();
            let uu = stripped.to_uppercase();

            if uu.contains("FINAL SINGLE POINT ENERGY") {
                if let Some(v) = stripped.split_whitespace().last().and_then(parse_f64) {
                    data.scf_energy = Some(v);
                }
            }
            if uu.contains("TOTAL CHARGE") {
                if let Some(v) = stripped.split_whitespace().last().and_then(|s| parse_i32(s)) {
                    data.charge = v;
                }
            }
            if uu.contains("MULTIPLICITY") {
                if let Some(v) = stripped.split_whitespace().last().and_then(|s| parse_i32(s)) {
                    data.mult = v;
                }
            }
            if uu.contains("SCF CONVERGED")
                || uu.contains("OPTIMIZATION CONVERGED")
                || uu.contains("HURRAY")
            {
                data.converged = true;
            }
            if uu.contains("RELAXED SURFACE SCAN") {
                data.is_scan = true;
            }
            if uu.contains("NUDGED ELASTIC BAND") || uu.contains(" NEB ") {
                data.is_neb = true;
            }
            if uu.contains("CURRENT TRAJECTORY WILL BE WRITTEN TO") {
                if let Some(cap) = RE_NEB_TRJ.captures(line) {
                    data.neb_trj_file = Some(cap[1].trim().to_string());
                } else {
                    data.neb_trj_file =
                        stripped.split_whitespace().last().map(|s| s.to_string());
                }
            }
            if uu.contains("CARTESIAN COORDINATES (ANGSTROEM)") {
                data.atoms.clear();
                data.coords.clear();
                let mut curr = i + 2;
                while curr < self.len() {
                    let lg = self.get(curr).trim();
                    if lg.is_empty() || lg.contains("---") {
                        break;
                    }
                    let parts: Vec<&str> = lg.split_whitespace().collect();
                    if parts.len() >= 4 {
                        if let (Some(x), Some(y), Some(z)) = (
                            parse_f64(parts[1]),
                            parse_f64(parts[2]),
                            parse_f64(parts[3]),
                        ) {
                            data.atoms.push(parts[0].to_string());
                            data.coords.push([x, y, z]);
                        }
                    }
                    curr += 1;
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // parse_gradients
    // -----------------------------------------------------------------------
    fn parse_gradients(&self, data: &mut ParsedData) {
        data.gradients.clear();
        data.all_gradients.clear();

        let mut gradient_starts: Vec<usize> = Vec::new();
        for (i, &line) in self.lines.iter().enumerate() {
            let uu = line.trim().to_uppercase();
            if uu.contains("CARTESIAN GRADIENT") && !uu.contains("NORM") {
                gradient_starts.push(i);
            }
        }

        for &start_idx in &gradient_starts {
            let mut block_grads: Vec<Gradient> = Vec::new();
            let mut curr = start_idx + 1;
            let mut found_data = false;

            while curr < self.len() && curr < start_idx + 15 {
                let line = self.get(curr).trim();
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 3 && parts[0].parse::<i32>().is_ok() {
                    found_data = true;
                    break;
                }
                curr += 1;
            }
            if !found_data {
                continue;
            }

            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.contains("-------") || line.contains("Difference to") || line.is_empty() {
                    break;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 5 {
                    let parsed = if parts.len() >= 6 && parts[2] == ":" {
                        parse_i32(parts[0]).and_then(|idx_raw| {
                            let vx = parse_f64(parts[3])?;
                            let vy = parse_f64(parts[4])?;
                            let vz = parse_f64(parts[5])?;
                            let idx = if idx_raw > 0 { idx_raw as usize - 1 } else { 0 };
                            Some((idx, parts[1].to_string(), vx, vy, vz))
                        })
                    } else {
                        parse_i32(parts[0]).and_then(|idx_raw| {
                            let vx = parse_f64(parts[2])?;
                            let vy = parse_f64(parts[3])?;
                            let vz = parse_f64(parts[4])?;
                            let idx = if idx_raw > 0 { idx_raw as usize - 1 } else { 0 };
                            Some((idx, parts[1].to_string(), vx, vy, vz))
                        })
                    };
                    if let Some((idx, sym, vx, vy, vz)) = parsed {
                        block_grads.push(Gradient {
                            atom_idx: idx,
                            atom_sym: sym,
                            vector: [vx, vy, vz],
                        });
                    }
                }
                curr += 1;
            }

            if !block_grads.is_empty() {
                data.all_gradients.push(GradientBlock {
                    line: start_idx,
                    grads: block_grads,
                });
            }
        }

        if let Some(last) = data.all_gradients.last() {
            data.gradients = last.grads.clone();
        }
    }

    // -----------------------------------------------------------------------
    // Helper: read coords from a given line index (look forward for CARTESIAN block)
    // -----------------------------------------------------------------------
    fn read_coords_from(&self, idx: usize) -> (Vec<String>, Vec<[f64; 3]>, bool) {
        let limit = 1000usize.min(self.len().saturating_sub(idx));
        for k in 0..limit {
            let line = self.get(idx + k).to_uppercase();
            if line.contains("CARTESIAN COORDINATES (ANGSTROEM)") {
                let mut c_idx = idx + k + 2;
                let mut atoms = Vec::new();
                let mut coords = Vec::new();
                while c_idx < self.len() {
                    let cl = self.get(c_idx).trim();
                    if cl.is_empty() || cl.contains("-------") {
                        break;
                    }
                    let parts: Vec<&str> = cl.split_whitespace().collect();
                    if parts.len() >= 4 {
                        if let (Some(x), Some(y), Some(z)) = (
                            parse_f64(parts[1]),
                            parse_f64(parts[2]),
                            parse_f64(parts[3]),
                        ) {
                            atoms.push(parts[0].to_string());
                            coords.push([x, y, z]);
                        }
                    }
                    c_idx += 1;
                }
                if !atoms.is_empty() {
                    return (atoms, coords, true);
                }
            }
        }
        (Vec::new(), Vec::new(), false)
    }

    // -----------------------------------------------------------------------
    // Helper: parse convergence info block starting at c_idx within next_marker
    // -----------------------------------------------------------------------
    fn parse_conv_info(
        &self,
        c_idx_start: usize,
        next_marker: usize,
    ) -> HashMap<String, ConvEntry> {
        let mut conv_info = HashMap::new();
        let mut c_idx = c_idx_start;
        while c_idx < next_marker && c_idx < c_idx_start + 30 {
            let cl = self.get(c_idx).trim();
            if cl.is_empty() || cl.contains("---") {
                c_idx += 1;
                continue;
            }
            let p: Vec<&str> = cl.split_whitespace().collect();
            if p.len() >= 4 {
                let s = p[p.len() - 1].to_uppercase();
                let t = p[p.len() - 2];
                let v = p[p.len() - 3];
                let n = p[..p.len() - 3].join(" ").to_lowercase();

                if s == "YES" || s == "NO" {
                    if !n.is_empty() && n != "item" {
                        conv_info.insert(
                            n,
                            ConvEntry {
                                value: v.to_string(),
                                tolerance: t.to_string(),
                                converged: s,
                            },
                        );
                    }
                } else if cl.to_lowercase().contains("max(") {
                    for cap in RE_MAX_STATS.captures_iter(cl) {
                        conv_info.insert(
                            cap[1].to_string(),
                            ConvEntry {
                                value: cap[2].to_string(),
                                tolerance: String::new(),
                                converged: "INFO".to_string(),
                            },
                        );
                    }
                }
            }
            c_idx += 1;
        }
        conv_info
    }

    // -----------------------------------------------------------------------
    // parse_trajectory
    // -----------------------------------------------------------------------
    fn parse_trajectory(&self, data: &mut ParsedData) {
        let mut current_scan_step: Option<i32> = None;

        for i in 0..self.len() {
            let uu_line = self.upper(i);

            // --- NEB PATH SUMMARY ---
            if uu_line.contains("PATH SUMMARY") && i > 0 && self.get(i - 1).contains("----") {
                let mut curr = i + 1;
                let mut header_found = false;
                while curr < self.len() && curr < i + 10 {
                    let l = self.get(curr);
                    if l.contains("Image") && l.contains("E(Eh)") {
                        header_found = true;
                        curr += 1;
                        break;
                    }
                    curr += 1;
                }
                if header_found {
                    while curr < self.len() {
                        let l_row = self.get(curr).trim();
                        if l_row.is_empty() {
                            break;
                        }
                        let parts: Vec<&str> = l_row.split_whitespace().collect();
                        if parts.len() >= 3 && parts[0].parse::<i32>().is_ok() {
                            if let (Some(img_idx), Some(dist), Some(en)) = (
                                parse_i32(parts[0]),
                                parse_f64(parts[1]),
                                parse_f64(parts[2]),
                            ) {
                                let step_num = data.scan_steps.len() as i32 + 1;
                                data.scan_steps.push(ScanStep {
                                    step_type: "neb_image".to_string(),
                                    scan_step_id: Some(img_idx),
                                    step: step_num,
                                    energy: en,
                                    scan_coord: Some(dist),
                                    dist: Some(dist),
                                    ..Default::default()
                                });
                            }
                        }
                        curr += 1;
                    }
                }
            }

            // --- RELAXED SURFACE SCAN STEP ---
            if uu_line.contains("RELAXED SURFACE SCAN STEP") {
                let step_idx = RE_SCAN_STEP
                    .captures(&uu_line)
                    .and_then(|c| parse_i32(&c[1]))
                    .unwrap_or(0);
                current_scan_step = Some(step_idx);

                let next_marker = (i + 1..self.len())
                    .find(|&m| {
                        self.lines[m]
                            .to_uppercase()
                            .contains("RELAXED SURFACE SCAN STEP")
                    })
                    .unwrap_or(self.len());

                let mut en = 0.0f64;
                let mut conv_info = HashMap::new();
                let mut coord_val: Option<f64> = None;

                for k in i..next_marker {
                    let uu = self.upper(k);
                    let raw = self.get(k);
                    if uu.contains("ACTUAL SCAN COORDINATE") {
                        coord_val = raw.split_whitespace().last().and_then(parse_f64);
                    }
                    if uu.contains("FINAL SINGLE POINT ENERGY") {
                        en = raw.split_whitespace().last().and_then(parse_f64).unwrap_or(en);
                    } else if uu.contains("TOTAL ENERGY") && raw.contains(':') && uu.contains("EH") {
                        if let Some(part) = raw.split(':').nth(1) {
                            en = part.split_whitespace().next().and_then(parse_f64).unwrap_or(en);
                        }
                    } else if uu.contains("CURRENT ENERGY") && raw.contains("....") {
                        if let Some(part) = raw.split("....").nth(1) {
                            en = part.split_whitespace().next().and_then(parse_f64).unwrap_or(en);
                        }
                    } else if uu.contains("GEOMETRY CONVERGENCE") || uu.contains("CONVERGENCE CRITERIA") {
                        conv_info = self.parse_conv_info(k + 1, next_marker);
                    }
                }

                // Attach gradients from within [i, next_marker)
                let step_grads: Vec<Gradient> = data
                    .all_gradients
                    .iter()
                    .filter(|g| g.line >= i && g.line < next_marker)
                    .last()
                    .map(|g| g.grads.clone())
                    .unwrap_or_default();

                let (atoms, coords, found) = self.read_coords_from(i);
                if found {
                    data.scan_steps.push(ScanStep {
                        step_type: "scan_step".to_string(),
                        scan_step_id: current_scan_step,
                        step: step_idx,
                        energy: en,
                        scan_coord: coord_val,
                        dist: None,
                        atoms,
                        coords,
                        convergence: conv_info,
                        gradients: step_grads,
                    });
                }
            }
            // --- OPTIMIZATION CYCLE ---
            else if uu_line.contains("OPTIMIZATION CYCLE") {
                let cycle_idx = RE_CYCLE
                    .captures(&uu_line)
                    .and_then(|c| parse_i32(&c[1]))
                    .unwrap_or(0);

                let next_marker = (i + 1..self.len())
                    .find(|&m| {
                        let u = self.lines[m].to_uppercase();
                        u.contains("OPTIMIZATION CYCLE")
                            || u.contains("OPTIMIZATION HAS CONVERGED")
                            || u.contains("OPTIMIZATION HAS RUN OUT OF CYCLES")
                            || u.contains("ORCA TERMINATED NORMALLY")
                    })
                    .unwrap_or(self.len());

                let mut en = 0.0f64;
                let mut conv_info = HashMap::new();

                for k in i..next_marker {
                    let uu = self.upper(k);
                    let raw = self.get(k);
                    if uu.contains("FINAL SINGLE POINT ENERGY") {
                        en = raw.split_whitespace().last().and_then(parse_f64).unwrap_or(en);
                    } else if uu.contains("TOTAL ENERGY") && raw.contains(':') && uu.contains("EH") {
                        if let Some(part) = raw.split(':').nth(1) {
                            en = part.split_whitespace().next().and_then(parse_f64).unwrap_or(en);
                        }
                    } else if uu.contains("CURRENT ENERGY") && raw.contains("....") {
                        if let Some(part) = raw.split("....").nth(1) {
                            en = part.split_whitespace().next().and_then(parse_f64).unwrap_or(en);
                        }
                    } else if uu.contains("GEOMETRY CONVERGENCE") || uu.contains("CONVERGENCE CRITERIA") {
                        conv_info = self.parse_conv_info(k + 1, next_marker);
                    }
                }

                let step_grads: Vec<Gradient> = data
                    .all_gradients
                    .iter()
                    .filter(|g| g.line >= i && g.line < next_marker)
                    .last()
                    .map(|g| g.grads.clone())
                    .unwrap_or_default();

                let (atoms, coords, found) = self.read_coords_from(i);
                if found {
                    data.scan_steps.push(ScanStep {
                        step_type: "opt_cycle".to_string(),
                        scan_step_id: current_scan_step,
                        step: cycle_idx,
                        energy: en,
                        scan_coord: None,
                        dist: None,
                        atoms,
                        coords,
                        convergence: conv_info,
                        gradients: step_grads,
                    });
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // parse_frequencies
    // -----------------------------------------------------------------------
    fn parse_frequencies(&self, data: &mut ParsedData) {
        data.frequencies.clear();

        // 1. Frequency values
        let mut freq_start: Option<usize> = None;
        for (i, &line) in self.lines.iter().enumerate() {
            if line.to_uppercase().contains("VIBRATIONAL FREQUENCIES") {
                freq_start = Some(i);
            }
        }
        if let Some(fs) = freq_start {
            let mut curr = fs + 1;
            while curr < self.len() && curr < fs + 10 {
                let l = self.get(curr);
                if l.contains(':') && (l.contains("cm**-1") || l.contains("cm-1")) {
                    break;
                }
                curr += 1;
            }
            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.contains("NORMAL MODES") {
                    break;
                }
                if line.contains("-------") && !data.frequencies.is_empty() {
                    break;
                }
                if line.is_empty() {
                    curr += 1;
                    continue;
                }
                if line.contains(':') && (line.contains("cm**-1") || line.contains("cm-1")) {
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 2 {
                        if let Some(v) = parse_f64(parts[1]) {
                            data.frequencies.push(FreqEntry {
                                freq: v,
                                ..Default::default()
                            });
                        }
                    }
                } else if !data.frequencies.is_empty() && !line.contains(':') {
                    break;
                }
                curr += 1;
            }
        }

        // 2. IR intensities
        let mut ir_start: Option<usize> = None;
        for (i, &line) in self.lines.iter().enumerate() {
            if line.to_uppercase().contains("IR SPECTRUM") {
                ir_start = Some(i);
            }
        }
        if let Some(is) = ir_start {
            let mut curr = is + 5;
            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.contains("The first frequency") {
                    break;
                }
                if line.contains("-----") && curr > is + 10 {
                    break;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() > 3 && parts[0].contains(':') {
                    let idx_str = parts[0].replace(':', "");
                    if let (Some(idx), Some(inten)) =
                        (parse_i32(&idx_str), parse_f64(parts[3]))
                    {
                        let idx = idx as usize;
                        if idx < data.frequencies.len() {
                            data.frequencies[idx].ir = inten;
                        }
                    }
                }
                curr += 1;
            }
        }

        // 3. Raman activities
        let mut raman_start: Option<usize> = None;
        for (i, &line) in self.lines.iter().enumerate() {
            if line.to_uppercase().contains("RAMAN SPECTRUM") {
                raman_start = Some(i);
            }
        }
        if let Some(rs) = raman_start {
            let mut curr = rs + 5;
            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.contains("The first frequency") {
                    break;
                }
                if line.contains("-----") && curr > rs + 10 {
                    break;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() > 2 && parts[0].contains(':') {
                    let idx_str = parts[0].replace(':', "");
                    if let (Some(idx), Some(act)) =
                        (parse_i32(&idx_str), parse_f64(parts[2]))
                    {
                        let idx = idx as usize;
                        if idx < data.frequencies.len() {
                            data.frequencies[idx].raman = act;
                        }
                    }
                }
                curr += 1;
            }
        }

        // 4. Normal modes
        let mut modes_start: Option<usize> = None;
        for (i, &line) in self.lines.iter().enumerate() {
            if line.to_uppercase().contains("NORMAL MODES") {
                modes_start = Some(i);
            }
        }
        if let (Some(ms), false) = (modes_start, data.atoms.is_empty()) {
            let n_atoms = data.atoms.len();
            let n_coords = n_atoms * 3;
            let mut curr = ms + 7;
            let mut mode_buffer: HashMap<usize, Vec<f64>> = HashMap::new();

            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.is_empty() {
                    curr += 1;
                    continue;
                }
                if line.to_uppercase().contains("IR SPECTRUM") || line.contains("--------") {
                    break;
                }
                // Try parsing as header (integers)
                let parts: Vec<&str> = line.split_whitespace().collect();
                let headers: Option<Vec<usize>> = parts
                    .iter()
                    .map(|s| s.parse::<usize>().ok())
                    .collect::<Option<Vec<_>>>();

                if let Some(headers) = headers {
                    let start_data = curr + 1;
                    for r in 0..n_coords {
                        let dl = self.get(start_data + r);
                        let dparts: Vec<&str> = dl.split_whitespace().collect();
                        if dparts.len() < 2 {
                            continue;
                        }
                        let values: Vec<f64> = dparts[1..]
                            .iter()
                            .filter_map(|s| parse_f64(s))
                            .collect();
                        for (c, &m_idx) in headers.iter().enumerate() {
                            if c < values.len() {
                                mode_buffer.entry(m_idx).or_default().push(values[c]);
                            }
                        }
                    }
                    curr = start_data + n_coords;
                } else {
                    curr += 1;
                }
            }

            for (m_idx, vec_flat) in mode_buffer {
                if m_idx < data.frequencies.len() {
                    let vecs: Vec<[f64; 3]> = vec_flat
                        .chunks(3)
                        .filter(|c| c.len() == 3)
                        .map(|c| [c[0], c[1], c[2]])
                        .collect();
                    data.frequencies[m_idx].vector = vecs;
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // parse_thermal
    // -----------------------------------------------------------------------
    fn parse_thermal(&self, data: &mut ParsedData) {
        let thermal_keys: Vec<(&str, &str)> = vec![
            ("ELECTRONIC ENERGY", "electronic_energy"),
            ("ZERO POINT ENERGY", "zpe"),
            ("THERMAL VIBRATIONAL CORRECTION", "corr_vib"),
            ("THERMAL ROTATIONAL CORRECTION", "corr_rot"),
            ("THERMAL TRANSLATIONAL CORRECTION", "corr_trans"),
            ("TOTAL THERMAL ENERGY", "thermal_energy"),
            ("TOTAL THERMAL CORRECTION", "corr_thermal_total"),
            ("NON-THERMAL (ZPE) CORRECTION", "corr_zpe"),
            ("TOTAL CORRECTION", "corr_total"),
            ("TOTAL ENTHALPY", "enthalpy"),
            ("THERMAL ENTHALPY CORRECTION", "thermal_enthalpy_corr"),
            ("ELECTRONIC ENTROPY", "s_el"),
            ("VIBRATIONAL ENTROPY", "s_vib"),
            ("ROTATIONAL ENTROPY", "s_rot"),
            ("TRANSLATIONAL ENTROPY", "s_trans"),
            ("FINAL ENTROPY TERM", "entropy"),
            ("FINAL GIBBS FREE ENERGY", "gibbs"),
            ("G-E(EL)", "gibbs_corr"),
        ];

        let mut start_line: Option<usize> = None;
        for (i, &line) in self.lines.iter().enumerate() {
            if line.to_uppercase().contains("THERMOCHEMISTRY AT") {
                start_line = Some(i);
            }
        }
        let start_line = match start_line {
            Some(s) => s,
            None => return,
        };

        let mut curr = start_line;
        while curr < self.len() {
            let line = self.get(curr);
            let uu = line.to_uppercase();

            if uu.contains("TIMINGS FOR INDIVIDUAL MODULES") {
                break;
            }
            if uu.contains("TEMPERATURE") && uu.contains("K") && uu.contains("...") {
                if let Some(cap) = RE_TEMP.captures(line) {
                    if let Some(v) = parse_f64(&cap[1]) {
                        data.thermal.insert("temperature".to_string(), v);
                    }
                }
            }
            for (key_upper, val_key) in &thermal_keys {
                if uu.contains(key_upper) {
                    let val = if line.contains("Eh") {
                        let pre_eh = line.split("Eh").next().unwrap_or("");
                        all_floats(pre_eh).last().copied()
                    } else {
                        all_floats(line).last().copied()
                    };
                    if let Some(v) = val {
                        data.thermal.insert(val_key.to_string(), v);
                    }
                }
            }
            curr += 1;
        }

        // Post-processing
        let enthalpy = data.thermal.get("enthalpy").copied();
        let el_en = data.thermal.get("electronic_energy").copied();
        if let (Some(h), Some(e)) = (enthalpy, el_en) {
            data.thermal.insert("enthalpy_corr".to_string(), h - e);
        }
        let gibbs = data.thermal.get("gibbs").copied();
        if let (Some(g), Some(e)) = (gibbs, el_en) {
            data.thermal.insert("gibbs_corr".to_string(), g - e);
        }

        let imag_count = data
            .frequencies
            .iter()
            .filter(|f| f.freq < 0.0)
            .count() as f64;
        data.thermal.insert("imaginary_freq_count".to_string(), imag_count);
    }

    // -----------------------------------------------------------------------
    // parse_orbital_energies
    // -----------------------------------------------------------------------
    fn parse_orbital_energies(&self, data: &mut ParsedData) {
        data.orbital_energies.clear();

        let mut start_indices: Vec<(usize, String)> = Vec::new();
        for (i, &line) in self.lines.iter().enumerate() {
            let uu = line.to_uppercase();
            let next_has_dashes = i + 1 < self.len() && self.get(i + 1).contains("---");
            if uu.contains("ORBITAL ENERGIES") && next_has_dashes {
                start_indices.push((i, "restricted".to_string()));
            } else if uu.contains("SPIN UP ORBITALS") && next_has_dashes {
                start_indices.push((i, "alpha".to_string()));
            } else if uu.contains("SPIN DOWN ORBITALS") && next_has_dashes {
                start_indices.push((i, "beta".to_string()));
            }
        }

        // Keep last occurrence of each spin type
        let mut final_map: HashMap<String, usize> = HashMap::new();
        for (idx, spin) in start_indices {
            final_map.insert(spin, idx);
        }
        let mut filtered: Vec<(usize, String)> = final_map.into_iter().map(|(s, i)| (i, s)).collect();
        filtered.sort_by_key(|(i, _)| *i);

        for (start_idx, spin) in filtered {
            let mut curr = start_idx + 2;
            // Find column header
            while curr < self.len() && curr < start_idx + 10 {
                let line = self.get(curr);
                if line.contains("NO") && (line.contains("OCC") || line.contains("E(Eh)") || line.contains("E(eV)")) {
                    curr += 1;
                    break;
                }
                curr += 1;
            }
            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.is_empty() || line.contains("---") || line.contains("****") || line.contains("MULLIKEN") {
                    break;
                }
                if line.starts_with('*') {
                    curr += 1;
                    continue;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 4 {
                    let first = parts[0].trim_start_matches('-');
                    if let (Some(orb_idx), Some(occ), Some(e_eh), Some(e_ev)) = (
                        parse_i32(first),
                        parse_f64(parts[1]),
                        parse_f64(parts[2]),
                        parse_f64(parts[3]),
                    ) {
                        let orb_type = if occ > 0.1 { "occupied" } else { "virtual" };
                        data.orbital_energies.push(OrbitalEnergy {
                            index: orb_idx,
                            occupation: occ,
                            energy_eh: e_eh,
                            energy_ev: e_ev,
                            spin: spin.clone(),
                            orb_type: orb_type.to_string(),
                        });
                    }
                }
                curr += 1;
            }
        }
    }

    // -----------------------------------------------------------------------
    // parse_mo_coeffs
    // -----------------------------------------------------------------------
    fn parse_mo_coeffs(&self, data: &mut ParsedData) {
        let mut start_indices: Vec<(usize, String)> = Vec::new();
        for (i, &line) in self.lines.iter().enumerate() {
            let uu = line.to_uppercase();
            let next_has_dashes = i + 1 < self.len() && self.get(i + 1).contains("---");
            if uu.contains("MOLECULAR ORBITALS") && next_has_dashes {
                start_indices.push((i, "restricted".to_string()));
            } else if uu.contains("SPIN UP ORBITALS") && next_has_dashes {
                start_indices.push((i, "alpha".to_string()));
            } else if uu.contains("SPIN DOWN ORBITALS") && next_has_dashes {
                start_indices.push((i, "beta".to_string()));
            }
        }

        // Keep last occurrence per spin type
        let mut final_map: HashMap<String, usize> = HashMap::new();
        for (idx, spin) in start_indices {
            final_map.insert(spin, idx);
        }
        let mut filtered: Vec<(usize, String)> = final_map.into_iter().map(|(s, i)| (i, s)).collect();
        filtered.sort_by_key(|(i, _)| *i);

        for (start_idx, spin) in filtered {
            let mut curr = start_idx + 2;
            let header_line = self.get(start_idx).to_uppercase();
            let mut current_spin = spin.clone();
            if spin == "restricted" && header_line.contains("(UHF)") {
                current_spin = "alpha".to_string();
            }

            let mut current_mos: Vec<i32> = Vec::new();
            let mut last_first_mo_idx: i32 = -1;

            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.is_empty() {
                    curr += 1;
                    continue;
                }
                if line.contains("TIMINGS") {
                    break;
                }
                if line.contains("--------") {
                    curr += 1;
                    continue;
                }
                if line.contains("ORBITALS") && curr + 1 < self.len() && self.get(curr + 1).contains("--------") {
                    break;
                }

                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.is_empty() {
                    curr += 1;
                    continue;
                }

                // Check if header line (all integers)
                let all_int = parts.iter().all(|p| p.parse::<i32>().is_ok());
                if all_int && !parts.is_empty() {
                    current_mos = parts.iter().filter_map(|p| parse_i32(p)).collect();

                    // Detect spin switch
                    if let Some(&first) = current_mos.first() {
                        if first <= last_first_mo_idx && last_first_mo_idx != -1 {
                            if current_spin == "alpha" || current_spin == "restricted" {
                                current_spin = "beta".to_string();
                            }
                        }
                        last_first_mo_idx = first;
                    }

                    // Initialize MO storage
                    for &idx in &current_mos {
                        let key = format!("{}_{}", idx, current_spin);
                        data.mo_coeffs.entry(key).or_insert_with(|| MoCoeff {
                            spin: current_spin.clone(),
                            id: idx,
                            ..Default::default()
                        });
                    }

                    // Parse energy / occ lines that follow
                    if curr + 2 < self.len() {
                        let next1: Vec<&str> = self.get(curr + 1).split_whitespace().collect();
                        let next2: Vec<&str> = self.get(curr + 2).split_whitespace().collect();
                        if next1.len() == current_mos.len() && next2.len() == current_mos.len() {
                            let energies: Option<Vec<f64>> =
                                next1.iter().map(|s| parse_f64(s)).collect();
                            let occs: Option<Vec<f64>> =
                                next2.iter().map(|s| parse_f64(s)).collect();
                            if let (Some(energies), Some(occs)) = (energies, occs) {
                                for (k, &mo_idx) in current_mos.iter().enumerate() {
                                    let key = format!("{}_{}", mo_idx, current_spin);
                                    if let Some(entry) = data.mo_coeffs.get_mut(&key) {
                                        entry.energy = energies[k];
                                        entry.occ = occs[k];
                                    }
                                }
                                curr += 2;
                            }
                        }
                    }
                    curr += 1;
                    continue;
                }

                // Coefficient line
                if parts.len() >= 2 {
                    let (atom_idx, sym, orb, val_strs): (i32, String, String, Vec<&str>) =
                        if let Some(cap) = RE_MERGED_ATOM.captures(parts[0]) {
                            if parts.len() >= 3 {
                                let ai = parse_i32(&cap[1]).unwrap_or(-1);
                                (ai, cap[2].to_string(), parts[1].to_string(), parts[2..].to_vec())
                            } else {
                                curr += 1;
                                continue;
                            }
                        } else if parts.len() >= 3 && parts[0].parse::<i32>().is_ok() {
                            let ai = parse_i32(parts[0]).unwrap_or(-1);
                            (ai, parts[1].to_string(), parts[2].to_string(), parts[3..].to_vec())
                        } else {
                            curr += 1;
                            continue;
                        };

                    if val_strs.len() == current_mos.len() {
                        for (k, v_str) in val_strs.iter().enumerate() {
                            if let Some(val) = parse_f64(v_str) {
                                let mo_idx = current_mos[k];
                                let key = format!("{}_{}", mo_idx, current_spin);
                                if let Some(entry) = data.mo_coeffs.get_mut(&key) {
                                    entry.coeffs.push(MoCoeffEntry {
                                        atom_idx,
                                        sym: sym.clone(),
                                        orb: orb.clone(),
                                        coeff: val,
                                    });
                                }
                            }
                        }
                    }
                }
                curr += 1;
            }
        }
    }

    // -----------------------------------------------------------------------
    // parse_charges
    // -----------------------------------------------------------------------
    fn parse_charges(&self, data: &mut ParsedData) {
        let mut mulliken_start: i64 = -1;
        let mut loewdin_start: i64 = -1;
        let mut hirshfeld_start: i64 = -1;
        let mut mayer_start: i64 = -1;
        let mut nbo_start: i64 = -1;
        let mut chelpg_start: i64 = -1;
        let mut mk_start: i64 = -1;
        let mut mbis_start: i64 = -1;
        let mut resp_start: i64 = -1;
        let mut fmo_start: i64 = -1;

        for (i, &line) in self.lines.iter().enumerate() {
            let uu = line.to_uppercase();
            if uu.contains("MULLIKEN ATOMIC CHARGES") { mulliken_start = i as i64; }
            else if uu.contains("LOEWDIN ATOMIC CHARGES") { loewdin_start = i as i64; }
            else if uu.contains("HIRSHFELD ANALYSIS") { hirshfeld_start = i as i64; }
            else if uu.contains("MAYER POPULATION ANALYSIS") { mayer_start = i as i64; }
            else if uu.contains("NATURAL POPULATIONS") { nbo_start = i as i64; }
            else if uu.contains("CHELPG ATOMIC CHARGES") { chelpg_start = i as i64; }
            else if uu.contains("MERZ-KOLLMAN ATOMIC CHARGES") || uu.contains("MK ATOMIC CHARGES") { mk_start = i as i64; }
            else if uu.contains("MBIS ANALYSIS") { mbis_start = i as i64; }
            else if uu.contains("RESP ATOMIC CHARGES") { resp_start = i as i64; }
            else if uu.contains("FRONTIER MOLECULAR ORBITAL POPULATION ANALYSIS") { fmo_start = i as i64; }
        }

        // Standard block parser
        let parse_standard = |start_idx: i64, header_lines: usize, hirshfeld: bool, mbis: bool| -> Vec<ChargeEntry> {
            let mut res = Vec::new();
            if start_idx < 0 { return res; }
            let mut curr = start_idx as usize + header_lines;
            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.is_empty() || line.contains("---") || line.contains("Sum of") {
                    if !res.is_empty() { break; }
                    curr += 1;
                    continue;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 3 {
                    let (idx_str, sym) = if let Some(cap) = RE_MERGED_ATOM.captures(parts[0]) {
                        (cap[1].to_string(), cap[2].to_string())
                    } else {
                        let idx = parts[0].trim_end_matches(':').to_string();
                        let sym = parts[1].trim_end_matches(':').to_string();
                        (idx, sym)
                    };
                    if let Some(atom_idx) = parse_i32(&idx_str) {
                        let val = if hirshfeld || mbis {
                            parse_f64(parts[2])
                        } else if parts.len() >= 4 && parts[2] == ":" {
                            parse_f64(parts[3])
                        } else {
                            parse_f64(parts[2])
                        };
                        if let Some(charge) = val {
                            let mut entry = ChargeEntry {
                                atom_idx,
                                atom_sym: sym,
                                charge,
                                ..Default::default()
                            };
                            if hirshfeld && parts.len() >= 4 {
                                entry.spin = parse_f64(parts[3]);
                            } else if mbis {
                                if parts.len() >= 4 { entry.population = parse_f64(parts[3]); }
                                if parts.len() >= 5 { entry.spin = parse_f64(parts[4]); }
                            }
                            res.push(entry);
                        }
                    }
                }
                curr += 1;
            }
            res
        };

        let mulliken = parse_standard(mulliken_start, 2, false, false);
        if !mulliken.is_empty() { data.charges.insert("Mulliken".to_string(), mulliken); }

        let loewdin = parse_standard(loewdin_start, 2, false, false);
        if !loewdin.is_empty() { data.charges.insert("Loewdin".to_string(), loewdin); }

        let hirshfeld = parse_standard(hirshfeld_start, 2, true, false);
        if !hirshfeld.is_empty() { data.charges.insert("Hirshfeld".to_string(), hirshfeld); }

        let chelpg = parse_standard(chelpg_start, 2, false, false);
        if !chelpg.is_empty() { data.charges.insert("CHELPG".to_string(), chelpg); }

        let mk = parse_standard(mk_start, 2, false, false);
        if !mk.is_empty() { data.charges.insert("MK".to_string(), mk); }

        let mbis = parse_standard(mbis_start, 3, false, true);
        if !mbis.is_empty() { data.charges.insert("MBIS".to_string(), mbis); }

        let resp = parse_standard(resp_start, 2, false, false);
        if !resp.is_empty() { data.charges.insert("RESP".to_string(), resp); }

        // Mayer
        if mayer_start >= 0 {
            let mut mayer_res: Vec<ChargeEntry> = Vec::new();
            let mut curr = mayer_start as usize + 1;
            while curr < self.len() && curr < mayer_start as usize + 15 {
                let l = self.get(curr);
                if l.contains("ATOM") && l.contains("QA") {
                    curr += 1;
                    break;
                }
                curr += 1;
            }
            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.is_empty() || line.contains("---") || line.contains("Mayer bond") {
                    break;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 4 {
                    if let (Some(idx), Some(qa)) = (parse_i32(parts[0]), parse_f64(parts.get(4).copied().unwrap_or(""))) {
                        let mut entry = ChargeEntry {
                            atom_idx: idx,
                            atom_sym: parts[1].to_string(),
                            charge: qa,
                            ..Default::default()
                        };
                        if parts.len() >= 8 {
                            entry.valency = parse_f64(parts[5]);
                            entry.bonded_valency = parse_f64(parts[6]);
                            entry.free_valency = parse_f64(parts[7]);
                        }
                        mayer_res.push(entry);
                    }
                }
                curr += 1;
            }
            if !mayer_res.is_empty() {
                if !data.charges.contains_key("Mulliken") {
                    data.charges.insert("Mulliken".to_string(), mayer_res.clone());
                }
                data.charges.insert("Mayer".to_string(), mayer_res);
            }
        }

        // NBO
        if nbo_start >= 0 {
            let mut nbo_charges: Vec<ChargeEntry> = Vec::new();
            let nbo_s = nbo_start as usize;

            // Try summary table first
            let mut summary_start: i64 = -1;
            for i in nbo_s..nbo_s.saturating_add(2000).min(self.len()) {
                if self.get(i).contains("Summary of Natural Population Analysis") {
                    summary_start = i as i64;
                    break;
                }
            }

            if summary_start >= 0 {
                let mut curr = summary_start as usize + 1;
                while curr < self.len() {
                    let l = self.get(curr);
                    if l.contains("Atom No") && l.contains("Charge") && l.contains("Core") {
                        curr += 1;
                        if curr < self.len() && self.get(curr).contains("----") { curr += 1; }
                        break;
                    }
                    curr += 1;
                }
                while curr < self.len() {
                    let line = self.get(curr).trim();
                    if line.contains("====") || line.contains("Total") { break; }
                    if line.is_empty() { curr += 1; continue; }
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 7 {
                        if let (Some(idx), Some(chg), Some(core), Some(val), Some(ryd), Some(tot)) = (
                            parse_i32(parts[1]),
                            parse_f64(parts[2]),
                            parse_f64(parts[3]),
                            parse_f64(parts[4]),
                            parse_f64(parts[5]),
                            parse_f64(parts[6]),
                        ) {
                            nbo_charges.push(ChargeEntry {
                                atom_idx: idx - 1,
                                atom_sym: parts[0].to_string(),
                                charge: chg,
                                core: Some(core),
                                valence: Some(val),
                                rydberg: Some(ryd),
                                total: Some(tot),
                                ..Default::default()
                            });
                        }
                    }
                    curr += 1;
                }
            }

            if nbo_charges.is_empty() {
                let mut curr = nbo_s + 1;
                while curr < self.len() {
                    let line = self.get(curr).trim();
                    if line.contains("---") { curr += 1; continue; }
                    if line.contains("================") || line.contains("Natural Electron Configuration") {
                        if !nbo_charges.is_empty() { break; }
                        curr += 1; continue;
                    }
                    if line.is_empty() {
                        if !nbo_charges.is_empty() { break; }
                        curr += 1; continue;
                    }
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 3 {
                        if let (Some(idx), Some(chg)) = (parse_i32(parts[1]), parse_f64(parts[2])) {
                            nbo_charges.push(ChargeEntry {
                                atom_idx: idx - 1,
                                atom_sym: parts[0].to_string(),
                                charge: chg,
                                ..Default::default()
                            });
                        }
                    }
                    curr += 1;
                }
            }

            if !nbo_charges.is_empty() {
                data.charges.insert("NBO".to_string(), nbo_charges);
            }
        }

        // FMO
        if fmo_start >= 0 {
            let mut fmo_data: Vec<ChargeEntry> = Vec::new();
            let mut curr = fmo_start as usize + 1;
            let mut table_start = false;

            while curr < self.len() && curr < fmo_start as usize + 40 {
                if self.get(curr).contains("--------") {
                    let prev = if curr > 0 { self.get(curr - 1) } else { "" };
                    let prev2 = if curr > 1 { self.get(curr - 2) } else { "" };
                    if prev.contains("Atom") || prev2.contains("HOMO") {
                        table_start = true;
                        curr += 1;
                        break;
                    }
                }
                curr += 1;
            }

            if table_start {
                while curr < self.len() {
                    let line = self.get(curr).trim();
                    if line.is_empty() || line.contains("--------") {
                        if !fmo_data.is_empty() { break; }
                        curr += 1; continue;
                    }
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 5 {
                        let atom_lbl = parts[0];
                        if atom_lbl.contains('-') {
                            let sub: Vec<&str> = atom_lbl.splitn(2, '-').collect();
                            if sub.len() == 2 {
                                if let (Some(idx), Some(hm), Some(hl), Some(lm), Some(ll)) = (
                                    parse_i32(sub[0]),
                                    parse_f64(parts[1]),
                                    parse_f64(parts[2]),
                                    parse_f64(parts[3]),
                                    parse_f64(parts[4]),
                                ) {
                                    fmo_data.push(ChargeEntry {
                                        atom_idx: idx,
                                        atom_sym: sub[1].to_string(),
                                        charge: hm,
                                        homo_mulliken: Some(hm),
                                        homo_loewdin: Some(hl),
                                        lumo_mulliken: Some(lm),
                                        lumo_loewdin: Some(ll),
                                        ..Default::default()
                                    });
                                }
                            }
                        }
                    }
                    curr += 1;
                }
                if !fmo_data.is_empty() {
                    data.charges.insert("FMO".to_string(), fmo_data);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // parse_dipole
    // -----------------------------------------------------------------------
    fn parse_dipole(&self, data: &mut ParsedData) {
        let candidates: Vec<usize> = self
            .lines
            .iter()
            .enumerate()
            .filter(|(_, l)| {
                let uu = l.to_uppercase();
                uu.contains("TOTAL DIPOLE MOMENT") && l.contains(':')
            })
            .map(|(i, _)| i)
            .collect();

        let idx = match candidates.last() {
            Some(&i) => i,
            None => return,
        };
        let line = self.get(idx);
        if let Some(rest) = line.split(':').nth(1) {
            let parts: Vec<&str> = rest.trim().split_whitespace().collect();
            if parts.len() >= 3 {
                if let (Some(x), Some(y), Some(z)) =
                    (parse_f64(parts[0]), parse_f64(parts[1]), parse_f64(parts[2]))
                {
                    let mag = if idx + 1 < self.len() {
                        let line2 = self.get(idx + 1);
                        if line2.contains("Magnitude") && line2.contains(':') {
                            line2.split(':').nth(1).and_then(|s| parse_f64(s.trim())).unwrap_or_else(|| (x*x+y*y+z*z).sqrt())
                        } else {
                            (x * x + y * y + z * z).sqrt()
                        }
                    } else {
                        (x * x + y * y + z * z).sqrt()
                    };
                    data.dipole = Some(Dipole { x, y, z, magnitude: mag });
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // parse_tddft
    // -----------------------------------------------------------------------
    fn parse_tddft(&self, data: &mut ParsedData) {
        let mut states_dict: HashMap<i32, TddftState> = HashMap::new();

        // Pass 1: Detailed excited state blocks
        let mut current_state_id: i32 = -1;
        for &line in &self.lines {
            let line_upper = line.to_uppercase();

            if let Some(cap) = RE_STATE_HEADER.captures(&line_upper) {
                if let Some(id) = parse_i32(&cap[1]) {
                    current_state_id = id;
                    let entry = get_or_insert_state(&mut states_dict, id);
                    if let Some(cap) = RE_ENERGY_EV.captures(line) {
                        entry.energy_ev = parse_f64(&cap[1]).unwrap_or(0.0);
                    }
                    if let Some(cap) = RE_ENERGY_CM.captures(line) {
                        entry.energy_cm = parse_f64(&cap[1]).unwrap_or(0.0);
                    }
                    if let Some(cap) = RE_ENERGY_NM.captures(line) {
                        entry.energy_nm = parse_f64(&cap[1]).unwrap_or(0.0);
                    }
                } else {
                    current_state_id = -1;
                }
            } else if current_state_id >= 0 {
                if !line.contains("-------") && !line_upper.contains("SPECTRUM") {
                    if line.contains("->") && line.contains(':') {
                        let parts: Vec<&str> = line.trim().splitn(2, ':').collect();
                        if parts.len() >= 2 {
                            let trans_desc = parts[0].trim();
                            let coeff = parts[1].trim();
                            let t_str = format!("{} (coeff: {})", trans_desc, coeff);
                            let entry = get_or_insert_state(&mut states_dict, current_state_id);
                            if !entry.transitions.contains(&t_str) {
                                entry.transitions.push(t_str);
                            }
                        }
                    }
                }
            }
        }

        // Pass 2: Summary tables
        let parse_summary_table = |states: &mut HashMap<i32, TddftState>, lines: &[&str], start_idx: usize, data_key: &str| {
            let mut curr = start_idx + 1;
            let mut header_found = false;

            while curr < lines.len() && curr < start_idx + 30 {
                let line = lines[curr].trim();
                if line.is_empty() || line.contains("--------") {
                    curr += 1;
                    continue;
                }
                let u_line = line.to_uppercase();
                if (u_line.contains("TRANSITION") || u_line.contains("STATE")) && (u_line.contains("ENERGY") || u_line.contains("WAVELENGTH")) {
                    header_found = true;
                    break;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 3 && parts[0].parse::<i32>().is_ok() {
                    header_found = true;
                    break;
                }
                curr += 1;
            }
            if !header_found { return; }

            let mut data_parsing_started = false;
            while curr < lines.len() {
                let line = lines[curr].trim();
                if line.is_empty() { curr += 1; continue; }
                if line.to_uppercase().contains("TOTAL") { break; }
                if line.contains("--------") {
                    if !data_parsing_started {
                        data_parsing_started = true;
                        curr += 1;
                        continue;
                    } else {
                        break;
                    }
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if !parts.is_empty() && parts[0].parse::<i32>().is_ok() {
                    data_parsing_started = true;
                }

                if parts.contains(&"->") {
                    if let Some(arrow_idx) = parts.iter().position(|&p| p == "->") {
                        if parts.len() > arrow_idx + 5 {
                            let target = parts[arrow_idx + 1];
                            if let Some(cap) = RE_TARGET_STATE.captures(target) {
                                if let Some(s_id) = parse_i32(&cap[1]) {
                                    let entry = states.entry(s_id).or_insert_with(|| TddftState { state: s_id, ..Default::default() });
                                    let strength = parse_f64(parts[arrow_idx + 5]).unwrap_or(0.0);
                                    match data_key {
                                        "osc_len" => entry.osc_len = strength,
                                        "osc" => entry.osc = strength,
                                        "osc_vel" => entry.osc_vel = strength,
                                        "rot_len" => entry.rot_len = strength,
                                        "rotatory_strength" => entry.rotatory_strength = strength,
                                        "rot_vel" => entry.rot_vel = strength,
                                        _ => {}
                                    }
                                    if let Some(ev) = parse_f64(parts[arrow_idx + 2]) {
                                        if entry.energy_ev == 0.0 { entry.energy_ev = ev; }
                                    }
                                    if let Some(cm) = parse_f64(parts[arrow_idx + 3]) {
                                        if entry.energy_cm == 0.0 { entry.energy_cm = cm; }
                                    }
                                    if let Some(nm) = parse_f64(parts[arrow_idx + 4]) {
                                        if entry.energy_nm == 0.0 { entry.energy_nm = nm; }
                                    }
                                }
                            }
                        }
                    }
                } else if parts.len() >= 4 && parts[0].parse::<i32>().is_ok() {
                    if let Some(s_id) = parse_i32(parts[0]) {
                        let entry = states.entry(s_id).or_insert_with(|| TddftState { state: s_id, ..Default::default() });
                        let strength = parse_f64(parts[3]).unwrap_or(0.0);
                        match data_key {
                            "osc_len" => entry.osc_len = strength,
                            "osc" => entry.osc = strength,
                            "osc_vel" => entry.osc_vel = strength,
                            "rot_len" => entry.rot_len = strength,
                            "rotatory_strength" => entry.rotatory_strength = strength,
                            "rot_vel" => entry.rot_vel = strength,
                            _ => {}
                        }
                        if entry.energy_ev == 0.0 { entry.energy_ev = parse_f64(parts[1]).unwrap_or(0.0); }
                        if entry.energy_nm == 0.0 { entry.energy_nm = parse_f64(parts[2]).unwrap_or(0.0); }
                    }
                }
                curr += 1;
            }
        };

        for (i, &line) in self.lines.iter().enumerate() {
            let line_upper = line.to_uppercase();
            if line_upper.contains("ABSORPTION SPECTRUM") {
                if line_upper.contains("ELECTRIC DIPOLE") {
                    parse_summary_table(&mut states_dict, &self.lines, i, "osc_len");
                    parse_summary_table(&mut states_dict, &self.lines, i, "osc");
                } else if line_upper.contains("VELOCITY DIPOLE") {
                    parse_summary_table(&mut states_dict, &self.lines, i, "osc_vel");
                }
            }
            if line_upper.contains("CD SPECTRUM") {
                if line_upper.contains("ELECTRIC DIPOLE") {
                    parse_summary_table(&mut states_dict, &self.lines, i, "rot_len");
                    parse_summary_table(&mut states_dict, &self.lines, i, "rotatory_strength");
                } else if line_upper.contains("VELOCITY DIPOLE") {
                    parse_summary_table(&mut states_dict, &self.lines, i, "rot_vel");
                }
            }
        }

        // Finalization
        for item in states_dict.values_mut() {
            if item.energy_ev == 0.0 {
                if item.energy_nm > 0.1 {
                    item.energy_ev = 1239.84193 / item.energy_nm;
                } else if item.energy_cm > 0.1 {
                    item.energy_ev = item.energy_cm / 8065.54425;
                }
            }
        }

        let mut valid: Vec<TddftState> = states_dict.into_values().filter(|s| s.energy_ev > 0.0).collect();
        valid.sort_by(|a, b| a.energy_ev.partial_cmp(&b.energy_ev).unwrap());
        data.tddft = valid;
    }

    // -----------------------------------------------------------------------
    // parse_nmr
    // -----------------------------------------------------------------------
    fn parse_nmr(&self, data: &mut ParsedData) {
        // Shielding summary
        let mut summary_start: i64 = -1;
        for (i, &line) in self.lines.iter().enumerate() {
            if line.to_uppercase().contains("CHEMICAL SHIELDING SUMMARY (PPM)") {
                summary_start = i as i64;
                break;
            }
        }

        if summary_start >= 0 {
            let mut curr = summary_start as usize + 1;
            while curr < self.len() && curr < summary_start as usize + 10 {
                let l_up = self.get(curr).to_uppercase();
                if (l_up.contains('N') && l_up.contains("SHIELDING"))
                    || (l_up.contains("NUCLEUS") && (l_up.contains("ISOTROPIC") || l_up.contains("ELEMENT")))
                {
                    curr += 1;
                    break;
                }
                curr += 1;
            }
            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.is_empty() { break; }
                if line.contains("---") {
                    if !data.nmr_shielding.is_empty() { break; }
                    curr += 1; continue;
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 3 {
                    if let (Some(idx), Some(val)) = (parse_i32(parts[0]), parse_f64(parts[2])) {
                        data.nmr_shielding.push(NmrShielding {
                            atom_idx: idx,
                            atom_sym: parts[1].to_string(),
                            shielding: val,
                        });
                    }
                }
                curr += 1;
            }
        }

        // Couplings
        let mut coupling_start: i64 = -1;
        for (i, &line) in self.lines.iter().enumerate() {
            if line.to_uppercase().contains("SUMMARY OF ISOTROPIC COUPLING CONSTANTS") {
                coupling_start = i as i64;
                break;
            }
        }

        if coupling_start >= 0 {
            let mut curr = coupling_start as usize + 1;
            while curr < self.len() {
                if self.get(curr).contains("----------------") {
                    curr += 1;
                    break;
                }
                curr += 1;
            }

            let mut current_col_indices: Vec<i32> = Vec::new();
            while curr < self.len() {
                let line = self.get(curr).trim();
                if line.is_empty() { curr += 1; continue; }
                if line.contains("Maximum memory used") || line.contains("Timings") || line.contains("ORCA TERMINATED") { break; }

                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.is_empty() { curr += 1; continue; }

                // Header detection: 3rd token is integer
                let is_header = if parts.len() >= 3 {
                    parts[2].parse::<i32>().is_ok()
                } else if parts.len() == 2 {
                    true // "0 C" only
                } else {
                    false
                };

                if is_header {
                    current_col_indices.clear();
                    let mut p_idx = 0;
                    while p_idx + 1 < parts.len() {
                        if parts[p_idx].parse::<i32>().is_ok() {
                            current_col_indices.push(parts[p_idx].parse().unwrap());
                            p_idx += 2;
                        } else {
                            p_idx += 1;
                        }
                    }
                } else if parts.len() >= 2 && parts[0].parse::<i32>().is_ok() {
                    if let Some(row_atom_idx) = parse_i32(parts[0]) {
                        let values = &parts[2..];
                        for (c_i, val_str) in values.iter().enumerate() {
                            if c_i < current_col_indices.len() {
                                let col_atom_idx = current_col_indices[c_i];
                                if let Some(val) = parse_f64(val_str) {
                                    if row_atom_idx < col_atom_idx {
                                        data.nmr_couplings.push(NmrCoupling {
                                            atom_idx1: row_atom_idx,
                                            atom_idx2: col_atom_idx,
                                            coupling: val,
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
                curr += 1;
            }
        }
    }

    // -----------------------------------------------------------------------
    // parse_basis_set
    // -----------------------------------------------------------------------
    fn parse_basis_set(&self, data: &mut ParsedData) {
        let start_idx = match self.lines.iter().position(|l| l.contains("BASIS SET IN INPUT FORMAT")) {
            Some(i) => i,
            None => return,
        };

        let mut curr = start_idx + 2;
        let mut basis_defs: HashMap<String, Vec<BasisShell>> = HashMap::new();
        let mut current_sym: Option<String> = None;
        let mut current_shells: Vec<(i32, Vec<f64>, Vec<f64>)> = Vec::new(); // (l, exps, coeffs)

        while curr < self.len() {
            let line = self.get(curr).trim();

            if line.contains("--------") && curr > start_idx + 10 { break; }
            if line.contains("AUXILIARY BASIS") { break; }

            if line.starts_with("NewGTO") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Some(sym) = &current_sym {
                        let shells = current_shells.drain(..).map(|(l, exps, coeffs)| BasisShell {
                            l, exps, coeffs, ..Default::default()
                        }).collect();
                        basis_defs.insert(sym.clone(), shells);
                    }
                    current_sym = Some(parts[1].to_string());
                    current_shells.clear();
                }
                curr += 1;
                continue;
            }

            if line.starts_with("end") {
                if let Some(sym) = &current_sym {
                    let shells = current_shells.drain(..).map(|(l, exps, coeffs)| BasisShell {
                        l, exps, coeffs, ..Default::default()
                    }).collect();
                    basis_defs.insert(sym.clone(), shells);
                }
                current_sym = None;
                curr += 1;
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            let l_map: HashMap<&str, i32> = [("S", 0), ("P", 1), ("D", 2), ("F", 3), ("G", 4)].iter().cloned().collect();
            if parts.len() >= 2 && l_map.contains_key(parts[0].to_uppercase().as_str()) {
                let sh_type = parts[0].to_uppercase();
                if let Some(&l_val) = l_map.get(sh_type.as_str()) {
                    if let Some(n_prim) = parse_i32(parts[1]) {
                        if n_prim <= 50 {
                            curr += 1;
                            let mut exps = Vec::new();
                            let mut coeffs = Vec::new();
                            for _ in 0..n_prim {
                                if curr >= self.len() { break; }
                                let pl = self.get(curr).trim();
                                let pp: Vec<&str> = pl.split_whitespace().collect();
                                if pp.len() >= 3 {
                                    if let (Some(e), Some(c)) = (parse_f64(pp[1]), parse_f64(pp[2])) {
                                        exps.push(e);
                                        coeffs.push(c);
                                    }
                                }
                                curr += 1;
                            }
                            if !exps.is_empty() {
                                current_shells.push((l_val, exps, coeffs));
                            }
                            continue;
                        }
                    }
                }
            }

            curr += 1;
            if curr > start_idx + 5000 { break; }
        }

        // Expand to actual atoms
        if data.atoms.is_empty() { return; }
        for (idx, (sym, coord)) in data.atoms.iter().zip(data.coords.iter()).enumerate() {
            if let Some(shells) = basis_defs.get(sym) {
                for s in shells {
                    data.basis_set_shells.push(BasisShell {
                        atom_idx: idx,
                        origin: *coord,
                        l: s.l,
                        exps: s.exps.clone(),
                        coeffs: s.coeffs.clone(),
                    });
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // parse_scf_trace
    // -----------------------------------------------------------------------
    fn parse_scf_trace(&self, data: &mut ParsedData) {
        let mut current_step_label = "Initial".to_string();
        let mut i = 0;

        while i < self.len() {
            let line = self.get(i);
            let uu = line.to_uppercase();

            if uu.contains("OPTIMIZATION CYCLE") {
                if let Some(cap) = RE_CYCLE.captures(&uu) {
                    current_step_label = format!("Cycle {}", &cap[1]);
                } else {
                    current_step_label = "Opt Cycle".to_string();
                }
            } else if uu.contains("SCAN STEP") {
                current_step_label = line.split_whitespace().last()
                    .map(|s| format!("Scan Step {}", s))
                    .unwrap_or_else(|| "Scan Step".to_string());
            } else if uu.contains("ORCA PROPERTIES") || uu.contains("ORCA PROPERTY") {
                current_step_label = "Property/Final".to_string();
            } else if uu.contains("OPTIMIZATION HAS CONVERGED") {
                current_step_label = "Post-Opt/Final".to_string();
            }

            if uu.contains("SCF ITERATIONS") || uu.contains("ORCA LEAN-SCF") || uu.contains("INCREMENTAL FOCK MATRIX") {
                let mut header_idx: i64 = -1;
                for k in 1..15usize {
                    if i + k >= self.len() { break; }
                    let uu_k = self.get(i + k).to_uppercase();
                    if uu_k.contains("ITER") && uu_k.contains("ENERGY") {
                        header_idx = (i + k) as i64;
                        break;
                    }
                }

                if header_idx >= 0 {
                    let mut trace: Vec<ScfIter> = Vec::new();
                    let mut idx = header_idx as usize + 1;
                    if idx < self.len() && self.get(idx).contains("---") {
                        idx += 1;
                    }

                    while idx < self.len() {
                        let l_scf = self.get(idx).trim();
                        if l_scf.is_empty() || l_scf.contains("---") || l_scf.contains("SUCCESS") || l_scf.contains("Energy Check") {
                            if !trace.is_empty() { break; }
                            idx += 1;
                            continue;
                        }
                        let parts: Vec<&str> = l_scf.split_whitespace().collect();
                        if parts.len() >= 2 {
                            if let (Some(it_no), Some(it_en)) = (parse_i32(parts[0]), parse_f64(parts[1])) {
                                trace.push(ScfIter { iter: it_no, energy: it_en });
                            }
                        }
                        idx += 1;
                    }

                    if !trace.is_empty() {
                        let same_count = data.scf_traces.iter()
                            .filter(|t| t.step.starts_with(&current_step_label))
                            .count();
                        let label = if same_count > 0 {
                            format!("{} ({})", current_step_label, same_count + 1)
                        } else {
                            current_step_label.clone()
                        };
                        data.scf_traces.push(ScfTrace { step: label, iterations: trace });
                        i = idx;
                        continue;
                    }
                }
            }
            i += 1;
        }
    }

    // -----------------------------------------------------------------------
    // parse_scan_results_table
    // -----------------------------------------------------------------------
    fn parse_scan_results_table(&self, data: &mut ParsedData) {
        let data_start = self.lines.iter().enumerate().rev()
            .find(|(_, l)| l.contains("Actual Energy"))
            .map(|(i, _)| i + 1);
        let data_start = match data_start {
            Some(s) => s,
            None => return,
        };

        let mut table_vals: Vec<(f64, f64)> = Vec::new();
        for i in data_start..self.len() {
            let line = self.get(i).trim();
            if line.is_empty() {
                if !table_vals.is_empty() { break; }
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                if let (Some(coord), Some(en)) = (parse_f64(parts[0]), parse_f64(parts[parts.len() - 1])) {
                    table_vals.push((coord, en));
                } else if !table_vals.is_empty() {
                    break;
                }
            }
        }

        if table_vals.is_empty() { return; }

        if data.scan_steps.is_empty() {
            for (idx, (coord, en)) in table_vals.iter().enumerate() {
                data.scan_steps.push(ScanStep {
                    step_type: "scan_step_summary".to_string(),
                    scan_step_id: Some(idx as i32),
                    step: idx as i32,
                    energy: *en,
                    scan_coord: Some(*coord),
                    ..Default::default()
                });
            }
        } else {
            let sids: Vec<i32> = data.scan_steps.iter()
                .filter_map(|s| s.scan_step_id)
                .collect();
            if !sids.is_empty() {
                let offset = *sids.iter().min().unwrap();
                for s in data.scan_steps.iter_mut() {
                    if let Some(sid) = s.scan_step_id {
                        let idx = (sid - offset) as usize;
                        if idx < table_vals.len() {
                            s.scan_coord = Some(table_vals[idx].0);
                        }
                    }
                }
            } else {
                for (idx, s) in data.scan_steps.iter_mut().enumerate() {
                    if idx < table_vals.len() {
                        s.scan_coord = Some(table_vals[idx].0);
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// TDDFT helper — standalone fn avoids closure lifetime issues
// ---------------------------------------------------------------------------

fn get_or_insert_state(states: &mut HashMap<i32, TddftState>, idx: i32) -> &mut TddftState {
    states.entry(idx).or_insert_with(move || TddftState {
        state: idx,
        ..Default::default()
    })
}

// ---------------------------------------------------------------------------
// XYZ multi-frame parser
// ---------------------------------------------------------------------------

fn do_parse_xyz_content<'py>(py: Python<'py>, content: &str) -> PyResult<Bound<'py, PyList>> {
    let lines: Vec<&str> = content.lines().collect();
    let result = PyList::empty_bound(py);
    let mut i = 0;
    let n_lines = lines.len();

    let re_ts = Regex::new(r"\bTS\b").unwrap();
    let re_ci = Regex::new(r"\bCI\b").unwrap();

    while i < n_lines {
        let line = lines[i].trim();
        if line.is_empty() { i += 1; continue; }

        let natoms = match line.parse::<usize>() {
            Ok(n) => n,
            Err(_) => { i += 1; continue; }
        };
        i += 1;
        if i >= n_lines { break; }

        let comment = lines[i].trim();
        let upper_comment = comment.to_uppercase();

        // Filter TS / CI (but not CI-NEB method)
        let is_excluded = re_ts.is_match(&upper_comment)
            || (re_ci.is_match(&upper_comment) && !upper_comment.contains("CI-NEB"));

        if is_excluded {
            i += 1 + natoms;
            continue;
        }

        // Extract energy
        let energy = RE_ENERGY_LABEL.captures(comment)
            .and_then(|c| parse_f64(&c[1]))
            .or_else(|| {
                RE_FLOATS_IN_LINE.find_iter(comment).last().and_then(|m| parse_f64(m.as_str()))
            })
            .unwrap_or(0.0);

        // Extract distance/coord
        let dist_val = RE_DIST_LABEL.captures(comment).and_then(|c| parse_f64(&c[1]));

        i += 1;
        let mut atoms: Vec<String> = Vec::new();
        let mut coords: Vec<[f64; 3]> = Vec::new();

        for _ in 0..natoms {
            if i >= n_lines { break; }
            let parts: Vec<&str> = lines[i].split_whitespace().collect();
            if parts.len() >= 4 {
                if let (Some(x), Some(y), Some(z)) = (parse_f64(parts[1]), parse_f64(parts[2]), parse_f64(parts[3])) {
                    atoms.push(parts[0].to_string());
                    coords.push([x, y, z]);
                }
            }
            i += 1;
        }

        let step = PyDict::new_bound(py);
        step.set_item("type", "neb_step")?;
        step.set_item("energy", energy)?;
        step.set_item("dist", dist_val.into_py(py))?;
        step.set_item("scan_coord", dist_val.into_py(py))?;

        let atoms_list = PyList::new_bound(py, &atoms);
        step.set_item("atoms", atoms_list)?;

        let coords_list = PyList::empty_bound(py);
        for c in &coords {
            let row = PyList::new_bound(py, c);
            coords_list.append(row)?;
        }
        step.set_item("coords", coords_list)?;
        result.append(step)?;
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Convert ParsedData to Python dict
// ---------------------------------------------------------------------------

fn to_python<'py>(py: Python<'py>, data: &ParsedData) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);

    // Scalars
    d.set_item("scf_energy", data.scf_energy.into_py(py))?;
    d.set_item("converged", data.converged)?;
    d.set_item("charge", data.charge)?;
    d.set_item("mult", data.mult)?;
    d.set_item("version", data.version.clone().into_py(py))?;
    d.set_item("is_scan", data.is_scan)?;
    d.set_item("is_neb", data.is_neb)?;
    d.set_item("neb_trj_file", data.neb_trj_file.clone().into_py(py))?;

    // atoms / coords
    let atoms_list = PyList::new_bound(py, &data.atoms);
    d.set_item("atoms", atoms_list)?;

    let coords_list = PyList::empty_bound(py);
    for c in &data.coords {
        let row = PyList::new_bound(py, c);
        coords_list.append(row)?;
    }
    d.set_item("coords", coords_list)?;

    // scf_traces
    let scf_traces_list = PyList::empty_bound(py);
    for trace in &data.scf_traces {
        let t = PyDict::new_bound(py);
        t.set_item("step", trace.step.clone())?;
        let iters = PyList::empty_bound(py);
        for it in &trace.iterations {
            let it_d = PyDict::new_bound(py);
            it_d.set_item("iter", it.iter)?;
            it_d.set_item("energy", it.energy)?;
            iters.append(it_d)?;
        }
        t.set_item("iterations", iters)?;
        scf_traces_list.append(t)?;
    }
    d.set_item("scf_traces", scf_traces_list)?;

    // frequencies
    let freq_list = PyList::empty_bound(py);
    for f in &data.frequencies {
        let fd = PyDict::new_bound(py);
        fd.set_item("freq", f.freq)?;
        fd.set_item("ir", f.ir)?;
        fd.set_item("raman", f.raman)?;
        let vl = PyList::empty_bound(py);
        for v in &f.vector {
            let row = PyList::new_bound(py, v);
            vl.append(row)?;
        }
        fd.set_item("vector", vl)?;
        freq_list.append(fd)?;
    }
    d.set_item("frequencies", freq_list)?;

    // orbital_energies + mos (backward compat)
    let oe_list = PyList::empty_bound(py);
    for oe in &data.orbital_energies {
        let od = PyDict::new_bound(py);
        od.set_item("index", oe.index)?;
        od.set_item("id", oe.index)?;
        od.set_item("occupation", oe.occupation)?;
        od.set_item("occ", oe.occupation)?;
        od.set_item("energy_eh", oe.energy_eh)?;
        od.set_item("energy_ev", oe.energy_ev)?;
        od.set_item("energy", oe.energy_eh)?;
        od.set_item("spin", oe.spin.clone())?;
        od.set_item("type", oe.orb_type.clone())?;
        oe_list.append(od)?;
    }
    d.set_item("orbital_energies", &oe_list)?;
    d.set_item("mos", oe_list)?;

    // charges
    let charges_dict = PyDict::new_bound(py);
    for (k, entries) in &data.charges {
        let el = PyList::empty_bound(py);
        for e in entries {
            let ed = PyDict::new_bound(py);
            ed.set_item("atom_idx", e.atom_idx)?;
            ed.set_item("atom_sym", e.atom_sym.clone())?;
            ed.set_item("charge", e.charge)?;
            if let Some(v) = e.spin { ed.set_item("spin", v)?; }
            if let Some(v) = e.population { ed.set_item("population", v)?; }
            if let Some(v) = e.valency { ed.set_item("valency", v)?; }
            if let Some(v) = e.bonded_valency { ed.set_item("bonded_valency", v)?; }
            if let Some(v) = e.free_valency { ed.set_item("free_valency", v)?; }
            if let Some(v) = e.core { ed.set_item("core", v)?; }
            if let Some(v) = e.valence { ed.set_item("valence", v)?; }
            if let Some(v) = e.rydberg { ed.set_item("rydberg", v)?; }
            if let Some(v) = e.total { ed.set_item("total", v)?; }
            if let Some(v) = e.homo_mulliken { ed.set_item("homo_mulliken", v)?; }
            if let Some(v) = e.homo_loewdin { ed.set_item("homo_loewdin", v)?; }
            if let Some(v) = e.lumo_mulliken { ed.set_item("lumo_mulliken", v)?; }
            if let Some(v) = e.lumo_loewdin { ed.set_item("lumo_loewdin", v)?; }
            el.append(ed)?;
        }
        charges_dict.set_item(k, el)?;
    }
    d.set_item("charges", charges_dict)?;

    // dipole
    match &data.dipole {
        None => {
            d.set_item("dipole", py.None())?;
            d.set_item("dipoles", py.None())?;
        }
        Some(dip) => {
            let dd = PyDict::new_bound(py);
            let vec_list = PyList::new_bound(py, &[dip.x, dip.y, dip.z]);
            dd.set_item("vector", vec_list)?;
            dd.set_item("magnitude", dip.magnitude)?;
            d.set_item("dipole", &dd)?;
            d.set_item("dipoles", dd)?;
        }
    }

    // tddft
    let tddft_list = PyList::empty_bound(py);
    for s in &data.tddft {
        let sd = PyDict::new_bound(py);
        sd.set_item("state", s.state)?;
        sd.set_item("energy_ev", s.energy_ev)?;
        sd.set_item("energy_nm", s.energy_nm)?;
        sd.set_item("energy_cm", s.energy_cm)?;
        sd.set_item("osc", s.osc)?;
        sd.set_item("osc_len", s.osc_len)?;
        sd.set_item("osc_vel", s.osc_vel)?;
        sd.set_item("rotatory_strength", s.rotatory_strength)?;
        sd.set_item("rot_len", s.rot_len)?;
        sd.set_item("rot_vel", s.rot_vel)?;
        let tl = PyList::new_bound(py, &s.transitions);
        sd.set_item("transitions", tl)?;
        tddft_list.append(sd)?;
    }
    d.set_item("tddft", tddft_list)?;
    // backward compat alias
    d.set_item("excitation_energies", PyList::empty_bound(py))?;

    // nmr_shielding
    let nmr_s_list = PyList::empty_bound(py);
    for ns in &data.nmr_shielding {
        let nd = PyDict::new_bound(py);
        nd.set_item("atom_idx", ns.atom_idx)?;
        nd.set_item("atom_sym", ns.atom_sym.clone())?;
        nd.set_item("shielding", ns.shielding)?;
        nmr_s_list.append(nd)?;
    }
    d.set_item("nmr_shielding", nmr_s_list)?;

    // nmr_couplings
    let nmr_c_list = PyList::empty_bound(py);
    for nc in &data.nmr_couplings {
        let nd = PyDict::new_bound(py);
        nd.set_item("atom_idx1", nc.atom_idx1)?;
        nd.set_item("atom_idx2", nc.atom_idx2)?;
        nd.set_item("coupling", nc.coupling)?;
        nmr_c_list.append(nd)?;
    }
    d.set_item("nmr_couplings", nmr_c_list)?;

    // thermal
    let thermal_dict = PyDict::new_bound(py);
    for (k, v) in &data.thermal {
        thermal_dict.set_item(k, v)?;
    }
    d.set_item("thermal", thermal_dict)?;

    // scan_steps
    let steps_list = PyList::empty_bound(py);
    for s in &data.scan_steps {
        let sd = PyDict::new_bound(py);
        sd.set_item("type", s.step_type.clone())?;
        sd.set_item("scan_step_id", s.scan_step_id.into_py(py))?;
        sd.set_item("step", s.step)?;
        sd.set_item("energy", s.energy)?;
        sd.set_item("scan_coord", s.scan_coord.into_py(py))?;
        sd.set_item("dist", s.dist.into_py(py))?;
        let al = PyList::new_bound(py, &s.atoms);
        sd.set_item("atoms", al)?;
        let cl = PyList::empty_bound(py);
        for c in &s.coords {
            let row = PyList::new_bound(py, c);
            cl.append(row)?;
        }
        sd.set_item("coords", cl)?;

        let conv_d = PyDict::new_bound(py);
        for (ck, cv) in &s.convergence {
            let ce = PyDict::new_bound(py);
            ce.set_item("value", cv.value.clone())?;
            ce.set_item("tolerance", cv.tolerance.clone())?;
            ce.set_item("converged", cv.converged.clone())?;
            conv_d.set_item(ck, ce)?;
        }
        sd.set_item("convergence", conv_d)?;

        let grad_l = PyList::empty_bound(py);
        for g in &s.gradients {
            let gd = PyDict::new_bound(py);
            gd.set_item("atom_idx", g.atom_idx)?;
            gd.set_item("atom_sym", g.atom_sym.clone())?;
            let vl = PyList::new_bound(py, &g.vector);
            gd.set_item("vector", vl)?;
            grad_l.append(gd)?;
        }
        sd.set_item("gradients", grad_l)?;
        steps_list.append(sd)?;
    }
    d.set_item("scan_steps", &steps_list)?;

    // all_steps (alias, same as scan_steps)
    let steps_alias = PyList::empty_bound(py);
    for item in steps_list.iter() { steps_alias.append(item)?; }
    d.set_item("all_steps", steps_alias)?;

    // gradients (last gradient block)
    let grad_list = PyList::empty_bound(py);
    for g in &data.gradients {
        let gd = PyDict::new_bound(py);
        gd.set_item("atom_idx", g.atom_idx)?;
        gd.set_item("atom_sym", g.atom_sym.clone())?;
        let vl = PyList::new_bound(py, &g.vector);
        gd.set_item("vector", vl)?;
        grad_list.append(gd)?;
    }
    d.set_item("gradients", grad_list)?;

    // all_gradients
    let all_grad_list = PyList::empty_bound(py);
    for gb in &data.all_gradients {
        let gbd = PyDict::new_bound(py);
        gbd.set_item("line", gb.line)?;
        let gl = PyList::empty_bound(py);
        for g in &gb.grads {
            let gd = PyDict::new_bound(py);
            gd.set_item("atom_idx", g.atom_idx)?;
            gd.set_item("atom_sym", g.atom_sym.clone())?;
            let vl = PyList::new_bound(py, &g.vector);
            gd.set_item("vector", vl)?;
            gl.append(gd)?;
        }
        gbd.set_item("grads", gl)?;
        all_grad_list.append(gbd)?;
    }
    d.set_item("all_gradients", all_grad_list)?;

    // basis_set_shells
    let bss_list = PyList::empty_bound(py);
    for shell in &data.basis_set_shells {
        let sd = PyDict::new_bound(py);
        sd.set_item("atom_idx", shell.atom_idx)?;
        let ol = PyList::new_bound(py, &shell.origin);
        sd.set_item("origin", ol)?;
        sd.set_item("l", shell.l)?;
        let el = PyList::new_bound(py, &shell.exps);
        sd.set_item("exps", el)?;
        let cl = PyList::new_bound(py, &shell.coeffs);
        sd.set_item("coeffs", cl)?;
        bss_list.append(sd)?;
    }
    d.set_item("basis_set_shells", bss_list)?;

    // mo_coeffs
    let mo_dict = PyDict::new_bound(py);
    for (key, mo) in &data.mo_coeffs {
        let md = PyDict::new_bound(py);
        md.set_item("id", mo.id)?;
        md.set_item("spin", mo.spin.clone())?;
        md.set_item("energy", mo.energy)?;
        md.set_item("occ", mo.occ)?;
        let cl = PyList::empty_bound(py);
        for c in &mo.coeffs {
            let cd = PyDict::new_bound(py);
            cd.set_item("atom_idx", c.atom_idx)?;
            cd.set_item("sym", c.sym.clone())?;
            cd.set_item("orb", c.orb.clone())?;
            cd.set_item("coeff", c.coeff)?;
            cl.append(cd)?;
        }
        md.set_item("coeffs", cl)?;
        mo_dict.set_item(key, md)?;
    }
    d.set_item("mo_coeffs", mo_dict)?;

    Ok(d)
}

// ---------------------------------------------------------------------------
// PyO3 exported functions
// ---------------------------------------------------------------------------

/// Parse a complete ORCA output file content string.
/// Returns a Python dict with all parsed data (same schema as OrcaParser.data).
#[pyfunction]
fn parse_all(py: Python, content: &str) -> PyResult<PyObject> {
    let parser = OrcaRustParser::new(content);
    let data = parser.parse();
    let d = to_python(py, &data)?;
    Ok(d.into())
}

/// Parse multi-frame XYZ content (e.g. NEB trajectory files).
/// Returns a Python list of step dicts.
#[pyfunction]
fn parse_xyz_content(py: Python, content: &str) -> PyResult<PyObject> {
    let result = do_parse_xyz_content(py, content)?;
    Ok(result.into())
}

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------

#[pymodule]
fn orca_parser_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_all, m)?)?;
    m.add_function(wrap_pyfunction!(parse_xyz_content, m)?)?;
    m.add("__version__", "1.0.0")?;
    Ok(())
}

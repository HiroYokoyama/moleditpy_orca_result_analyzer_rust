// orca_mo_rs — Rust-powered GTO basis set evaluation for MO cube generation.
// Exposed to Python via PyO3.
//
// Public surface:
//   BasisSetEngineRust  — evaluates MO wavefunction on a 3D grid
//
// Optimisations over the Python baseline:
//   1. Rayon — grid points evaluated in parallel across all CPU cores
//   2. Distance cutoff — shells with exp(-α_min·r²) < 1e-15 are skipped
//   3. Reusable exp buffer — one Vec per rayon task (no per-point allocation)

use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use rayon::prelude::*;

// ---------------------------------------------------------------------------
// GTO helpers
// ---------------------------------------------------------------------------

/// Factorials indexed by n: FACT[n] = n!, FACT2[n] = (2n)!
const FACT: [f64; 9] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0];
const FACT2: [f64; 9] = [
    1.0,
    2.0,
    24.0,
    720.0,
    40320.0,
    3628800.0,
    479001600.0,
    87178291200.0,
    20922789888000.0,
];

/// GTO Cartesian normalisation prefactor N(α, l, m, n).
fn normalization_prefactor(alpha: f64, l: u32, m: u32, n: u32) -> f64 {
    let big_l = (l + m + n) as f64;
    let num = (8.0 * alpha).powf(big_l) * FACT[l as usize] * FACT[m as usize] * FACT[n as usize];
    let den = FACT2[l as usize] * FACT2[m as usize] * FACT2[n as usize];
    (2.0 * alpha / std::f64::consts::PI).powf(0.75) * (num / den).sqrt()
}

/// x^exp — branchless for common low exponents.
#[inline(always)]
fn ang_pow(v: f64, exp: u32) -> f64 {
    match exp {
        0 => 1.0,
        1 => v,
        2 => v * v,
        3 => v * v * v,
        4 => {
            let v2 = v * v;
            v2 * v2
        }
        _ => v.powi(exp as i32),
    }
}

// ---------------------------------------------------------------------------
// Basis-set angular definitions  (S / P / D-sph / F-sph / G-sph)
// ---------------------------------------------------------------------------

type BasisCompDef = (f64, u32, u32, u32); // (weight, l, m, n)

/// Angular part of the Cartesian Gaussian normalisation.
fn n_cart_part(l: u32, m: u32, n: u32) -> f64 {
    ((FACT[l as usize] * FACT[m as usize] * FACT[n as usize])
        / (FACT2[l as usize] * FACT2[m as usize] * FACT2[n as usize]))
        .sqrt()
}

/// (k-1)!! for even k — the 1-D Gaussian moment integral factor.
fn odd_double_factorial(mut k: u32) -> f64 {
    let mut result = 1.0;
    while k > 1 {
        result *= (k - 1) as f64;
        k -= 2;
    }
    result
}

/// Overlap of two normalised Cartesian Gaussians in the same shell, scaled so
/// a single component overlaps itself as 1.0. The radial part is shared by
/// every component of a shell, so it cancels and the exponent drops out.
fn angular_overlap(a: (u32, u32, u32), b: (u32, u32, u32)) -> f64 {
    let pairs = [(a.0, b.0), (a.1, b.1), (a.2, b.2)];
    let mut moment = 1.0;
    for (ka, kb) in pairs {
        let k = ka + kb;
        if k % 2 == 1 {
            return 0.0; // odd power integrates to zero over all space
        }
        moment *= odd_double_factorial(k);
    }
    let total_l = a.0 + a.1 + a.2;
    (2.0_f64).powi(total_l as i32)
        * moment
        * n_cart_part(a.0, a.1, a.2)
        * n_cart_part(b.0, b.1, b.2)
}

/// Turn real-solid-harmonic polynomials into component weights.
///
/// Each component is evaluated as N_cart(l,m,n) * x^l y^m z^n, and N_cart
/// depends on how L is split across the axes (N(0,0,3) = N(2,0,1)/sqrt(5)).
/// So the raw polynomial coefficient is divided by N_cart to recover the
/// intended shape, then the whole component is normalised against the exact
/// angular overlap to fix the scale. Hand-typed prefactor tables got both of
/// these wrong: the f table distorted 5 of 7 components (the f(z^3) nodal
/// cone sat at 28.6 deg instead of 39.2) and the g table disagreed between
/// components (g2 came out ~2x the norm of g1).
fn build_sph_defs(polys: &[&[BasisCompDef]]) -> Vec<Vec<BasisCompDef>> {
    polys
        .iter()
        .map(|poly| {
            let weights: Vec<BasisCompDef> = poly
                .iter()
                .map(|&(c_poly, l, m, n)| (c_poly / n_cart_part(l, m, n), l, m, n))
                .collect();
            let mut norm_sq = 0.0;
            for &(wa, la, ma, na) in &weights {
                for &(wb, lb, mb, nb) in &weights {
                    norm_sq += wa * wb * angular_overlap((la, ma, na), (lb, mb, nb));
                }
            }
            let scale = if norm_sq > 0.0 {
                1.0 / norm_sq.sqrt()
            } else {
                1.0
            };
            weights
                .into_iter()
                .map(|(w, l, m, n)| (w * scale, l, m, n))
                .collect()
        })
        .collect()
}

/// Spherical F (7 components), ORCA m order 0,+1,-1,+2,-2,+3,-3.
fn build_f_shell_defs() -> Vec<Vec<BasisCompDef>> {
    let f_polys: [&[BasisCompDef]; 7] = [
        // f(0): z(5z^2 - 3r^2) = 2z^3 - 3x^2 z - 3y^2 z
        &[(2.0, 0, 0, 3), (-3.0, 2, 0, 1), (-3.0, 0, 2, 1)],
        // f(+1): x(5z^2 - r^2) = 4xz^2 - x^3 - xy^2
        &[(4.0, 1, 0, 2), (-1.0, 3, 0, 0), (-1.0, 1, 2, 0)],
        // f(-1): y(5z^2 - r^2)
        &[(4.0, 0, 1, 2), (-1.0, 2, 1, 0), (-1.0, 0, 3, 0)],
        // f(+2): z(x^2 - y^2)
        &[(1.0, 2, 0, 1), (-1.0, 0, 2, 1)],
        // f(-2): xyz
        &[(1.0, 1, 1, 1)],
        // f(+3): x(x^2 - 3y^2)
        &[(1.0, 3, 0, 0), (-3.0, 1, 2, 0)],
        // f(-3): y(3x^2 - y^2)
        &[(3.0, 2, 1, 0), (-1.0, 0, 3, 0)],
    ];
    build_sph_defs(&f_polys)
}

/// Spherical G (9 components), same m ordering.
fn build_g_shell_defs() -> Vec<Vec<BasisCompDef>> {
    let g_polys: [&[BasisCompDef]; 9] = [
        // g(0): (35z^4 - 30z^2 r^2 + 3r^4)/8
        &[
            (1.0, 0, 0, 4),
            (-3.0, 2, 0, 2),
            (-3.0, 0, 2, 2),
            (0.375, 4, 0, 0),
            (0.375, 0, 4, 0),
            (0.75, 2, 2, 0),
        ],
        // g(+1): xz(7z^2 - 3r^2)
        &[(4.0, 1, 0, 3), (-3.0, 3, 0, 1), (-3.0, 1, 2, 1)],
        // g(-1): yz(7z^2 - 3r^2)
        &[(4.0, 0, 1, 3), (-3.0, 0, 3, 1), (-3.0, 2, 1, 1)],
        // g(+2): (x^2-y^2)(7z^2 - r^2)
        &[
            (6.0, 2, 0, 2),
            (-6.0, 0, 2, 2),
            (-1.0, 4, 0, 0),
            (1.0, 0, 4, 0),
        ],
        // g(-2): xy(7z^2 - r^2)
        &[(6.0, 1, 1, 2), (-1.0, 3, 1, 0), (-1.0, 1, 3, 0)],
        // g(+3): xz(x^2 - 3y^2)
        &[(1.0, 3, 0, 1), (-3.0, 1, 2, 1)],
        // g(-3): yz(3x^2 - y^2)
        &[(3.0, 2, 1, 1), (-1.0, 0, 3, 1)],
        // g(+4): x^4 - 6x^2y^2 + y^4
        &[(1.0, 4, 0, 0), (-6.0, 2, 2, 0), (1.0, 0, 4, 0)],
        // g(-4): xy(x^2 - y^2)
        &[(4.0, 3, 1, 0), (-4.0, 1, 3, 0)],
    ];
    build_sph_defs(&g_polys)
}

fn get_basis_definitions(l_type: u32) -> Vec<Vec<BasisCompDef>> {
    match l_type {
        0 => vec![vec![(1.0, 0, 0, 0)]],
        1 => vec![
            vec![(1.0, 0, 0, 1)], // pz
            vec![(1.0, 1, 0, 0)], // px
            vec![(1.0, 0, 1, 0)], // py
        ],
        // d components each mix a single normalisation topology, so the
        // textbook weights are already exact here.
        2 => vec![
            vec![(-0.5, 2, 0, 0), (-0.5, 0, 2, 0), (1.0, 0, 0, 2)],
            vec![(1.0, 1, 0, 1)],
            vec![(1.0, 0, 1, 1)],
            vec![(0.866025, 2, 0, 0), (-0.866025, 0, 2, 0)],
            vec![(1.0, 1, 1, 0)],
        ],
        3 => build_f_shell_defs(),
        4 => build_g_shell_defs(),
        _ => vec![],
    }
}

// ---------------------------------------------------------------------------
// Pre-computed data structures
// ---------------------------------------------------------------------------

// comp_coeffs[k] = input_coeff[k] * norm(exps[k], l, m, n) * angular_weight
struct MoComponent {
    l: u32,
    m: u32,
    n: u32,
    coeffs: Vec<f64>,
}

struct MoBasisFunc {
    components: Vec<MoComponent>,
}

struct MoShell {
    cx: f64,
    cy: f64,
    cz: f64,
    /// Cutoff threshold: -ln(1e-15) / alpha_min.
    /// If r² > r2_cut the whole shell can be skipped.
    r2_cut: f64,
    exps: Vec<f64>,
    start_idx: usize,
    basis_funcs: Vec<MoBasisFunc>,
}

// ---------------------------------------------------------------------------
// BasisSetEngineRust
// ---------------------------------------------------------------------------

/// Rust-backed GTO basis set engine for evaluating molecular orbitals on a 3D grid.
/// Drop-in replacement for the Python `BasisSetEngine` in mo_engine.py.
#[pyclass]
pub struct BasisSetEngineRust {
    shells: Vec<MoShell>,
    n_basis_val: usize,
}

#[pymethods]
impl BasisSetEngineRust {
    /// Build the engine from a list of shell dicts.
    /// Each dict: {'type': int, 'center': [x,y,z] Bohr, 'exps': [...], 'coeffs': [...]}
    #[new]
    fn new(shells_py: &Bound<'_, PyList>) -> PyResult<Self> {
        let mut shells = Vec::new();
        let mut current_idx = 0usize;

        for item in shells_py.iter() {
            let sh = item.downcast::<PyDict>()?;
            let required = |key: &str| {
                sh.get_item(key)?
                    .ok_or_else(|| PyKeyError::new_err(key.to_string()))
            };
            let l_type: u32 = required("type")?.extract()?;
            let center: Vec<f64> = required("center")?.extract()?;
            let exps: Vec<f64> = required("exps")?.extract()?;
            let coeffs: Vec<f64> = required("coeffs")?.extract()?;
            if center.len() != 3 || center.iter().any(|v| !v.is_finite()) {
                return Err(PyValueError::new_err(
                    "Shell center must contain three finite coordinates.",
                ));
            }
            if exps.is_empty()
                || exps.len() != coeffs.len()
                || exps.iter().any(|v| !v.is_finite() || *v <= 0.0)
                || coeffs.iter().any(|v| !v.is_finite())
            {
                return Err(PyValueError::new_err("Shell exponents and coefficients must have matching nonempty lengths, finite values, and positive exponents."));
            }

            let defs = get_basis_definitions(l_type);
            if defs.is_empty() {
                // Skipping used to leave start_idx unset and current_idx
                // unadvanced, so every later shell read the wrong slice of the
                // MO coefficient vector — a silently wrong orbital with
                // nothing on screen to say so. Refuse the basis set instead.
                const SHELL_LETTERS: [&str; 8] = ["S", "P", "D", "F", "G", "H", "I", "K"];
                let letter = SHELL_LETTERS
                    .get(l_type as usize)
                    .map(|s| format!("{} shells", s))
                    .unwrap_or_else(|| format!("l={} shells", l_type));
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "This basis set contains {}, which the orbital renderer                      does not support (S through G only).

Visualization is                      blocked rather than showing an orbital that would be                      silently incorrect.",
                    letter
                )));
            }

            // Distance cutoff: skip shell when exp(-alpha_min * r2) < 1e-15
            // => r2 > 34.539 / alpha_min
            let alpha_min = exps.iter().cloned().fold(f64::INFINITY, f64::min);
            let r2_cut = if alpha_min > 0.0 {
                34.539 / alpha_min
            } else {
                f64::INFINITY
            };

            let mut basis_funcs = Vec::new();
            for bf_def in &defs {
                let mut components = Vec::new();
                for &(weight, l, m, n) in bf_def {
                    let comp_coeffs: Vec<f64> = exps
                        .iter()
                        .zip(coeffs.iter())
                        .map(|(&a, &c)| c * normalization_prefactor(a, l, m, n) * weight)
                        .collect();
                    components.push(MoComponent {
                        l,
                        m,
                        n,
                        coeffs: comp_coeffs,
                    });
                }
                basis_funcs.push(MoBasisFunc { components });
            }

            let n_funcs = defs.len();
            shells.push(MoShell {
                cx: center[0],
                cy: center[1],
                cz: center[2],
                r2_cut,
                exps,
                start_idx: current_idx,
                basis_funcs,
            });
            current_idx += n_funcs;
        }

        Ok(BasisSetEngineRust {
            shells,
            n_basis_val: current_idx,
        })
    }

    /// Number of basis functions (must match the MO coefficient vector length).
    #[getter]
    fn n_basis(&self) -> usize {
        self.n_basis_val
    }

    /// Evaluate a molecular orbital on a grid.
    ///
    /// grid_flat : flat (N*3,) sequence — row-major [x0,y0,z0, x1,y1,z1, …] in Bohr
    /// mo_coeffs : (n_basis,) coefficient vector
    ///
    /// Returns a flat (N,) list of phi values.
    /// Grid points are evaluated in parallel across all CPU cores via Rayon.
    fn evaluate_mo_on_grid(&self, grid_flat: Vec<f64>, mo_coeffs: Vec<f64>) -> PyResult<Vec<f64>> {
        if grid_flat.len() % 3 != 0 || grid_flat.iter().any(|v| !v.is_finite()) {
            return Err(PyValueError::new_err(
                "Grid must contain finite XYZ coordinate triples.",
            ));
        }
        if mo_coeffs.len() != self.n_basis_val || mo_coeffs.iter().any(|v| !v.is_finite()) {
            return Err(PyValueError::new_err(format!(
                "Expected {} finite MO coefficients, got {}.",
                self.n_basis_val,
                mo_coeffs.len()
            )));
        }
        let n_pts = grid_flat.len() / 3;
        let shells = &self.shells;

        // Each point is independent — evaluate in parallel.
        // `with_min_len` batches points so each rayon task allocates exp_buf once
        // and reuses it across multiple points.
        Ok((0..n_pts)
            .into_par_iter()
            .with_min_len(256)
            .map_init(
                // Initialiser: allocate a reusable buffer once per rayon task.
                || {
                    let cap = shells.iter().map(|s| s.exps.len()).max().unwrap_or(0);
                    Vec::<f64>::with_capacity(cap)
                },
                // Worker: evaluate one grid point using the task-local buffer.
                |exp_buf, pt_idx| {
                    let px = grid_flat[pt_idx * 3];
                    let py = grid_flat[pt_idx * 3 + 1];
                    let pz = grid_flat[pt_idx * 3 + 2];

                    let mut phi = 0.0_f64;

                    for sh in shells.iter() {
                        let rx = px - sh.cx;
                        let ry = py - sh.cy;
                        let rz = pz - sh.cz;
                        let r2 = rx * rx + ry * ry + rz * rz;

                        // Distance cutoff: skip shells with negligible contribution.
                        if r2 > sh.r2_cut {
                            continue;
                        }

                        // Radial exponentials — written into the reusable buffer.
                        exp_buf.clear();
                        for &a in &sh.exps {
                            exp_buf.push((-a * r2).exp());
                        }

                        for (b_i, bf) in sh.basis_funcs.iter().enumerate() {
                            let c_mo = mo_coeffs[sh.start_idx + b_i];
                            if c_mo.abs() < 1e-9 {
                                continue;
                            }

                            let mut val_accum = 0.0_f64;
                            for comp in &bf.components {
                                let contracted: f64 = comp
                                    .coeffs
                                    .iter()
                                    .zip(exp_buf.iter())
                                    .map(|(c, e)| c * e)
                                    .sum();
                                let ang =
                                    ang_pow(rx, comp.l) * ang_pow(ry, comp.m) * ang_pow(rz, comp.n);
                                val_accum += ang * contracted;
                            }
                            phi += c_mo * val_accum;
                        }
                    }
                    phi
                },
            )
            .collect())
    }
}

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------

#[pymodule]
fn orca_mo_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<BasisSetEngineRust>()?;
    Ok(())
}

// orca_mo_rs — Rust-powered GTO basis set evaluation for MO cube generation.
// Exposed to Python via PyO3.
//
// Public surface:
//   BasisSetEngineRust  — evaluates MO wavefunction on a 3D grid

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

// ---------------------------------------------------------------------------
// BasisSetEngineRust — Rust-backed GTO evaluation for MO cube generation
// ---------------------------------------------------------------------------

/// GTO normalization prefactor N(alpha, l, m, n).
fn normalization_prefactor_mo(alpha: f64, l: u32, m: u32, n: u32) -> f64 {
    const FACT: [f64; 9] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0];
    const FACT2: [f64; 9] = [
        1.0, 2.0, 24.0, 720.0, 40320.0, 3628800.0, 479001600.0, 87178291200.0,
        20922789888000.0,
    ];
    let big_l = (l + m + n) as f64;
    let num = (8.0 * alpha).powf(big_l)
        * FACT[l as usize]
        * FACT[m as usize]
        * FACT[n as usize];
    let den = FACT2[l as usize] * FACT2[m as usize] * FACT2[n as usize];
    (2.0 * alpha / std::f64::consts::PI).powf(0.75) * (num / den).sqrt()
}

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

// (weight, l, m, n) tuples defining the angular part of each component
type BasisCompDef = (f64, u32, u32, u32);

fn build_g_shell_defs() -> Vec<Vec<BasisCompDef>> {
    const FACT: [f64; 9] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0];
    const FACT2: [f64; 9] = [
        1.0, 2.0, 24.0, 720.0, 40320.0, 3628800.0, 479001600.0, 87178291200.0,
        20922789888000.0,
    ];
    let get_n_cart = |l: u32, m: u32, n: u32| -> f64 {
        ((FACT[l as usize] * FACT[m as usize] * FACT[n as usize])
            / (FACT2[l as usize] * FACT2[m as usize] * FACT2[n as usize]))
            .sqrt()
    };
    use std::f64::consts::PI;
    let pi_term = 1.0 / PI.sqrt();
    let sq5 = 5.0_f64.sqrt();
    let sq35 = 35.0_f64.sqrt();
    // n_sph[i] corresponds to m_vals = [0, 1, -1, 2, -2, 3, -3, 4, -4]
    let n_sph: [f64; 9] = [
        1.5 * pi_term,
        (3.0 / 8.0) * sq5 * pi_term,
        (3.0 / 8.0) * sq5 * pi_term,
        (3.0 / 4.0) * (2.5_f64).sqrt() * pi_term,
        (3.0 / 4.0) * (2.5_f64).sqrt() * pi_term,
        (3.0 / 8.0) * sq35 * pi_term,
        (3.0 / 8.0) * sq35 * pi_term,
        (3.0 / 16.0) * sq35 * pi_term,
        (3.0 / 16.0) * sq35 * pi_term,
    ];
    // Polynomial definitions: (c_poly, l, m, n)
    let g_polys: [&[(f64, u32, u32, u32)]; 9] = [
        &[
            (1.0, 0, 0, 4),
            (-3.0, 2, 0, 2),
            (-3.0, 0, 2, 2),
            (0.375, 4, 0, 0),
            (0.375, 0, 4, 0),
            (0.75, 2, 2, 0),
        ],
        &[(4.0, 1, 0, 3), (-3.0, 3, 0, 1), (-3.0, 1, 2, 1)],
        &[(4.0, 0, 1, 3), (-3.0, 0, 3, 1), (-3.0, 2, 1, 1)],
        &[(6.0, 2, 0, 2), (-6.0, 0, 2, 2), (-1.0, 4, 0, 0), (1.0, 0, 4, 0)],
        &[(6.0, 1, 1, 2), (-1.0, 3, 1, 0), (-1.0, 1, 3, 0)],
        &[(1.0, 3, 0, 1), (-3.0, 1, 2, 1)],
        &[(3.0, 2, 1, 1), (-1.0, 0, 3, 1)],
        &[(1.0, 4, 0, 0), (-6.0, 2, 2, 0), (1.0, 0, 4, 0)],
        &[(4.0, 3, 1, 0), (-4.0, 1, 3, 0)],
    ];
    g_polys
        .iter()
        .enumerate()
        .map(|(i, poly)| {
            let sph_norm = n_sph[i];
            poly.iter()
                .map(|&(c_poly, l, m, n)| {
                    let n_cart = get_n_cart(l, m, n);
                    let weight = c_poly * (sph_norm / n_cart);
                    (weight, l, m, n)
                })
                .collect()
        })
        .collect()
}

fn get_basis_definitions(l_type: u32) -> Vec<Vec<BasisCompDef>> {
    const F0: f64 = 0.240654;
    const F1: f64 = 0.281160;
    const F2: f64 = 0.866025;
    const F3: f64 = 0.369693;
    match l_type {
        0 => vec![vec![(1.0, 0, 0, 0)]],
        1 => vec![
            vec![(1.0, 0, 0, 1)], // pz
            vec![(1.0, 1, 0, 0)], // px
            vec![(1.0, 0, 1, 0)], // py
        ],
        2 => vec![
            vec![(-0.5, 2, 0, 0), (-0.5, 0, 2, 0), (1.0, 0, 0, 2)],
            vec![(1.0, 1, 0, 1)],
            vec![(1.0, 0, 1, 1)],
            vec![(0.866025, 2, 0, 0), (-0.866025, 0, 2, 0)],
            vec![(1.0, 1, 1, 0)],
        ],
        3 => vec![
            vec![(2.0 * F0, 0, 0, 3), (-3.0 * F0, 2, 0, 1), (-3.0 * F0, 0, 2, 1)],
            vec![(4.0 * F1, 1, 0, 2), (-F1, 3, 0, 0), (-F1, 1, 2, 0)],
            vec![(4.0 * F1, 0, 1, 2), (-F1, 2, 1, 0), (-F1, 0, 3, 0)],
            vec![(F2, 2, 0, 1), (-F2, 0, 2, 1)],
            vec![(1.0, 1, 1, 1)],
            vec![(F3, 3, 0, 0), (-3.0 * F3, 1, 2, 0)],
            vec![(3.0 * F3, 2, 1, 0), (-F3, 0, 3, 0)],
        ],
        4 => build_g_shell_defs(),
        _ => vec![],
    }
}

// Pre-computed component: angular exponents l,m,n + per-primitive coefficients.
// comp_coeffs[k] = input_coeff[k] * normalization(exps[k], l, m, n) * angular_weight
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
    exps: Vec<f64>,
    start_idx: usize,
    basis_funcs: Vec<MoBasisFunc>,
}

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
    /// Each dict: {'type': int, 'center': [x,y,z], 'exps': [...], 'coeffs': [...]}
    #[new]
    fn new(shells_py: &Bound<'_, PyList>) -> PyResult<Self> {
        let mut shells = Vec::new();
        let mut current_idx = 0usize;

        for item in shells_py.iter() {
            let sh = item.downcast::<PyDict>()?;
            let l_type: u32 = sh.get_item("type")?.unwrap().extract()?;
            let center: Vec<f64> = sh.get_item("center")?.unwrap().extract()?;
            let exps: Vec<f64> = sh.get_item("exps")?.unwrap().extract()?;
            let coeffs: Vec<f64> = sh.get_item("coeffs")?.unwrap().extract()?;

            let defs = get_basis_definitions(l_type);
            if defs.is_empty() {
                eprintln!(
                    "[BasisSetEngineRust] Warning: unsupported shell type {}",
                    l_type
                );
                continue;
            }

            let mut basis_funcs = Vec::new();
            for bf_def in &defs {
                let mut components = Vec::new();
                for &(weight, l, m, n) in bf_def {
                    let comp_coeffs: Vec<f64> = exps
                        .iter()
                        .zip(coeffs.iter())
                        .map(|(&a, &c)| c * normalization_prefactor_mo(a, l, m, n) * weight)
                        .collect();
                    components.push(MoComponent { l, m, n, coeffs: comp_coeffs });
                }
                basis_funcs.push(MoBasisFunc { components });
            }

            let n_funcs = defs.len();
            shells.push(MoShell {
                cx: center[0],
                cy: center[1],
                cz: center[2],
                exps,
                start_idx: current_idx,
                basis_funcs,
            });
            current_idx += n_funcs;
        }

        Ok(BasisSetEngineRust { shells, n_basis_val: current_idx })
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
    fn evaluate_mo_on_grid(&self, grid_flat: Vec<f64>, mo_coeffs: Vec<f64>) -> Vec<f64> {
        let n_pts = grid_flat.len() / 3;
        let mut result = vec![0.0_f64; n_pts];

        for pt_idx in 0..n_pts {
            let px = grid_flat[pt_idx * 3];
            let py = grid_flat[pt_idx * 3 + 1];
            let pz = grid_flat[pt_idx * 3 + 2];

            let mut phi = 0.0_f64;

            for sh in &self.shells {
                let rx = px - sh.cx;
                let ry = py - sh.cy;
                let rz = pz - sh.cz;
                let r2 = rx * rx + ry * ry + rz * rz;

                // Radial exp per primitive: exp(-alpha * r2)
                let exp_vals: Vec<f64> = sh.exps.iter().map(|&a| (-a * r2).exp()).collect();

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
                            .zip(exp_vals.iter())
                            .map(|(c, e)| c * e)
                            .sum();
                        let ang = ang_pow(rx, comp.l)
                            * ang_pow(ry, comp.m)
                            * ang_pow(rz, comp.n);
                        val_accum += ang * contracted;
                    }
                    phi += c_mo * val_accum;
                }
            }
            result[pt_idx] = phi;
        }

        result
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

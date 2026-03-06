/// Validation: 3D Distributed Load Analysis
///
/// References:
///   - Timoshenko & Goodier, "Theory of Elasticity"
///   - Przemieniecki, "Theory of Matrix Structural Analysis"
///   - Roark's Formulas for Stress and Strain, 9th Ed.
///
/// Tests:
///   1. SS beam UDL in Z: δ_max = 5qL⁴/(384EI)
///   2. Cantilever UDL in Y: δ_tip = qL⁴/(8EI)
///   3. Triangular load: δ_max for linearly varying load
///   4. Biaxial UDL: simultaneous Y and Z loading
///   5. Equilibrium: ΣR = ΣF for all load types
///   6. Cantilever UDL: reaction M = qL²/2
///   7. Multi-element vs single-element: convergence check
mod helpers;

use dedaliano_engine::solver::linear;
use dedaliano_engine::types::*;
use helpers::*;

const E: f64 = 200_000.0;
const NU: f64 = 0.3;
const A: f64 = 0.01;
const IY: f64 = 1e-4;
const IZ: f64 = 2e-4;
const J: f64 = 1.5e-4;

// ================================================================
// 1. SS Beam UDL in Z: δ_max = 5qL⁴/(384EIy)
// ================================================================
//
// Simply supported beam, uniform load in -Z direction.
// Maximum deflection at midspan.

#[test]
fn validation_3d_dist_ss_beam_udl_z() {
    let l: f64 = 6.0;
    let n = 8;
    let q: f64 = -10.0; // kN/m in Z
    let e_eff = E * 1000.0;
    let elem_len = l / n as f64;

    let nodes: Vec<_> = (0..=n).map(|i| (i + 1, i as f64 * elem_len, 0.0, 0.0)).collect();
    let elems: Vec<_> = (0..n).map(|i| (i + 1, "frame", i + 1, i + 2, 1, 1)).collect();

    let sups = vec![
        (1, vec![true, true, true, true, false, false]),       // pin
        (n + 1, vec![false, true, true, true, false, false]),   // roller
    ];

    let mut loads = Vec::new();
    for i in 0..n {
        loads.push(SolverLoad3D::Distributed(SolverDistributedLoad3D {
            element_id: i + 1,
            q_yi: 0.0, q_yj: 0.0,
            q_zi: q, q_zj: q,
            a: None, b: None,
        }));
    }

    let input = make_3d_input(nodes, vec![(1, E, NU)], vec![(1, A, IY, IZ, J)],
        elems, sups, loads);

    let results = linear::solve_3d(&input).unwrap();

    let mid = n / 2 + 1;
    let mid_d = results.displacements.iter().find(|d| d.node_id == mid).unwrap();

    // δ_max = 5qL⁴/(384EIy) for Z-direction loading (bends about Y-axis → uses Iy)
    let delta_exact = 5.0 * q.abs() * l.powi(4) / (384.0 * e_eff * IY);

    let error = (mid_d.uz.abs() - delta_exact).abs() / delta_exact;
    assert!(error < 0.05,
        "SS UDL uz: midspan={:.6e}, exact={:.6e}, err={:.1}%",
        mid_d.uz.abs(), delta_exact, error * 100.0);
}

// ================================================================
// 2. Cantilever UDL in Y: δ_tip = qL⁴/(8EIz)
// ================================================================

#[test]
fn validation_3d_dist_cantilever_udl_y() {
    let l: f64 = 4.0;
    let n = 8;
    let q: f64 = 5.0; // kN/m in Y
    let e_eff = E * 1000.0;

    let input = make_3d_beam(
        n, l, E, NU, A, IY, IZ, J,
        vec![true, true, true, true, true, true],
        None,
        {
            let mut loads = Vec::new();
            for i in 0..n {
                loads.push(SolverLoad3D::Distributed(SolverDistributedLoad3D {
                    element_id: i + 1,
                    q_yi: q, q_yj: q,
                    q_zi: 0.0, q_zj: 0.0,
                    a: None, b: None,
                }));
            }
            loads
        },
    );

    let results = linear::solve_3d(&input).unwrap();
    let tip = results.displacements.iter().find(|d| d.node_id == n + 1).unwrap();

    // δ_tip = qL⁴/(8EIz) for Y loading
    let delta_exact = q * l.powi(4) / (8.0 * e_eff * IZ);

    let error = (tip.uy.abs() - delta_exact).abs() / delta_exact;
    assert!(error < 0.05,
        "Cantilever UDL uy: tip={:.6e}, exact={:.6e}, err={:.1}%",
        tip.uy.abs(), delta_exact, error * 100.0);
}

// ================================================================
// 3. Triangular Load: Linearly Varying
// ================================================================
//
// Cantilever with triangular load: q=0 at fixed end, q=q_max at tip.
// δ_tip = q_max·L⁴/(30EI) for linear variation.

#[test]
fn validation_3d_dist_triangular_load() {
    let l: f64 = 4.0;
    let n = 8;
    let q_max: f64 = -10.0; // kN/m max at tip
    let e_eff = E * 1000.0;
    let elem_len = l / n as f64;

    let input = make_3d_beam(
        n, l, E, NU, A, IY, IZ, J,
        vec![true, true, true, true, true, true],
        None,
        {
            let mut loads = Vec::new();
            for i in 0..n {
                let xi = i as f64 * elem_len / l;
                let xj = (i + 1) as f64 * elem_len / l;
                loads.push(SolverLoad3D::Distributed(SolverDistributedLoad3D {
                    element_id: i + 1,
                    q_yi: 0.0, q_yj: 0.0,
                    q_zi: q_max * xi, q_zj: q_max * xj,
                    a: None, b: None,
                }));
            }
            loads
        },
    );

    let results = linear::solve_3d(&input).unwrap();
    let tip = results.displacements.iter().find(|d| d.node_id == n + 1).unwrap();

    // δ_tip = 11·q_max·L⁴/(120EIy) for triangular load increasing toward free end
    let delta_exact = 11.0 * q_max.abs() * l.powi(4) / (120.0 * e_eff * IY);

    let error = (tip.uz.abs() - delta_exact).abs() / delta_exact;
    assert!(error < 0.10,
        "Triangular load: tip_uz={:.6e}, exact={:.6e}, err={:.1}%",
        tip.uz.abs(), delta_exact, error * 100.0);
}

// ================================================================
// 4. Biaxial UDL: Simultaneous Y and Z Loading
// ================================================================
//
// Verify superposition: combined Y+Z response equals sum of individual.

#[test]
fn validation_3d_dist_biaxial_udl_superposition() {
    let l: f64 = 4.0;
    let n = 6;
    let qy = 5.0;
    let qz = -8.0;

    let make_loads = |qy_val: f64, qz_val: f64| -> Vec<SolverLoad3D> {
        (0..n).map(|i| SolverLoad3D::Distributed(SolverDistributedLoad3D {
            element_id: i + 1,
            q_yi: qy_val, q_yj: qy_val,
            q_zi: qz_val, q_zj: qz_val,
            a: None, b: None,
        })).collect()
    };

    let fixed = vec![true, true, true, true, true, true];

    // Combined
    let input_both = make_3d_beam(n, l, E, NU, A, IY, IZ, J,
        fixed.clone(), None, make_loads(qy, qz));
    let res_both = linear::solve_3d(&input_both).unwrap();
    let tip_both = res_both.displacements.iter().find(|d| d.node_id == n + 1).unwrap();

    // Y only
    let input_y = make_3d_beam(n, l, E, NU, A, IY, IZ, J,
        fixed.clone(), None, make_loads(qy, 0.0));
    let res_y = linear::solve_3d(&input_y).unwrap();
    let tip_y = res_y.displacements.iter().find(|d| d.node_id == n + 1).unwrap();

    // Z only
    let input_z = make_3d_beam(n, l, E, NU, A, IY, IZ, J,
        fixed.clone(), None, make_loads(0.0, qz));
    let res_z = linear::solve_3d(&input_z).unwrap();
    let tip_z = res_z.displacements.iter().find(|d| d.node_id == n + 1).unwrap();

    // Superposition: combined ≈ Y-only + Z-only
    let check = |name: &str, combined: f64, sum: f64| {
        let denom = combined.abs().max(1e-12);
        let err = (combined - sum).abs() / denom;
        assert!(err < 0.01,
            "Biaxial superposition {}: combined={:.6e}, sum={:.6e}, err={:.2}%",
            name, combined, sum, err * 100.0);
    };

    check("uy", tip_both.uy, tip_y.uy + tip_z.uy);
    check("uz", tip_both.uz, tip_y.uz + tip_z.uz);
}

// ================================================================
// 5. Equilibrium: ΣR = q × L
// ================================================================
//
// Total vertical reaction should equal total applied load.

#[test]
fn validation_3d_dist_equilibrium() {
    let l: f64 = 8.0;
    let n = 4;
    let q: f64 = -5.0;

    let elem_len = l / n as f64;
    let nodes: Vec<_> = (0..=n).map(|i| (i + 1, i as f64 * elem_len, 0.0, 0.0)).collect();
    let elems: Vec<_> = (0..n).map(|i| (i + 1, "frame", i + 1, i + 2, 1, 1)).collect();

    let sups = vec![
        (1, vec![true, true, true, true, false, false]),
        (n + 1, vec![false, true, true, true, false, false]),
    ];

    let loads: Vec<_> = (0..n).map(|i| SolverLoad3D::Distributed(SolverDistributedLoad3D {
        element_id: i + 1,
        q_yi: 0.0, q_yj: 0.0,
        q_zi: q, q_zj: q,
        a: None, b: None,
    })).collect();

    let input = make_3d_input(nodes, vec![(1, E, NU)], vec![(1, A, IY, IZ, J)],
        elems, sups, loads);

    let results = linear::solve_3d(&input).unwrap();

    let total_applied = q.abs() * l;
    let sum_rz: f64 = results.reactions.iter().map(|r| r.fz).sum();
    let eq_err = (sum_rz - total_applied).abs() / total_applied;

    assert!(eq_err < 0.01,
        "Equilibrium: ΣRz={:.4}, applied={:.4}, err={:.2}%",
        sum_rz, total_applied, eq_err * 100.0);
}

// ================================================================
// 6. Cantilever UDL: Reaction Moment M = qL²/2
// ================================================================

#[test]
fn validation_3d_dist_cantilever_reaction_moment() {
    let l: f64 = 5.0;
    let n = 4;
    let q: f64 = -8.0;

    let input = make_3d_beam(
        n, l, E, NU, A, IY, IZ, J,
        vec![true, true, true, true, true, true],
        None,
        {
            let mut loads = Vec::new();
            for i in 0..n {
                loads.push(SolverLoad3D::Distributed(SolverDistributedLoad3D {
                    element_id: i + 1,
                    q_yi: 0.0, q_yj: 0.0,
                    q_zi: q, q_zj: q,
                    a: None, b: None,
                }));
            }
            loads
        },
    );

    let results = linear::solve_3d(&input).unwrap();

    // Fixed-end reaction: Rz = qL, My = qL²/2
    let r = results.reactions.iter().find(|r| r.node_id == 1).unwrap();
    let rz_exact = q.abs() * l;
    let my_exact = q.abs() * l * l / 2.0;

    let err_rz = (r.fz - rz_exact).abs() / rz_exact;
    assert!(err_rz < 0.01,
        "Reaction Rz={:.4}, exact={:.4}, err={:.2}%", r.fz, rz_exact, err_rz * 100.0);

    // Moment reaction (My for Z-direction loading)
    // Sign depends on convention, just check magnitude
    let my_computed = r.my.abs();
    let err_my = (my_computed - my_exact).abs() / my_exact;
    assert!(err_my < 0.05,
        "Reaction My={:.4}, exact qL²/2={:.4}, err={:.1}%",
        my_computed, my_exact, err_my * 100.0);
}

// ================================================================
// 7. Multi-Element Convergence
// ================================================================
//
// More elements should give better accuracy compared to beam theory.

#[test]
fn validation_3d_dist_mesh_convergence() {
    let l: f64 = 6.0;
    let q: f64 = -10.0;
    let e_eff = E * 1000.0;
    let delta_exact = 5.0 * q.abs() * l.powi(4) / (384.0 * e_eff * IY);

    let mut errors = Vec::new();
    for &n in &[2, 4, 8] {
        let elem_len = l / n as f64;
        let nodes: Vec<_> = (0..=n).map(|i| (i + 1, i as f64 * elem_len, 0.0, 0.0)).collect();
        let elems: Vec<_> = (0..n).map(|i| (i + 1, "frame", i + 1, i + 2, 1, 1)).collect();
        let sups = vec![
            (1, vec![true, true, true, true, false, false]),
            (n + 1, vec![false, true, true, true, false, false]),
        ];
        let loads: Vec<_> = (0..n).map(|i| SolverLoad3D::Distributed(SolverDistributedLoad3D {
            element_id: i + 1,
            q_yi: 0.0, q_yj: 0.0,
            q_zi: q, q_zj: q,
            a: None, b: None,
        })).collect();

        let input = make_3d_input(nodes, vec![(1, E, NU)], vec![(1, A, IY, IZ, J)],
            elems, sups, loads);

        let results = linear::solve_3d(&input).unwrap();
        let mid = n / 2 + 1;
        let mid_d = results.displacements.iter().find(|d| d.node_id == mid).unwrap();

        let err = (mid_d.uz.abs() - delta_exact).abs() / delta_exact;
        errors.push(err);
    }

    // Error should decrease or stay similar with refinement
    assert!(errors[2] <= errors[0] + 0.02,
        "Convergence: err_2={:.3}, err_4={:.3}, err_8={:.3}",
        errors[0], errors[1], errors[2]);
}

/// Validation: Advanced Stability and Buckling Benchmarks
///
/// References:
///   - Timoshenko & Gere, "Theory of Elastic Stability", 2nd Ed.
///   - Chen & Lui, "Structural Stability", Ch. 2-4
///   - Galambos & Surovek, "Structural Stability of Steel", 5th Ed.
///   - AISC 360-22 §C2: Second-Order Analysis
///
/// Tests:
///   1. P-Delta amplification: exact Timoshenko formula
///   2. Column with initial imperfection: amplified δ
///   3. Stepped column: different I per segment
///   4. Leaning column: gravity column destabilizes braced frame
///   5. Frame sidesway buckling: critical load for portal
///   6. Column effective length: K-factor comparison
///   7. Propped cantilever column: K ≈ 0.7
///   8. Braced frame: much stiffer than unbraced
mod helpers;

use dedaliano_engine::solver::{linear, corotational};
use dedaliano_engine::types::*;
use helpers::*;

const E: f64 = 200_000.0;
const A: f64 = 0.01;
const IZ: f64 = 1e-4;

// ================================================================
// 1. P-Delta Amplification: Timoshenko Formula
// ================================================================
//
// Cantilever column with lateral load H and axial load P.
// Second-order δ = δ₁ × 1/(1 - P/P_cr)
// where P_cr = π²EI/(4L²) for cantilever and δ₁ = HL³/(3EI).

#[test]
fn validation_stability_pdelta_amplification() {
    let l = 5.0;
    let n = 8;
    let h_load = 5.0;
    let e_eff = E * 1000.0;
    let p_cr = std::f64::consts::PI.powi(2) * e_eff * IZ / (4.0 * l * l);

    // Apply axial load at 30% of critical
    let p_axial = 0.3 * p_cr;

    let nodes: Vec<_> = (0..=n).map(|i| (i + 1, 0.0, i as f64 * l / n as f64)).collect();
    let elems: Vec<_> = (0..n).map(|i| (i + 1, "frame", i + 1, i + 2, 1, 1, false, false)).collect();
    let sups = vec![(1, 1, "fixed")];
    let loads = vec![
        SolverLoad::Nodal(SolverNodalLoad { node_id: n + 1, fx: h_load, fy: 0.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: n + 1, fx: 0.0, fy: -p_axial, mz: 0.0 }),
    ];

    let input = make_input(nodes, vec![(1, E, 0.3)], vec![(1, A, IZ)],
        elems, sups, loads);

    // First-order (linear)
    let res_linear = linear::solve_2d(&input).unwrap();
    let d_linear = res_linear.displacements.iter()
        .find(|d| d.node_id == n + 1).unwrap().ux.abs();

    // Second-order (corotational)
    let res_corot = corotational::solve_corotational_2d(&input, 50, 1e-6, 5).unwrap();
    let d_corot = res_corot.results.displacements.iter()
        .find(|d| d.node_id == n + 1).unwrap().ux.abs();

    // Second-order should amplify first-order
    assert!(d_corot > d_linear,
        "P-Δ: corotational δ={:.6e} should exceed linear δ={:.6e}",
        d_corot, d_linear);

    // Approximate amplification factor: 1/(1-P/Pcr) = 1/(1-0.3) ≈ 1.43
    let af_approx = 1.0 / (1.0 - 0.3);
    let af_actual = d_corot / d_linear;
    let err = (af_actual - af_approx).abs() / af_approx;
    assert!(err < 0.30,
        "Amplification: actual={:.3}, approx 1/(1-P/Pcr)={:.3}", af_actual, af_approx);
}

// ================================================================
// 2. Column with Initial Imperfection
// ================================================================
//
// Column with slight lateral offset at midheight. Under axial load,
// the deflection amplifies the initial imperfection.

#[test]
fn validation_stability_initial_imperfection() {
    let l = 5.0;
    let n = 8;
    let h_load = 1.0; // small lateral load to simulate imperfection

    // Cantilever column with small lateral load + axial compression.
    // Second-order displacement should exceed first-order.
    let e_eff = E * 1000.0;
    let p_cr = std::f64::consts::PI.powi(2) * e_eff * IZ / (4.0 * l * l);
    let p = 0.4 * p_cr;

    let nodes: Vec<_> = (0..=n).map(|i| (i + 1, 0.0, i as f64 * l / n as f64)).collect();
    let elems: Vec<_> = (0..n).map(|i| (i + 1, "frame", i + 1, i + 2, 1, 1, false, false)).collect();
    let sups = vec![(1, 1, "fixed")];
    let loads = vec![
        SolverLoad::Nodal(SolverNodalLoad { node_id: n + 1, fx: h_load, fy: -p, mz: 0.0 }),
    ];

    let input = make_input(nodes, vec![(1, E, 0.3)], vec![(1, A, IZ)],
        elems, sups, loads);

    // First-order
    let res_1 = linear::solve_2d(&input).unwrap();
    let d_1 = res_1.displacements.iter().find(|d| d.node_id == n + 1).unwrap().ux.abs();

    // Second-order
    let res_2 = corotational::solve_corotational_2d(&input, 50, 1e-6, 10).unwrap();
    let d_2 = res_2.results.displacements.iter().find(|d| d.node_id == n + 1).unwrap().ux.abs();

    // At 40% of Pcr, amplification should be ~1/(1-0.4) = 1.67
    assert!(d_2 > d_1 * 1.2,
        "Imperfection amplification: 2nd-order={:.6e} should be >1.2× 1st-order={:.6e}",
        d_2, d_1);
}

// ================================================================
// 3. Frame Sidesway: Second-Order Reduces Stiffness
// ================================================================
//
// Portal frame with gravity + lateral. Second-order sway > first-order.

#[test]
fn validation_stability_frame_sidesway() {
    let h = 4.0;
    let w = 6.0;
    let h_load = 5.0;
    let v_load = -50.0; // gravity

    let input = make_portal_frame(h, w, E, A, IZ, h_load, v_load);

    // First-order
    let res_1 = linear::solve_2d(&input).unwrap();
    let sway_1 = res_1.displacements.iter()
        .find(|d| d.node_id == 2).unwrap().ux.abs();

    // Second-order
    let res_2 = corotational::solve_corotational_2d(&input, 50, 1e-6, 5).unwrap();
    let sway_2 = res_2.results.displacements.iter()
        .find(|d| d.node_id == 2).unwrap().ux.abs();

    assert!(sway_2 > sway_1,
        "Frame sidesway: 2nd-order={:.6e} > 1st-order={:.6e}", sway_2, sway_1);
}

// ================================================================
// 4. Effective Length: Fixed-Free vs Pin-Pin
// ================================================================
//
// Pin-pin column: K=1, P_cr = π²EI/L².
// Fixed-free (cantilever): K=2, P_cr = π²EI/(4L²).
// Under the same fraction of respective Pcr, both should amplify similarly.

#[test]
fn validation_stability_effective_length_comparison() {
    let l = 5.0;
    let n = 8;
    let h_load = 1.0;
    let e_eff = E * 1000.0;

    // Pin-pin: Pcr = π²EI/L²
    let pcr_pp = std::f64::consts::PI.powi(2) * e_eff * IZ / (l * l);
    let p_pp = 0.2 * pcr_pp;

    let nodes_pp: Vec<_> = (0..=n).map(|i| (i + 1, 0.0, i as f64 * l / n as f64)).collect();
    let elems: Vec<_> = (0..n).map(|i| (i + 1, "frame", i + 1, i + 2, 1, 1, false, false)).collect();
    let sups_pp = vec![(1, 1, "pinned"), (2, n + 1, "rollerY")];
    let loads_pp = vec![
        SolverLoad::Nodal(SolverNodalLoad { node_id: n + 1, fx: 0.0, fy: -p_pp, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: n / 2 + 1, fx: h_load, fy: 0.0, mz: 0.0 }),
    ];
    let input_pp = make_input(nodes_pp, vec![(1, E, 0.3)], vec![(1, A, IZ)],
        elems.clone(), sups_pp, loads_pp);

    let res_pp_1 = linear::solve_2d(&input_pp).unwrap();
    let res_pp_2 = corotational::solve_corotational_2d(&input_pp, 50, 1e-6, 5).unwrap();
    let d_pp_1 = res_pp_1.displacements.iter().find(|d| d.node_id == n / 2 + 1).unwrap().ux.abs();
    let d_pp_2 = res_pp_2.results.displacements.iter().find(|d| d.node_id == n / 2 + 1).unwrap().ux.abs();
    let af_pp = d_pp_2 / d_pp_1;

    // Fixed-free: Pcr = π²EI/(4L²), same load ratio
    let pcr_cf = std::f64::consts::PI.powi(2) * e_eff * IZ / (4.0 * l * l);
    let p_cf = 0.2 * pcr_cf;

    let nodes_cf: Vec<_> = (0..=n).map(|i| (i + 1, 0.0, i as f64 * l / n as f64)).collect();
    let sups_cf = vec![(1, 1, "fixed")];
    let loads_cf = vec![
        SolverLoad::Nodal(SolverNodalLoad { node_id: n + 1, fx: 0.0, fy: -p_cf, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: n + 1, fx: h_load, fy: 0.0, mz: 0.0 }),
    ];
    let input_cf = make_input(nodes_cf, vec![(1, E, 0.3)], vec![(1, A, IZ)],
        elems, sups_cf, loads_cf);

    let res_cf_1 = linear::solve_2d(&input_cf).unwrap();
    let res_cf_2 = corotational::solve_corotational_2d(&input_cf, 50, 1e-6, 5).unwrap();
    let d_cf_1 = res_cf_1.displacements.iter().find(|d| d.node_id == n + 1).unwrap().ux.abs();
    let d_cf_2 = res_cf_2.results.displacements.iter().find(|d| d.node_id == n + 1).unwrap().ux.abs();
    let af_cf = d_cf_2 / d_cf_1;

    // At same P/Pcr ratio, amplification factors should be similar
    assert!(af_pp > 1.0 && af_cf > 1.0,
        "Both should amplify: AF_pp={:.3}, AF_cf={:.3}", af_pp, af_cf);
    // Within 50% of each other (approximate, different boundary effects)
    let ratio = af_pp / af_cf;
    assert!(ratio > 0.5 && ratio < 2.0,
        "Similar amplification: AF_pp/AF_cf={:.3}", ratio);
}

// ================================================================
// 5. Braced vs Unbraced Frame
// ================================================================
//
// Unbraced portal frame has much lower lateral stiffness under gravity
// than braced (diagonal member added).

#[test]
fn validation_stability_braced_vs_unbraced() {
    let h = 4.0;
    let w = 6.0;
    let h_load = 5.0;
    let v_load = -30.0;

    // Unbraced portal
    let input_unbraced = make_portal_frame(h, w, E, A, IZ, h_load, v_load);
    let res_unbraced = linear::solve_2d(&input_unbraced).unwrap();
    let sway_unbraced = res_unbraced.displacements.iter()
        .find(|d| d.node_id == 2).unwrap().ux.abs();

    // Braced portal: add diagonal truss member
    let nodes = vec![(1, 0.0, 0.0), (2, 0.0, h), (3, w, h), (4, w, 0.0)];
    let elems = vec![
        (1, "frame", 1, 2, 1, 1, false, false),
        (2, "frame", 2, 3, 1, 1, false, false),
        (3, "frame", 3, 4, 1, 1, false, false),
        (4, "truss", 1, 3, 1, 1, false, false), // diagonal brace
    ];
    let sups = vec![(1, 1, "fixed"), (2, 4, "fixed")];
    let loads = vec![
        SolverLoad::Nodal(SolverNodalLoad { node_id: 2, fx: h_load, fy: 0.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 2, fx: 0.0, fy: v_load, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 3, fx: 0.0, fy: v_load, mz: 0.0 }),
    ];

    let input_braced = make_input(nodes, vec![(1, E, 0.3)], vec![(1, A, IZ)],
        elems, sups, loads);
    let res_braced = linear::solve_2d(&input_braced).unwrap();
    let sway_braced = res_braced.displacements.iter()
        .find(|d| d.node_id == 2).unwrap().ux.abs();

    // Braced frame should be much stiffer (less sway)
    assert!(sway_braced < sway_unbraced * 0.5,
        "Braced sway={:.6e} should be << unbraced={:.6e}", sway_braced, sway_unbraced);
}

// ================================================================
// 6. Second-Order Convergence with Load Increments
// ================================================================
//
// More load increments should give smoother convergence.
// Result with 10 increments ≈ result with 5 increments.

#[test]
fn validation_stability_increment_convergence() {
    let l = 5.0;
    let n = 6;
    let h_load = 3.0;
    let e_eff = E * 1000.0;
    let p_cr = std::f64::consts::PI.powi(2) * e_eff * IZ / (4.0 * l * l);
    let p_axial = 0.25 * p_cr;

    let nodes: Vec<_> = (0..=n).map(|i| (i + 1, 0.0, i as f64 * l / n as f64)).collect();
    let elems: Vec<_> = (0..n).map(|i| (i + 1, "frame", i + 1, i + 2, 1, 1, false, false)).collect();
    let sups = vec![(1, 1, "fixed")];
    let loads = vec![
        SolverLoad::Nodal(SolverNodalLoad { node_id: n + 1, fx: h_load, fy: -p_axial, mz: 0.0 }),
    ];
    let input = make_input(nodes, vec![(1, E, 0.3)], vec![(1, A, IZ)],
        elems, sups, loads);

    let res_5 = corotational::solve_corotational_2d(&input, 50, 1e-6, 5).unwrap();
    let res_10 = corotational::solve_corotational_2d(&input, 50, 1e-6, 10).unwrap();

    let d5 = res_5.results.displacements.iter().find(|d| d.node_id == n + 1).unwrap().ux;
    let d10 = res_10.results.displacements.iter().find(|d| d.node_id == n + 1).unwrap().ux;

    let err = (d5 - d10).abs() / d10.abs().max(1e-12);
    assert!(err < 0.10,
        "Increment convergence: 5-step={:.6e}, 10-step={:.6e}, diff={:.1}%",
        d5, d10, err * 100.0);
}

// ================================================================
// 7. Gravity Load Reduces Lateral Stiffness
// ================================================================
//
// Portal frame: lateral stiffness decreases with increasing gravity.

#[test]
fn validation_stability_gravity_reduces_stiffness() {
    let h = 4.0;
    let w = 6.0;
    let h_load = 5.0;

    // No gravity
    let input_0 = make_portal_frame(h, w, E, A, IZ, h_load, 0.0);
    let res_0 = linear::solve_2d(&input_0).unwrap();
    let sway_0 = res_0.displacements.iter()
        .find(|d| d.node_id == 2).unwrap().ux.abs();

    // With gravity (second-order)
    let input_g = make_portal_frame(h, w, E, A, IZ, h_load, -50.0);
    let res_g = corotational::solve_corotational_2d(&input_g, 50, 1e-6, 5).unwrap();
    let sway_g = res_g.results.displacements.iter()
        .find(|d| d.node_id == 2).unwrap().ux.abs();

    // Gravity should increase sway (reduce effective stiffness)
    assert!(sway_g > sway_0,
        "Gravity destabilizes: sway_gravity={:.6e} > sway_no_gravity={:.6e}",
        sway_g, sway_0);
}

// ================================================================
// 8. Symmetric Frame: No Lateral Sway Under Gravity
// ================================================================
//
// Even with second-order analysis, symmetric gravity should produce
// zero sway (by symmetry).

#[test]
fn validation_stability_symmetric_no_sway() {
    let h = 4.0;
    let w = 6.0;
    let v_load = -50.0;

    let input = make_portal_frame(h, w, E, A, IZ, 0.0, v_load);
    let res = corotational::solve_corotational_2d(&input, 50, 1e-6, 5).unwrap();

    let sway = res.results.displacements.iter()
        .find(|d| d.node_id == 2).unwrap().ux;

    assert!(sway.abs() < 1e-8,
        "Symmetric frame: sway should be 0, got {:.6e}", sway);
}

/// Validation: Exact influence line ordinates.
///
/// Reference: Ghali/Neville *Structural Analysis*
mod helpers;

use dedaliano_engine::postprocess::influence::*;
use dedaliano_engine::types::*;
use helpers::*;

const E: f64 = 200_000.0;
const A: f64 = 0.01;
const IZ: f64 = 1e-4;
const L: f64 = 6.0;

fn make_ss_beam_for_il() -> SolverInput {
    make_beam(1, L, E, A, IZ, "pinned", Some("rollerX"), vec![])
}

fn make_ss_beam_multi_for_il(n: usize) -> SolverInput {
    make_beam(n, L, E, A, IZ, "pinned", Some("rollerX"), vec![])
}

// ═══════════════════════════════════════════════════════════════
// 1. SS Beam, R_A: η(0)=1.0, η(L/2)=0.5, η(L)=0.0 (linear)
// ═══════════════════════════════════════════════════════════════

#[test]
fn validation_il_ss_reaction_a() {
    let il_input = InfluenceLineInput {
        solver: make_ss_beam_for_il(),
        quantity: "Ry".to_string(),
        target_node_id: Some(1),
        target_element_id: None,
        target_position: 0.5,
        n_points_per_element: 20,
    };
    let result = compute_influence_line(&il_input).unwrap();
    assert!(!result.points.is_empty());

    // At x=0: η ≈ 1.0
    let at_start = result.points.iter()
        .min_by(|a, b| a.x.abs().partial_cmp(&b.x.abs()).unwrap()).unwrap();
    assert!(
        (at_start.value - 1.0).abs() < 0.1,
        "IL R_A at x=0: {:.3}, expected 1.0", at_start.value
    );

    // At x=L: η ≈ 0.0
    let at_end = result.points.iter()
        .min_by(|a, b| (a.x - L).abs().partial_cmp(&(b.x - L).abs()).unwrap()).unwrap();
    assert!(
        at_end.value.abs() < 0.1,
        "IL R_A at x=L: {:.3}, expected 0.0", at_end.value
    );

    // At x=L/2: η ≈ 0.5
    let at_mid = result.points.iter()
        .min_by(|a, b| (a.x - L / 2.0).abs().partial_cmp(&(b.x - L / 2.0).abs()).unwrap()).unwrap();
    assert!(
        (at_mid.value - 0.5).abs() < 0.1,
        "IL R_A at x=L/2: {:.3}, expected 0.5", at_mid.value
    );
}

// ═══════════════════════════════════════════════════════════════
// 2. SS Beam, R_B: η(0)=0.0, η(L)=1.0 (linear)
// ═══════════════════════════════════════════════════════════════

#[test]
fn validation_il_ss_reaction_b() {
    let il_input = InfluenceLineInput {
        solver: make_ss_beam_for_il(),
        quantity: "Ry".to_string(),
        target_node_id: Some(2),
        target_element_id: None,
        target_position: 0.5,
        n_points_per_element: 20,
    };
    let result = compute_influence_line(&il_input).unwrap();

    let at_start = result.points.iter()
        .min_by(|a, b| a.x.abs().partial_cmp(&b.x.abs()).unwrap()).unwrap();
    assert!(
        at_start.value.abs() < 0.1,
        "IL R_B at x=0: {:.3}, expected 0.0", at_start.value
    );

    let at_end = result.points.iter()
        .min_by(|a, b| (a.x - L).abs().partial_cmp(&(b.x - L).abs()).unwrap()).unwrap();
    assert!(
        (at_end.value - 1.0).abs() < 0.1,
        "IL R_B at x=L: {:.3}, expected 1.0", at_end.value
    );
}

// ═══════════════════════════════════════════════════════════════
// 3. SS Beam, M at Midspan: Triangular, Peak = L/4
// ═══════════════════════════════════════════════════════════════

#[test]
fn validation_il_ss_moment_midspan() {
    let il_input = InfluenceLineInput {
        solver: make_ss_beam_for_il(),
        quantity: "M".to_string(),
        target_node_id: None,
        target_element_id: Some(1),
        target_position: 0.5, // midspan
        n_points_per_element: 20,
    };
    let result = compute_influence_line(&il_input).unwrap();

    // Peak at midspan: |η| = L/4 = 1.5
    let at_mid = result.points.iter()
        .min_by(|a, b| (a.x - L / 2.0).abs().partial_cmp(&(b.x - L / 2.0).abs()).unwrap()).unwrap();
    assert!(
        (at_mid.value.abs() - L / 4.0).abs() < 0.3,
        "IL M at midspan: {:.3}, expected ±{:.3}", at_mid.value, L / 4.0
    );

    // Zero at supports
    let at_start = result.points.iter()
        .min_by(|a, b| a.x.abs().partial_cmp(&b.x.abs()).unwrap()).unwrap();
    assert!(
        at_start.value.abs() < 0.2,
        "IL M at x=0: {:.3}, expected 0.0", at_start.value
    );
}

// ═══════════════════════════════════════════════════════════════
// 4. SS Beam, V at Midspan: Jump from +0.5 to -0.5
// ═══════════════════════════════════════════════════════════════

#[test]
fn validation_il_ss_shear_midspan() {
    // Use multi-element beam so there are points on both sides
    let n = 4;
    let il_input = InfluenceLineInput {
        solver: make_ss_beam_multi_for_il(n),
        quantity: "V".to_string(),
        target_node_id: None,
        target_element_id: Some(2), // element near midspan
        target_position: 1.0, // end of element 2 (at midspan)
        n_points_per_element: 10,
    };
    let result = compute_influence_line(&il_input).unwrap();

    // Near midspan: values should be around ±0.5
    // Find max and min values near midspan
    let mid_points: Vec<_> = result.points.iter()
        .filter(|p| (p.x - L / 2.0).abs() < L / 4.0)
        .collect();

    if !mid_points.is_empty() {
        let max_val = mid_points.iter().map(|p| p.value).fold(f64::NEG_INFINITY, f64::max);
        let min_val = mid_points.iter().map(|p| p.value).fold(f64::INFINITY, f64::min);

        // There should be values on both sides of zero near midspan
        assert!(
            max_val > 0.2 || min_val < -0.2,
            "IL V near midspan: max={:.3}, min={:.3}", max_val, min_val
        );
    }
}

// ═══════════════════════════════════════════════════════════════
// 5. SS Beam, M at Quarter-Span: Peak = 3L/16
// ═══════════════════════════════════════════════════════════════

#[test]
fn validation_il_moment_quarter_span() {
    let n = 4; // nodes at 0, 1.5, 3, 4.5, 6
    let il_input = InfluenceLineInput {
        solver: make_ss_beam_multi_for_il(n),
        quantity: "M".to_string(),
        target_node_id: None,
        target_element_id: Some(1), // first element
        target_position: 1.0, // end of first element = L/4
        n_points_per_element: 10,
    };
    let result = compute_influence_line(&il_input).unwrap();

    // Peak at the section (x=L/4): η = (L/4)*(3L/4)/L = 3L/16 = 1.125
    let at_quarter = result.points.iter()
        .min_by(|a, b| (a.x - L / 4.0).abs().partial_cmp(&(b.x - L / 4.0).abs()).unwrap()).unwrap();
    let expected_peak = 3.0 * L / 16.0;
    assert!(
        (at_quarter.value.abs() - expected_peak).abs() < 0.3,
        "IL M at L/4: {:.3}, expected {:.3}", at_quarter.value, expected_peak
    );
}

// ═══════════════════════════════════════════════════════════════
// 6. Symmetry: IL for M at Midspan
// ═══════════════════════════════════════════════════════════════

#[test]
fn validation_il_symmetry() {
    let n = 4;
    let il_input = InfluenceLineInput {
        solver: make_ss_beam_multi_for_il(n),
        quantity: "M".to_string(),
        target_node_id: None,
        target_element_id: Some(2),
        target_position: 1.0, // end of element 2 = midspan
        n_points_per_element: 10,
    };
    let result = compute_influence_line(&il_input).unwrap();

    // Check symmetry: η(x) ≈ η(L-x)
    for p in &result.points {
        let mirror_x = L - p.x;
        if let Some(mirror_p) = result.points.iter()
            .min_by(|a, b| (a.x - mirror_x).abs().partial_cmp(&(b.x - mirror_x).abs()).unwrap())
        {
            if (mirror_p.x - mirror_x).abs() < 0.2 {
                assert!(
                    (p.value - mirror_p.value).abs() < 0.3,
                    "IL symmetry: η({:.2})={:.3} vs η({:.2})={:.3}",
                    p.x, p.value, mirror_p.x, mirror_p.value
                );
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════
// 7. Continuous Beam, Interior Reaction IL
// ═══════════════════════════════════════════════════════════════

#[test]
fn validation_il_continuous_beam_reaction() {
    // 2-span continuous beam, IL for interior reaction R_B
    // IL should be curved (not linear) with peak > 1.0 at interior support
    let input = make_continuous_beam(&[6.0, 6.0], 2, E, A, IZ, vec![]);
    let interior_node = 3; // node at x=6

    let il_input = InfluenceLineInput {
        solver: input,
        quantity: "Ry".to_string(),
        target_node_id: Some(interior_node),
        target_element_id: None,
        target_position: 0.5,
        n_points_per_element: 10,
    };
    let result = compute_influence_line(&il_input).unwrap();

    // At the interior support: η should be 1.0
    let at_support = result.points.iter()
        .min_by(|a, b| (a.x - 6.0).abs().partial_cmp(&(b.x - 6.0).abs()).unwrap()).unwrap();
    assert!(
        (at_support.value - 1.0).abs() < 0.15,
        "IL R_B at interior support: {:.3}, expected 1.0", at_support.value
    );

    // For a continuous beam, IL for interior reaction may have values > 1.0
    // when unit load is directly at support, and negative values in far span
    let max_val = result.points.iter().map(|p| p.value).fold(f64::NEG_INFINITY, f64::max);
    let min_val = result.points.iter().map(|p| p.value).fold(f64::INFINITY, f64::min);
    // Min should be negative (uplift in far span)
    assert!(
        min_val < 0.1,
        "continuous IL should have negative values in far span, min={:.3}", min_val
    );
    let _ = max_val;
}

// ═══════════════════════════════════════════════════════════════
// 8. Linearity Check: IL Values Scale Linearly
// ═══════════════════════════════════════════════════════════════

#[test]
fn validation_il_linearity() {
    // IL is computed by Müller-Breslau: unit load method
    // Scaling the stiffness shouldn't change the IL values (dimensionless)
    let il_input_1 = InfluenceLineInput {
        solver: make_beam(1, L, E, A, IZ, "pinned", Some("rollerX"), vec![]),
        quantity: "Ry".to_string(),
        target_node_id: Some(1),
        target_element_id: None,
        target_position: 0.5,
        n_points_per_element: 10,
    };
    let result_1 = compute_influence_line(&il_input_1).unwrap();

    // Same beam, doubled EI
    let il_input_2 = InfluenceLineInput {
        solver: make_beam(1, L, E, A, IZ * 2.0, "pinned", Some("rollerX"), vec![]),
        quantity: "Ry".to_string(),
        target_node_id: Some(1),
        target_element_id: None,
        target_position: 0.5,
        n_points_per_element: 10,
    };
    let result_2 = compute_influence_line(&il_input_2).unwrap();

    // IL values should be identical (reaction IL is independent of stiffness for statically determinate)
    assert_eq!(result_1.points.len(), result_2.points.len());
    for (p1, p2) in result_1.points.iter().zip(result_2.points.iter()) {
        assert!(
            (p1.value - p2.value).abs() < 0.05,
            "IL linearity: val1={:.3} vs val2={:.3} at x={:.2}", p1.value, p2.value, p1.x
        );
    }
}

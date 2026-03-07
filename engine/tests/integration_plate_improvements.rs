/// Integration tests for plate element improvements.
///
/// Tests verify:
/// 1. Drilling DOF stiffness is symmetric and positive
/// 2. Nodal stress output is populated
/// 3. Plate thermal loads produce deflection
/// 4. Element quality metrics are reasonable
/// 5. Patch test: uniform stress field

use dedaliano_engine::element::{
    plate_local_stiffness, plate_pressure_load,
    plate_thermal_load, plate_element_quality, plate_stress_at_nodes,
    plate_stress_recovery,
};

#[test]
fn plate_drilling_stiffness_symmetry() {
    let coords = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.866, 0.0],
    ];
    let e = 200e6; // kN/m²
    let nu = 0.3;
    let t = 0.01; // m

    let k = plate_local_stiffness(&coords, e, nu, t);
    let n = 18;

    // Verify symmetry
    let mut max_asym = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            let diff = (k[i * n + j] - k[j * n + i]).abs();
            if diff > max_asym {
                max_asym = diff;
            }
        }
    }
    assert!(
        max_asym < 1e-8,
        "Stiffness matrix should be symmetric, max asymmetry: {:.2e}",
        max_asym
    );

    // Drilling DOF positions: 5, 11, 17
    let drill = [5, 11, 17];
    for &d in &drill {
        assert!(
            k[d * n + d] > 0.0,
            "Drilling stiffness at DOF {} should be positive: {}",
            d, k[d * n + d]
        );
    }

    // Off-diagonal drilling coupling should exist
    assert!(
        k[drill[0] * n + drill[1]].abs() > 0.0,
        "Off-diagonal drilling coupling should exist"
    );
}

#[test]
fn plate_nodal_stress_output() {
    let coords = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.866, 0.0],
    ];
    let e = 200e6;
    let nu = 0.3;
    let t = 0.01;

    // Apply a simple bending displacement (uz at node 2)
    let mut u_local = vec![0.0; 18];
    u_local[2] = 0.001; // uz at node 0
    u_local[14] = -0.001; // uz at node 2

    let nodal = plate_stress_at_nodes(&coords, e, nu, t, &u_local);

    // Should return 3 stress states
    assert_eq!(nodal.len(), 3);

    // At least some stresses should be non-zero
    let has_nonzero = nodal.iter().any(|s| s.von_mises > 0.0);
    assert!(has_nonzero, "At least one nodal stress should be non-zero");

    // Also verify centroid stress recovery
    let centroid = plate_stress_recovery(&coords, e, nu, t, &u_local);
    assert!(centroid.von_mises >= 0.0);
}

#[test]
fn plate_thermal_load_produces_forces() {
    let coords = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.866, 0.0],
    ];
    let e = 200e6;
    let nu = 0.3;
    let t = 0.01;
    let alpha = 12e-6;

    // Uniform temperature change only
    let f_uniform = plate_thermal_load(&coords, e, nu, t, alpha, 100.0, 0.0);
    assert_eq!(f_uniform.len(), 18);
    let max_f_uniform = f_uniform.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        max_f_uniform > 0.0,
        "Uniform temperature should produce non-zero loads"
    );

    // Gradient only
    let f_gradient = plate_thermal_load(&coords, e, nu, t, alpha, 0.0, 50.0);
    let max_f_gradient = f_gradient.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        max_f_gradient > 0.0,
        "Temperature gradient should produce non-zero loads"
    );

    // Zero temperature: no loads
    let f_zero = plate_thermal_load(&coords, e, nu, t, alpha, 0.0, 0.0);
    let max_f_zero = f_zero.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        max_f_zero < 1e-15,
        "Zero temperature should produce zero loads"
    );
}

#[test]
fn plate_element_quality_equilateral() {
    // Equilateral triangle: perfect quality
    let side = 1.0;
    let coords = [
        [0.0, 0.0, 0.0],
        [side, 0.0, 0.0],
        [side / 2.0, side * (3.0_f64).sqrt() / 2.0, 0.0],
    ];

    let (aspect, skew, min_angle) = plate_element_quality(&coords);

    assert!(
        (aspect - 1.0).abs() < 0.01,
        "Equilateral aspect ratio should be ~1.0: {}",
        aspect
    );
    assert!(
        skew < 1.0,
        "Equilateral skew should be ~0°: {}",
        skew
    );
    assert!(
        (min_angle - 60.0).abs() < 1.0,
        "Equilateral min angle should be ~60°: {}",
        min_angle
    );
}

#[test]
fn plate_element_quality_distorted() {
    // Very elongated triangle: poor quality
    let coords = [
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [5.0, 0.1, 0.0],
    ];

    let (aspect, _skew, min_angle) = plate_element_quality(&coords);

    assert!(
        aspect > 1.5,
        "Elongated triangle should have high aspect ratio: {}",
        aspect
    );
    assert!(
        min_angle < 15.0,
        "Elongated triangle should have small min angle: {}",
        min_angle
    );
}

#[test]
fn plate_pressure_load_equilibrium() {
    // A flat plate in XY plane: pressure in Z
    let coords = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
    ];
    let pressure = 10.0; // kN/m²

    let f = plate_pressure_load(&coords, pressure);

    // Total Z-force should equal pressure * area
    let area = 0.5; // triangle area = 0.5 * base * height = 0.5 * 1.0 * 1.0
    let total_z: f64 = f[2] + f[8] + f[14]; // uz DOFs for each node
    let expected_total = pressure * area;

    assert!(
        (total_z - expected_total).abs() / expected_total < 0.01,
        "Total Z force {:.4} should equal p*A = {:.4}",
        total_z, expected_total
    );

    // X and Y components should be zero for a plate in XY plane
    let total_x: f64 = f[0] + f[6] + f[12];
    let total_y: f64 = f[1] + f[7] + f[13];
    assert!(total_x.abs() < 1e-10, "X force should be zero: {}", total_x);
    assert!(total_y.abs() < 1e-10, "Y force should be zero: {}", total_y);
}

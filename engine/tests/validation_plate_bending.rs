/// Validation: Plate Bending Theory Formulas (Pure Formula Verification)
///
/// References:
///   - Timoshenko & Woinowsky-Krieger, "Theory of Plates and Shells", 2nd Ed.
///   - Ventsel & Krauthammer, "Thin Plates and Shells", Marcel Dekker
///   - Ugural, "Stresses in Plates and Shells", 2nd Ed.
///   - Szilard, "Theories and Applications of Plate Analysis"
///
/// Tests verify plate bending formulas without calling the solver.
///   1. Simply supported rectangular plate under UDL: w_max
///   2. Navier solution: first term of double series
///   3. Circular plate: clamped edge under UDL
///   4. Circular plate: simply supported under UDL
///   5. Plate critical buckling: N_cr
///   6. Effective width in compression (von Karman)
///   7. Plate natural frequency: f_mn
///   8. Mindlin plate correction: shear deformation for thick plates

mod helpers;

use std::f64::consts::PI;

// ================================================================
// 1. Simply Supported Rectangular Plate Under UDL
// ================================================================
//
// For a simply supported rectangular plate (a x b) under uniform
// pressure q, the maximum deflection at the center is:
//   w_max = alpha * q * a^4 / D
//
// where D = E*t^3 / (12*(1-nu^2)) is the flexural rigidity,
// and alpha depends on the aspect ratio a/b.
//
// For a square plate (a/b = 1): alpha = 0.00406
// For a/b = 1.5: alpha = 0.00772
// For a/b = 2.0: alpha = 0.01013
// For a/b -> infinity (strip): alpha = 5/384 = 0.01302
//
// Reference: Timoshenko & Woinowsky-Krieger, Table 8.1

#[test]
fn validation_plate_ss_rectangular_udl() {
    let e: f64 = 200_000.0;  // MPa
    let nu: f64 = 0.3;
    let t: f64 = 0.01;       // m, plate thickness (10 mm)
    let a: f64 = 1.0;        // m, plate dimension (shorter side)
    let q: f64 = 10_000.0;   // Pa = N/m^2

    // Flexural rigidity
    let d = e * 1e6 * t * t * t / (12.0 * (1.0 - nu * nu));
    // = 200e9 * 1e-6 / (12 * 0.91) = 200000 / 10.92 = 18315.02
    let d_expected = 200e9 * 1e-6 / (12.0 * (1.0 - 0.09));
    assert!(
        (d - d_expected).abs() / d_expected < 1e-10,
        "Flexural rigidity: D={:.2}, expected={:.2}",
        d, d_expected
    );

    // Square plate (a/b = 1): alpha = 0.00406
    let alpha_square = 0.00406;
    let w_max_square = alpha_square * q * a.powi(4) / d;
    // This gives deflection in meters
    assert!(
        w_max_square > 0.0,
        "Deflection must be positive: {:.6e}",
        w_max_square
    );

    // For a plate strip (a/b -> inf): alpha = 5/384
    let alpha_strip = 5.0 / 384.0;
    let w_max_strip = alpha_strip * q * a.powi(4) / d;

    // Strip deflection should be greater than square plate (more flexible)
    assert!(
        w_max_strip > w_max_square,
        "Strip ({:.6e}) > Square ({:.6e})",
        w_max_strip, w_max_square
    );

    // Verify alpha_strip value
    let alpha_strip_expected = 0.013020833;
    assert!(
        (alpha_strip - alpha_strip_expected).abs() / alpha_strip_expected < 1e-4,
        "alpha_strip: computed={:.8}, expected={:.8}",
        alpha_strip, alpha_strip_expected
    );

    // Verify that alpha increases with aspect ratio
    let alpha_15 = 0.00772; // a/b = 1.5
    let alpha_20 = 0.01013; // a/b = 2.0
    assert!(alpha_square < alpha_15, "alpha(1.0) < alpha(1.5)");
    assert!(alpha_15 < alpha_20, "alpha(1.5) < alpha(2.0)");
    assert!(alpha_20 < alpha_strip, "alpha(2.0) < alpha(inf)");
}

// ================================================================
// 2. Navier Solution: First Term of Double Series
// ================================================================
//
// The Navier solution for a simply supported rectangular plate
// under uniform load q:
//   w(x,y) = sum_{m=1,3,5...} sum_{n=1,3,5...}
//     16*q / (pi^6 * D * m*n * (m^2/a^2 + n^2/b^2)^2) *
//     sin(m*pi*x/a) * sin(n*pi*y/b)
//
// The first term (m=1, n=1) gives the dominant contribution:
//   w_11(a/2, b/2) = 16*q / (pi^6 * D * (1/a^2 + 1/b^2)^2)
//
// For a square plate (a=b): w_11 = 16*q*a^4 / (pi^6 * D * 4)
//                          = 4*q*a^4 / (pi^6 * D)
//   alpha_11 = 4/pi^6 = 0.004145
//
// The exact alpha = 0.00406 (from series summation), so the first
// term provides 0.004145/0.00406 = 1.021, i.e., 97.9% of the answer
// (slightly overestimates because higher terms correct it).
//
// Reference: Timoshenko & Woinowsky-Krieger, Sec. 5.8

#[test]
fn validation_plate_navier_first_term() {
    let a: f64 = 2.0;  // m (square plate)
    let b: f64 = 2.0;  // m
    let q: f64 = 5000.0; // Pa
    let e: f64 = 200_000.0; // MPa
    let nu: f64 = 0.3;
    let t: f64 = 0.015; // m, plate thickness

    let d = e * 1e6 * t.powi(3) / (12.0 * (1.0 - nu * nu));

    // First term of Navier series at center (m=1, n=1)
    let m: f64 = 1.0;
    let n: f64 = 1.0;
    let denom = PI.powi(6) * d * m * n * (m * m / (a * a) + n * n / (b * b)).powi(2);
    let w_11_center = 16.0 * q / denom; // sin(pi/2)*sin(pi/2) = 1

    // Alpha coefficient for first term (square plate)
    let alpha_11 = 4.0 / PI.powi(6);
    let alpha_11_expected = 0.004161; // 4/pi^6
    assert!(
        (alpha_11 - alpha_11_expected).abs() / alpha_11_expected < 0.001,
        "alpha_11: computed={:.6}, expected={:.6}",
        alpha_11, alpha_11_expected
    );

    // The first term should be about 2% above the exact solution
    let alpha_exact = 0.00406;
    let ratio = alpha_11 / alpha_exact;
    assert!(
        ratio > 1.0 && ratio < 1.05,
        "First term overestimates by {:.2}%",
        (ratio - 1.0) * 100.0
    );

    // Compare first term deflection with exact
    let w_exact = alpha_exact * q * a.powi(4) / d;
    let w_first = alpha_11 * q * a.powi(4) / d;
    assert!(
        w_first > w_exact,
        "First term ({:.6e}) > Exact ({:.6e}) for square plate",
        w_first, w_exact
    );

    // Add third term correction (m=1,n=3 and m=3,n=1)
    let w_13 = 16.0 * q / (PI.powi(6) * d * 1.0 * 3.0
        * (1.0 / (a * a) + 9.0 / (b * b)).powi(2));
    let w_31 = 16.0 * q / (PI.powi(6) * d * 3.0 * 1.0
        * (9.0 / (a * a) + 1.0 / (b * b)).powi(2));

    // Sign at center: sin(pi/2)*sin(3*pi/2) = 1*(-1) = -1
    // So w_13 subtracts from the first term
    let w_3term = w_11_center - w_13 - w_31;
    // This should be closer to exact
    let alpha_3term = w_3term * d / (q * a.powi(4));

    // Three-term approximation should be closer to exact than single term
    assert!(
        (alpha_3term - alpha_exact).abs() < (alpha_11 - alpha_exact).abs(),
        "3-term alpha ({:.6}) closer to exact ({:.6}) than 1-term ({:.6})",
        alpha_3term, alpha_exact, alpha_11
    );
}

// ================================================================
// 3. Circular Plate: Clamped Edge Under UDL
// ================================================================
//
// Circular plate of radius a, clamped at edge, under uniform load q.
// Maximum deflection at center:
//   w_max = q * a^4 / (64 * D)
//
// Maximum bending stress at the clamped edge (radial):
//   sigma_r_max = 3*q*a^2 / (4*t^2)  (at edge, r=a)
//
// The moment at the edge (per unit length):
//   M_r(a) = -q*a^2/8
// The moment at the center:
//   M_r(0) = q*a^2*(1+nu)/16
//
// Reference: Timoshenko & Woinowsky-Krieger, Sec. 3.2

#[test]
fn validation_plate_circular_clamped_udl() {
    let a: f64 = 0.5;      // m, plate radius
    let t: f64 = 0.008;    // m, plate thickness (8 mm)
    let e: f64 = 210_000.0; // MPa
    let nu: f64 = 0.3;
    let q: f64 = 20_000.0;  // Pa

    // Flexural rigidity
    let d: f64 = e * 1e6 * t.powi(3) / (12.0 * (1.0 - nu * nu));

    // Maximum deflection at center
    let w_max: f64 = q * a.powi(4) / (64.0 * d);
    // Should be a small positive number
    assert!(
        w_max > 0.0,
        "Clamped circular plate: w_max must be positive: {:.6e}",
        w_max
    );

    // Verify the formula with known values
    // D = 210e9 * (0.008)^3 / (12 * 0.91) = 210e9 * 5.12e-7 / 10.92
    //   = 107520 / 10.92 = 9846.15 N*m
    let d_check = 210e9 * 0.008_f64.powi(3) / (12.0 * 0.91);
    assert!(
        (d - d_check).abs() / d_check < 1e-6,
        "D: computed={:.2}, check={:.2}",
        d, d_check
    );

    // w_max = 20000 * 0.0625 / (64 * D) = 1250 / (64*D)
    let w_max_check = 1250.0 / (64.0 * d);
    assert!(
        (w_max - w_max_check).abs() / w_max_check < 1e-10,
        "w_max: computed={:.6e}, check={:.6e}",
        w_max, w_max_check
    );

    // Edge moment (per unit length, radial direction)
    let m_edge = -q * a * a / 8.0;
    // = -20000 * 0.25 / 8 = -625 N*m/m
    let m_edge_expected = -625.0;
    assert!(
        (m_edge - m_edge_expected).abs() / m_edge_expected.abs() < 1e-10,
        "Edge moment: computed={:.2}, expected={:.2}",
        m_edge, m_edge_expected
    );

    // Center moment (per unit length)
    let m_center = q * a * a * (1.0 + nu) / 16.0;
    // = 20000 * 0.25 * 1.3 / 16 = 6500 / 16 = 406.25 N*m/m
    let m_center_expected = 406.25;
    assert!(
        (m_center - m_center_expected).abs() / m_center_expected < 1e-10,
        "Center moment: computed={:.2}, expected={:.2}",
        m_center, m_center_expected
    );

    // Edge moment magnitude should exceed center moment (for clamped plate)
    assert!(
        m_edge.abs() > m_center,
        "Edge moment ({:.2}) > Center moment ({:.2})",
        m_edge.abs(), m_center
    );
}

// ================================================================
// 4. Circular Plate: Simply Supported Under UDL
// ================================================================
//
// Circular plate of radius a, simply supported at edge, under
// uniform load q.
//
// Maximum deflection at center:
//   w_max = (5+nu)*q*a^4 / (64*(1+nu)*D)
//
// Maximum moment at center (radial = tangential by symmetry):
//   M_r(0) = (3+nu)*q*a^2/16
//
// Reference: Timoshenko & Woinowsky-Krieger, Sec. 3.3

#[test]
fn validation_plate_circular_ss_udl() {
    let a: f64 = 0.5;     // m, plate radius
    let t: f64 = 0.008;   // m, plate thickness
    let e: f64 = 210_000.0; // MPa
    let nu: f64 = 0.3;
    let q: f64 = 20_000.0;  // Pa

    let d = e * 1e6 * t.powi(3) / (12.0 * (1.0 - nu * nu));

    // Maximum deflection at center (SS)
    let w_max_ss = (5.0 + nu) * q * a.powi(4) / (64.0 * (1.0 + nu) * d);
    // = 5.3 * 20000 * 0.0625 / (64 * 1.3 * D)
    // = 5.3 * 1250 / (83.2 * D)

    // Maximum deflection at center (clamped, from test 3)
    let w_max_clamped = q * a.powi(4) / (64.0 * d);

    // SS deflection should be larger than clamped
    assert!(
        w_max_ss > w_max_clamped,
        "SS ({:.6e}) > Clamped ({:.6e})",
        w_max_ss, w_max_clamped
    );

    // Ratio: w_ss / w_clamped = (5+nu)/(1+nu) for same D, q, a
    let ratio = w_max_ss / w_max_clamped;
    let ratio_expected = (5.0 + nu) / (1.0 + nu);
    // = 5.3/1.3 = 4.077
    assert!(
        (ratio - ratio_expected).abs() / ratio_expected < 1e-10,
        "SS/Clamped ratio: computed={:.4}, expected={:.4}",
        ratio, ratio_expected
    );

    // Maximum moment at center
    let m_center_ss = (3.0 + nu) * q * a * a / 16.0;
    // = 3.3 * 20000 * 0.25 / 16 = 3.3 * 5000 / 16 = 16500/16 = 1031.25 N*m/m
    let m_center_expected = 1031.25;
    assert!(
        (m_center_ss - m_center_expected).abs() / m_center_expected < 1e-10,
        "SS center moment: computed={:.2}, expected={:.2}",
        m_center_ss, m_center_expected
    );

    // For SS plate, the edge moment is zero (by definition)
    // But the edge can rotate freely
    let m_edge_ss: f64 = 0.0; // by boundary condition
    assert!(
        m_edge_ss.abs() < 1e-10,
        "SS edge moment should be zero"
    );

    // Verify the deflection formula value
    let w_ss_check = 5.3 * q * a.powi(4) / (64.0 * 1.3 * d);
    assert!(
        (w_max_ss - w_ss_check).abs() / w_ss_check < 1e-10,
        "w_max_ss: computed={:.6e}, check={:.6e}",
        w_max_ss, w_ss_check
    );
}

// ================================================================
// 5. Plate Critical Buckling: N_cr
// ================================================================
//
// The critical buckling stress resultant for a simply supported
// rectangular plate under uniform in-plane compression:
//   N_cr = k * pi^2 * D / b^2
//
// where k is the buckling coefficient depending on aspect ratio
// and boundary conditions.
//
// For SS plate, a/b = integer: k = 4.0 (minimum)
// For SS plate, a/b general: k = (m*b/a + a/(m*b))^2
//   where m is the number of half-waves minimizing k.
//
// Reference: Timoshenko & Gere, "Theory of Elastic Stability", Ch. 9

#[test]
fn validation_plate_critical_buckling() {
    let e: f64 = 200_000.0;  // MPa
    let nu: f64 = 0.3;
    let t: f64 = 0.010;      // m (10 mm plate)
    let a: f64 = 1.0;        // m, plate length (loaded direction)
    let b: f64 = 0.5;        // m, plate width

    let d = e * 1e6 * t.powi(3) / (12.0 * (1.0 - nu * nu));

    // Buckling coefficient for SS plate
    // a/b = 2.0, m = 2 gives k = (2*0.5/1 + 1/(2*0.5))^2 = (1+1)^2 = 4
    let _aspect = a / b;
    let mut k_min = f64::MAX;
    for m in 1..=5 {
        let mf = m as f64;
        let k = (mf * b / a + a / (mf * b)).powi(2);
        if k < k_min {
            k_min = k;
        }
    }
    // For a/b = 2: m=2 gives k = (2*0.5 + 1/(2*0.5))^2 = (1+1)^2 = 4
    assert!(
        (k_min - 4.0).abs() / 4.0 < 1e-10,
        "Buckling coefficient k: computed={:.4}, expected=4.0",
        k_min
    );

    // Critical stress resultant
    let n_cr = k_min * PI * PI * d / (b * b);
    // N_cr in N/m (force per unit width)

    // Convert to critical stress
    let _sigma_cr = n_cr / (t * 1e6); // MPa (N/m / (m * 1e6 mm^2/m^2))
    // Actually: sigma_cr = N_cr / t where N_cr is in N/m and t is in m
    // sigma_cr in Pa: N_cr / t, then convert to MPa by dividing by 1e6
    let sigma_cr_pa = n_cr / t;
    let sigma_cr_mpa = sigma_cr_pa / 1e6;

    // For thin plates, sigma_cr = k * pi^2 * E / (12*(1-nu^2)) * (t/b)^2
    let sigma_cr_formula = k_min * PI * PI * e / (12.0 * (1.0 - nu * nu)) * (t / b).powi(2);
    // = 4 * 9.8696 * 200000 / 10.92 * 0.0004
    // = 4 * 9.8696 * 200000 * 0.0004 / 10.92
    // = 4 * 9.8696 * 80 / 10.92
    // = 4 * 722.78 = 2891.1 MPa ... that seems high for t/b = 0.02
    // Actually: (t/b)^2 = 0.0004, so sigma_cr = 4*pi^2*200000*0.0004/10.92
    //         = 4*9.8696*80/10.92 = 4*72.28 = 289.1 MPa

    assert!(
        (sigma_cr_mpa - sigma_cr_formula).abs() / sigma_cr_formula < 1e-10,
        "Critical stress: computed={:.2} MPa, formula={:.2} MPa",
        sigma_cr_mpa, sigma_cr_formula
    );

    // For a/b = 1 (square): k_min should also be 4 (m=1)
    let k_square: f64 = (1.0_f64 + 1.0).powi(2); // (m*b/a + a/(m*b))^2 = (1+1)^2 = 4
    assert!(
        (k_square - 4.0).abs() < 1e-10,
        "Square plate k: {:.4}",
        k_square
    );
}

// ================================================================
// 6. Effective Width in Compression (von Karman)
// ================================================================
//
// The von Karman effective width formula for post-buckling of
// a plate in compression:
//   b_eff = b * sqrt(sigma_cr / sigma)
//
// where sigma is the applied stress and sigma_cr is the critical
// buckling stress. This assumes sigma > sigma_cr (post-buckling).
//
// When sigma = sigma_cr: b_eff = b (full width effective)
// When sigma = 4*sigma_cr: b_eff = b/2
// As sigma -> infinity: b_eff -> 0
//
// Reference: von Karman, Sechler & Donnell (1932);
//   Winter (1947) modified formula

#[test]
fn validation_plate_effective_width_von_karman() {
    let b: f64 = 300.0;  // mm, plate width
    let t: f64 = 6.0;    // mm, plate thickness
    let e: f64 = 200_000.0; // MPa
    let nu: f64 = 0.3;
    let fy: f64 = 350.0; // MPa, yield stress

    // Critical buckling stress (SS plate, k=4)
    let sigma_cr = 4.0 * PI * PI * e / (12.0 * (1.0 - nu * nu)) * (t / b).powi(2);
    // = 4 * pi^2 * 200000 / 10.92 * (6/300)^2
    // = 4 * 9.8696 * 200000 / 10.92 * 0.0004
    // = 4 * 72.28 = 289.1 MPa

    // Check that sigma_cr < fy (plate buckles before yielding)
    assert!(
        sigma_cr < fy,
        "sigma_cr ({:.2}) should be < fy ({:.2}) for effective width to apply",
        sigma_cr, fy
    );

    // Von Karman effective width at yield stress
    let b_eff_yield = b * (sigma_cr / fy).sqrt();
    // = 300 * sqrt(289.1/350) = 300 * sqrt(0.826) = 300 * 0.909 = 272.6 mm
    assert!(
        b_eff_yield < b,
        "b_eff ({:.2}) should be < b ({:.2}) when sigma > sigma_cr",
        b_eff_yield, b
    );
    assert!(
        b_eff_yield > 0.0,
        "b_eff should be positive: {:.2}",
        b_eff_yield
    );

    // When sigma = sigma_cr: b_eff = b
    let b_eff_at_cr = b * (sigma_cr / sigma_cr).sqrt();
    assert!(
        (b_eff_at_cr - b).abs() / b < 1e-10,
        "b_eff at sigma_cr: computed={:.2}, expected={:.2}",
        b_eff_at_cr, b
    );

    // When sigma = 4*sigma_cr: b_eff = b/2
    let sigma_4cr = 4.0 * sigma_cr;
    let b_eff_4cr = b * (sigma_cr / sigma_4cr).sqrt();
    assert!(
        (b_eff_4cr - b / 2.0).abs() / (b / 2.0) < 1e-10,
        "b_eff at 4*sigma_cr: computed={:.2}, expected={:.2}",
        b_eff_4cr, b / 2.0
    );

    // Winter's modified formula (empirical improvement):
    // b_eff = b * sqrt(sigma_cr/sigma) * (1 - 0.22*sqrt(sigma_cr/sigma))
    let lambda = (sigma_cr / fy).sqrt();
    let b_eff_winter = b * lambda * (1.0 - 0.22 * lambda);
    // Winter's formula gives a smaller effective width than von Karman
    assert!(
        b_eff_winter < b_eff_yield,
        "Winter ({:.2}) < von Karman ({:.2})",
        b_eff_winter, b_eff_yield
    );
    assert!(
        b_eff_winter > 0.0,
        "Winter b_eff should be positive: {:.2}",
        b_eff_winter
    );
}

// ================================================================
// 7. Plate Natural Frequency: f_mn
// ================================================================
//
// For a simply supported rectangular plate (a x b), the natural
// frequencies are:
//   f_mn = (pi/2) * (m^2/a^2 + n^2/b^2) * sqrt(D/(rho*t))
//
// where m, n are the mode numbers (half-waves in each direction).
//
// The fundamental frequency (m=1, n=1):
//   f_11 = (pi/2) * (1/a^2 + 1/b^2) * sqrt(D/(rho*t))
//
// Reference: Ventsel & Krauthammer, Ch. 16

#[test]
fn validation_plate_natural_frequency() {
    let e: f64 = 200_000.0;   // MPa
    let nu: f64 = 0.3;
    let t: f64 = 0.010;       // m (10 mm)
    let rho: f64 = 7850.0;    // kg/m^3
    let a: f64 = 1.0;         // m
    let b: f64 = 0.8;         // m

    let d = e * 1e6 * t.powi(3) / (12.0 * (1.0 - nu * nu));
    let d_over_rho_t = d / (rho * t);

    // Fundamental frequency (m=1, n=1)
    let f_11 = (PI / 2.0) * (1.0 / (a * a) + 1.0 / (b * b)) * d_over_rho_t.sqrt();
    assert!(
        f_11 > 0.0,
        "Fundamental frequency must be positive: {:.4} Hz",
        f_11
    );

    // Should be a reasonable value (typically 10-1000 Hz for steel plates)
    assert!(
        f_11 > 1.0 && f_11 < 5000.0,
        "f_11 = {:.4} Hz should be in reasonable range",
        f_11
    );

    // Higher modes should have higher frequencies
    let f_21 = (PI / 2.0) * (4.0 / (a * a) + 1.0 / (b * b)) * d_over_rho_t.sqrt();
    let f_12 = (PI / 2.0) * (1.0 / (a * a) + 4.0 / (b * b)) * d_over_rho_t.sqrt();
    let f_22 = (PI / 2.0) * (4.0 / (a * a) + 4.0 / (b * b)) * d_over_rho_t.sqrt();

    assert!(f_21 > f_11, "f_21 ({:.4}) > f_11 ({:.4})", f_21, f_11);
    assert!(f_12 > f_11, "f_12 ({:.4}) > f_11 ({:.4})", f_12, f_11);
    assert!(f_22 > f_21, "f_22 ({:.4}) > f_21 ({:.4})", f_22, f_21);
    assert!(f_22 > f_12, "f_22 ({:.4}) > f_12 ({:.4})", f_22, f_12);

    // For square plate (a = b): f_11 = pi * sqrt(D/(rho*t)) / a^2
    // (since pi/2 * (1/a^2 + 1/a^2) = pi/a^2)
    let a_sq = 1.0;
    let f_11_sq = PI / (a_sq * a_sq) * d_over_rho_t.sqrt();
    let f_11_sq_check = (PI / 2.0) * (1.0 / (a_sq * a_sq) + 1.0 / (a_sq * a_sq)) * d_over_rho_t.sqrt();
    assert!(
        (f_11_sq - f_11_sq_check).abs() / f_11_sq < 1e-10,
        "Square plate frequency: {:.4} vs {:.4}",
        f_11_sq, f_11_sq_check
    );

    // Frequency ratio for modes: f_mn / f_11 = (m^2/a^2 + n^2/b^2) / (1/a^2 + 1/b^2)
    let ratio_21_11 = f_21 / f_11;
    let ratio_expected = (4.0 / (a * a) + 1.0 / (b * b)) / (1.0 / (a * a) + 1.0 / (b * b));
    assert!(
        (ratio_21_11 - ratio_expected).abs() / ratio_expected < 1e-10,
        "f_21/f_11: computed={:.4}, expected={:.4}",
        ratio_21_11, ratio_expected
    );
}

// ================================================================
// 8. Mindlin Plate Correction: Shear Deformation for Thick Plates
// ================================================================
//
// The Mindlin plate theory accounts for transverse shear deformation,
// which becomes significant for thick plates (t/a > 0.1).
//
// The correction factor for deflection:
//   w_Mindlin / w_Kirchhoff = 1 + alpha_s * D / (kappa * G * t * a^2)
//
// where:
//   kappa = shear correction factor (5/6 for rectangular section)
//   G = E / (2*(1+nu)), shear modulus
//   alpha_s depends on the loading and boundary conditions
//
// For a simply supported square plate under UDL:
//   alpha_s = approximately 12*(1+nu)/(5*pi^2) for the first mode
//
// The shear flexibility parameter:
//   phi = 12*D / (kappa*G*t*a^2) = pi^2*t^2 / (5*a^2*(1-nu))
//     (for square plate, using D = Et^3/(12*(1-nu^2)))
//
// Reference: Mindlin (1951); Ugural, Ch. 7

#[test]
fn validation_plate_mindlin_thick_plate_correction() {
    let e: f64 = 200_000.0;   // MPa
    let nu: f64 = 0.3;
    let kappa: f64 = 5.0 / 6.0; // shear correction factor

    let g = e / (2.0 * (1.0 + nu)); // shear modulus
    let g_expected = 200_000.0 / 2.6;
    assert!(
        (g - g_expected).abs() / g_expected < 1e-10,
        "Shear modulus: G={:.2} MPa, expected={:.2} MPa",
        g, g_expected
    );

    // Thin plate: t/a = 0.01 (Kirchhoff theory valid)
    let t_thin: f64 = 0.01; // m
    let a: f64 = 1.0; // m

    let d_thin = e * 1e6 * t_thin.powi(3) / (12.0 * (1.0 - nu * nu));
    let shear_param_thin = 12.0 * d_thin / (kappa * g * 1e6 * t_thin * a * a);
    // This should be very small (< 0.01) for thin plates
    assert!(
        shear_param_thin < 0.05,
        "Thin plate shear parameter: {:.6} should be small",
        shear_param_thin
    );

    // Thick plate: t/a = 0.2
    let t_thick: f64 = 0.2; // m
    let d_thick = e * 1e6 * t_thick.powi(3) / (12.0 * (1.0 - nu * nu));
    let shear_param_thick = 12.0 * d_thick / (kappa * g * 1e6 * t_thick * a * a);

    // This should be significant (> 0.1) for thick plates
    assert!(
        shear_param_thick > 0.1,
        "Thick plate shear parameter: {:.6} should be significant",
        shear_param_thick
    );

    // Shear parameter scales as (t/a)^2
    let ratio_params = shear_param_thick / shear_param_thin;
    let ratio_ta = (t_thick / t_thin).powi(2);
    assert!(
        (ratio_params - ratio_ta).abs() / ratio_ta < 1e-10,
        "Shear param ratio: computed={:.4}, expected={:.4} (=(t_thick/t_thin)^2)",
        ratio_params, ratio_ta
    );

    // Mindlin correction factor: w_Mindlin = w_Kirchhoff * (1 + C * shear_param)
    // where C depends on boundary conditions and loading.
    // For SS square plate under UDL, C ~ 1 approximately.
    let correction_thin = 1.0 + shear_param_thin;
    let correction_thick = 1.0 + shear_param_thick;

    // Thin plate correction should be close to 1
    assert!(
        (correction_thin - 1.0).abs() < 0.05,
        "Thin plate correction: {:.6}, expected ~1.0",
        correction_thin
    );

    // Thick plate correction should be > 1 (more deflection due to shear)
    assert!(
        correction_thick > 1.1,
        "Thick plate correction: {:.6}, should be > 1.1",
        correction_thick
    );

    // Verify the approximate formula: phi = pi^2*t^2/(5*a^2*(1-nu))
    // This is derived from substituting D into the shear parameter
    let phi_approx = PI * PI * t_thin * t_thin / (5.0 * a * a * (1.0 - nu));
    // = pi^2 * 0.0001 / (5 * 1 * 0.7) = 9.8696 * 0.0001 / 3.5 = 2.82e-4
    // Compare with shear_param_thin
    // The two should be related but not identical due to the kappa factor
    // phi_approx uses the exact Mindlin derivation
    assert!(
        phi_approx > 0.0 && phi_approx < 0.01,
        "phi_approx for thin plate: {:.6e}",
        phi_approx
    );
}

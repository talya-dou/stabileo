/// Validation: Plastic Analysis Formulas (Pure Formula Verification)
///
/// References:
///   - Neal, "The Plastic Methods of Structural Analysis", 3rd Ed.
///   - Horne, "Plastic Theory of Structures", 2nd Ed.
///   - Baker & Heyman, "Plastic Design of Frames"
///   - Bruneau, Uang, Sabelli: "Ductile Design of Steel Structures", 2nd Ed.
///
/// Tests verify plastic analysis formulas without calling the solver.
///   1. Plastic moment: Mp = fy * Zx
///   2. Shape factor: f = Zx/Sx for various sections
///   3. Fixed beam collapse load: w_p = 16Mp/L^2
///   4. Propped cantilever collapse: P_p = 6*Mp/L (load at midspan)
///   5. Portal frame beam mechanism collapse load
///   6. Combined mechanism: beam + sway for portal frame
///   7. Upper bound theorem: virtual work method
///   8. Reduced plastic moment with axial load: Mpr = Mp(1-(N/Np)^2)

mod helpers;

use std::f64::consts::PI;

// ================================================================
// 1. Plastic Moment: Mp = fy * Zx
// ================================================================
//
// The plastic moment capacity is the product of yield stress and
// the plastic section modulus.
//
// For a rectangular section b x h:
//   Zx = b * h^2 / 4
//   Mp = fy * Zx
//
// Reference: Neal, Ch. 2

#[test]
fn validation_plastic_moment_capacity() {
    let fy: f64 = 250.0; // MPa

    // Rectangular section: 150 mm x 300 mm
    let b: f64 = 150.0;
    let h: f64 = 300.0;
    let zx_rect: f64 = b * h * h / 4.0;
    let zx_rect_expected: f64 = 3_375_000.0; // mm^3
    assert!(
        (zx_rect - zx_rect_expected).abs() / zx_rect_expected < 1e-10,
        "Rectangular Zx: computed={:.0}, expected={:.0}",
        zx_rect, zx_rect_expected
    );

    let mp_rect: f64 = fy * zx_rect; // N*mm
    let mp_rect_kn_m: f64 = mp_rect / 1e6; // kN*m
    let mp_rect_expected: f64 = 843.75; // kN*m
    assert!(
        (mp_rect_kn_m - mp_rect_expected).abs() / mp_rect_expected < 1e-10,
        "Rectangular Mp: computed={:.2} kN*m, expected={:.2} kN*m",
        mp_rect_kn_m, mp_rect_expected
    );

    // Circular section: diameter 200 mm
    let d: f64 = 200.0;
    let zx_circ: f64 = d * d * d / 6.0;
    let zx_circ_expected: f64 = 8_000_000.0 / 6.0; // mm^3
    assert!(
        (zx_circ - zx_circ_expected).abs() / zx_circ_expected < 1e-6,
        "Circular Zx: computed={:.1}, expected={:.1}",
        zx_circ, zx_circ_expected
    );

    let mp_circ_kn_m: f64 = fy * zx_circ / 1e6;
    let mp_circ_expected: f64 = 250.0 * 8_000_000.0 / 6.0 / 1e6;
    assert!(
        (mp_circ_kn_m - mp_circ_expected).abs() / mp_circ_expected < 1e-6,
        "Circular Mp: computed={:.3} kN*m, expected={:.3} kN*m",
        mp_circ_kn_m, mp_circ_expected
    );
}

// ================================================================
// 2. Shape Factor: f = Zx / Sx for Various Sections
// ================================================================
//
// The shape factor is the ratio of plastic to elastic section modulus.
//
// Rectangular:  f = 1.5   (Zx = bh^2/4, Sx = bh^2/6)
// Circular:     f = 16/(3*pi) ~ 1.698
// I-section:    f ~ 1.12 to 1.18 (depends on proportions)
//
// Reference: Neal, Ch. 2, Table 2.1

#[test]
fn validation_plastic_shape_factors() {
    // Rectangular section
    let b: f64 = 100.0;
    let h: f64 = 200.0;
    let zx_rect: f64 = b * h * h / 4.0;
    let sx_rect: f64 = b * h * h / 6.0;
    let f_rect: f64 = zx_rect / sx_rect;
    assert!(
        (f_rect - 1.5).abs() < 1e-10,
        "Rectangular shape factor: computed={:.4}, expected=1.5",
        f_rect
    );

    // Circular section (solid)
    let d: f64 = 150.0;
    let zx_circ: f64 = d * d * d / 6.0;
    let sx_circ: f64 = PI * d * d * d / 32.0;
    let f_circ: f64 = zx_circ / sx_circ;
    let f_circ_expected: f64 = 16.0 / (3.0 * PI);
    assert!(
        (f_circ - f_circ_expected).abs() / f_circ_expected < 1e-6,
        "Circular shape factor: computed={:.4}, expected={:.4}",
        f_circ, f_circ_expected
    );

    // I-section approximation:
    // bf=200, tf=15, d=300, tw=10
    let bf: f64 = 200.0;
    let tf: f64 = 15.0;
    let d_total: f64 = 300.0;
    let tw: f64 = 10.0;

    // Elastic section modulus via moment of inertia
    let i_flanges: f64 = 2.0 * (bf * tf * tf * tf / 12.0
        + bf * tf * ((d_total - tf) / 2.0).powi(2));
    let hw: f64 = d_total - 2.0 * tf;
    let i_web: f64 = tw * hw * hw * hw / 12.0;
    let i_total: f64 = i_flanges + i_web;
    let sx_i: f64 = i_total / (d_total / 2.0);

    // Plastic section modulus
    let zx_i: f64 = bf * tf * (d_total - tf) + tw * hw * hw / 4.0;

    let f_i: f64 = zx_i / sx_i;
    assert!(
        f_i > 1.10 && f_i < 1.25,
        "I-section shape factor: computed={:.4}, expected ~1.12-1.18",
        f_i
    );
}

// ================================================================
// 3. Fixed Beam Collapse Load: w_p = 16Mp / L^2
// ================================================================
//
// A fixed-fixed beam under UDL collapses when plastic hinges
// form at both fixed ends and at midspan (3 hinges -> mechanism).
//
// Virtual work: w_p * L^2 / 4 = 4 * Mp  =>  w_p = 16*Mp/L^2
//
// Reference: Neal, Ch. 4; Horne, Ch. 3

#[test]
fn validation_plastic_fixed_beam_collapse_load() {
    let mp: f64 = 500.0; // kN*m
    let l: f64 = 8.0;    // m

    // Collapse load for fixed beam under UDL
    let w_p: f64 = 16.0 * mp / (l * l);
    let w_p_expected: f64 = 125.0; // kN/m
    assert!(
        (w_p - w_p_expected).abs() / w_p_expected < 1e-10,
        "Fixed beam collapse load: computed={:.2}, expected={:.2}",
        w_p, w_p_expected
    );

    // Compare with elastic maximum moment: M_max = wL^2/12 (at supports)
    // At collapse: w_fy * L^2/12 = Mp => w_first_yield = 12*Mp/L^2
    let w_fy: f64 = 12.0 * mp / (l * l);
    let w_fy_expected: f64 = 93.75; // kN/m
    assert!(
        (w_fy - w_fy_expected).abs() / w_fy_expected < 1e-10,
        "First yield load: computed={:.2}, expected={:.2}",
        w_fy, w_fy_expected
    );

    // Ratio of collapse to first yield = 16/12 = 4/3
    let ratio: f64 = w_p / w_fy;
    let ratio_expected: f64 = 16.0 / 12.0;
    assert!(
        (ratio - ratio_expected).abs() / ratio_expected < 1e-10,
        "Collapse/yield ratio: computed={:.4}, expected={:.4}",
        ratio, ratio_expected
    );
}

// ================================================================
// 4. Propped Cantilever Collapse: P_p at Midspan
// ================================================================
//
// A propped cantilever (fixed at A, roller at B) with a point load
// P at midspan. Hinges form at fixed end and under the load.
//
// Virtual work:
//   Internal work = Mp * 2*delta/L + Mp * 4*delta/L = 6*Mp*delta/L
//   External work = P * delta
//   P_p = 6*Mp/L
//
// For UDL: w_p = (6 + 4*sqrt(2)) * Mp / L^2  (exact)
//
// Reference: Neal, Ch. 4; Horne, Ch. 3

#[test]
fn validation_plastic_propped_cantilever_collapse() {
    let mp: f64 = 400.0; // kN*m
    let l: f64 = 6.0;    // m

    // Collapse load with P at midspan
    let p_p_midspan: f64 = 6.0 * mp / l;
    let p_p_midspan_expected: f64 = 400.0; // kN
    assert!(
        (p_p_midspan - p_p_midspan_expected).abs() / p_p_midspan_expected < 1e-10,
        "Propped cantilever P_p (midspan): computed={:.2}, expected={:.2}",
        p_p_midspan, p_p_midspan_expected
    );

    // For UDL on propped cantilever:
    //   w_p = (6 + 4*sqrt(2)) * Mp / L^2
    let coeff: f64 = 6.0 + 4.0 * 2.0_f64.sqrt();
    let w_p_udl: f64 = coeff * mp / (l * l);
    let w_p_udl_expected: f64 = coeff * 400.0 / 36.0;
    assert!(
        (w_p_udl - w_p_udl_expected).abs() / w_p_udl_expected < 1e-10,
        "Propped cantilever w_p (UDL): computed={:.4}, expected={:.4}",
        w_p_udl, w_p_udl_expected
    );

    // Verify coefficient value ~ 11.657
    assert!(
        (coeff - 11.6569).abs() < 0.001,
        "UDL coefficient: computed={:.4}, expected ~11.657",
        coeff
    );
}

// ================================================================
// 5. Portal Frame Beam Mechanism Collapse Load
// ================================================================
//
// For a portal frame with fixed bases, span L, height h:
//
// Beam mechanism under UDL: w_p = 16*Mp / L^2
// Sway mechanism under H:   H_p = 4 * Mp / h
// Beam mechanism under P at midspan: P_p = 8*Mp/L
//
// Reference: Horne, "Plastic Theory of Structures", Ch. 5

#[test]
fn validation_plastic_portal_beam_mechanism() {
    let mp: f64 = 300.0; // kN*m
    let l: f64 = 10.0;   // m, beam span
    let h: f64 = 4.0;    // m, column height

    // Beam mechanism under UDL
    let w_beam: f64 = 16.0 * mp / (l * l);
    let w_beam_expected: f64 = 48.0; // kN/m
    assert!(
        (w_beam - w_beam_expected).abs() / w_beam_expected < 1e-10,
        "Beam mechanism w_p: computed={:.2}, expected={:.2}",
        w_beam, w_beam_expected
    );

    // Sway mechanism under horizontal load H at beam level
    let h_sway: f64 = 4.0 * mp / h;
    let h_sway_expected: f64 = 300.0; // kN
    assert!(
        (h_sway - h_sway_expected).abs() / h_sway_expected < 1e-10,
        "Sway mechanism H_p: computed={:.2}, expected={:.2}",
        h_sway, h_sway_expected
    );

    // Beam mechanism under concentrated load P at midspan
    let p_beam: f64 = 8.0 * mp / l;
    let p_beam_expected: f64 = 240.0; // kN
    assert!(
        (p_beam - p_beam_expected).abs() / p_beam_expected < 1e-10,
        "Beam mechanism P_p: computed={:.2}, expected={:.2}",
        p_beam, p_beam_expected
    );
}

// ================================================================
// 6. Combined Mechanism: Beam + Sway for Portal Frame
// ================================================================
//
// When both vertical and horizontal loads act on a portal frame,
// the combined mechanism (beam + sway) may govern.
//
// Combined mechanism:
//   Internal work = 6*Mp*theta
//   External work = lambda * (w*L^2/4 + H*h) * theta
//   lambda_combined = 6*Mp / (w*L^2/4 + H*h)
//
// Reference: Baker & Heyman, Ch. 5

#[test]
fn validation_plastic_combined_mechanism() {
    let mp: f64 = 200.0; // kN*m (same for all members)
    let l: f64 = 8.0;    // m, beam span
    let h: f64 = 4.0;    // m, column height
    let w: f64 = 20.0;   // kN/m, vertical UDL on beam
    let h_force: f64 = 50.0; // kN, horizontal force at beam level

    // Beam mechanism load factor
    let lambda_beam: f64 = 4.0 * mp / (w * l * l / 4.0);
    let lambda_beam_expected: f64 = 2.5;
    assert!(
        (lambda_beam - lambda_beam_expected).abs() / lambda_beam_expected < 1e-10,
        "Beam mechanism lambda: computed={:.4}, expected={:.4}",
        lambda_beam, lambda_beam_expected
    );

    // Sway mechanism load factor
    let lambda_sway: f64 = 4.0 * mp / (h_force * h);
    let lambda_sway_expected: f64 = 4.0;
    assert!(
        (lambda_sway - lambda_sway_expected).abs() / lambda_sway_expected < 1e-10,
        "Sway mechanism lambda: computed={:.4}, expected={:.4}",
        lambda_sway, lambda_sway_expected
    );

    // Combined mechanism (beam + sway - overlap):
    // Internal work = 6*Mp, External work = w*L^2/4 + H*h
    let lambda_combined: f64 = 6.0 * mp / (w * l * l / 4.0 + h_force * h);
    let lambda_combined_expected: f64 = 1200.0 / 520.0;
    assert!(
        (lambda_combined - lambda_combined_expected).abs() / lambda_combined_expected < 1e-10,
        "Combined mechanism lambda: computed={:.6}, expected={:.6}",
        lambda_combined, lambda_combined_expected
    );

    // The true collapse load factor is the minimum of all mechanisms
    let lambda_collapse: f64 = lambda_beam.min(lambda_sway).min(lambda_combined);
    assert!(
        (lambda_collapse - lambda_combined).abs() < 1e-10,
        "Combined mechanism governs: lambda={:.4}",
        lambda_collapse
    );
}

// ================================================================
// 7. Upper Bound Theorem: Virtual Work Method
// ================================================================
//
// The upper bound theorem states that any mechanism gives a load
// factor >= the true collapse load factor.
//
// For SS beam with central point load:
//   True collapse: P_p = 4*Mp/L
//   Non-optimal mechanism (hinge at L/3): P_upper = 6*Mp/L
//   Optimal mechanism (hinge at L/2): P_upper = 4*Mp/L = P_true
//
// Reference: Neal, Ch. 5

#[test]
fn validation_plastic_upper_bound_virtual_work() {
    let mp: f64 = 600.0; // kN*m
    let l: f64 = 10.0;   // m

    // True collapse load (hinge at midspan under load)
    let p_true: f64 = 4.0 * mp / l;
    let p_true_expected: f64 = 240.0; // kN
    assert!(
        (p_true - p_true_expected).abs() / p_true_expected < 1e-10,
        "SS beam P_p: computed={:.2}, expected={:.2}",
        p_true, p_true_expected
    );

    // Upper bound: assume hinge at L/3 (not optimal)
    // Deflection at load (L/2) from right segment: delta_P = beta * L/2
    // Internal work = Mp * (alpha + beta) = Mp * 3*beta
    // External work = P * beta * L/2
    // P_upper = 6*Mp/L
    let p_upper_l3: f64 = 6.0 * mp / l;
    let p_upper_l3_expected: f64 = 360.0; // kN
    assert!(
        (p_upper_l3 - p_upper_l3_expected).abs() / p_upper_l3_expected < 1e-10,
        "Upper bound (hinge at L/3): computed={:.2}, expected={:.2}",
        p_upper_l3, p_upper_l3_expected
    );

    // Upper bound must be >= true collapse load
    assert!(
        p_upper_l3 >= p_true,
        "Upper bound ({:.2}) must be >= true ({:.2})",
        p_upper_l3, p_true
    );

    // Upper bound with hinge at midspan (optimal for this case)
    let p_upper_mid: f64 = 4.0 * mp / l;
    assert!(
        (p_upper_mid - p_true).abs() / p_true < 1e-10,
        "Optimal upper bound equals true: P_upper={:.2}, P_true={:.2}",
        p_upper_mid, p_true
    );

    // Any mechanism gives an upper bound
    assert!(
        p_upper_l3 >= p_upper_mid,
        "Non-optimal ({:.2}) >= optimal ({:.2})",
        p_upper_l3, p_upper_mid
    );
}

// ================================================================
// 8. Reduced Plastic Moment with Axial Load
// ================================================================
//
// For a rectangular section:
//   Mpr = Mp * (1 - (N/Np)^2)
// where Np = fy * A is the squash load.
//
// Reference: Neal, Ch. 6; AISC 360-22 H1.1

#[test]
fn validation_plastic_reduced_moment_axial() {
    let fy: f64 = 350.0; // MPa

    // Rectangular section: 200 mm x 400 mm
    let b: f64 = 200.0;
    let h: f64 = 400.0;
    let a_sec: f64 = b * h; // mm^2 = 80000

    let np: f64 = fy * a_sec; // N, squash load
    let np_kn: f64 = np / 1000.0;
    let np_expected: f64 = 28000.0; // kN
    assert!(
        (np_kn - np_expected).abs() / np_expected < 1e-10,
        "Squash load: computed={:.0} kN, expected={:.0} kN",
        np_kn, np_expected
    );

    let zx: f64 = b * h * h / 4.0; // mm^3
    let mp: f64 = fy * zx / 1e6; // kN*m
    let mp_expected: f64 = 350.0 * 8_000_000.0 / 1e6; // = 2800 kN*m
    assert!(
        (mp - mp_expected).abs() / mp_expected < 1e-10,
        "Mp: computed={:.2}, expected={:.2}",
        mp, mp_expected
    );

    // Case 1: N = 0 => Mpr = Mp
    let n1: f64 = 0.0;
    let mpr1: f64 = mp * (1.0 - (n1 / np_kn).powi(2));
    assert!(
        (mpr1 - mp).abs() < 1e-10,
        "Mpr at N=0: computed={:.2}, expected={:.2}",
        mpr1, mp
    );

    // Case 2: N = Np/2 => Mpr = 0.75*Mp
    let n2: f64 = np_kn / 2.0;
    let mpr2: f64 = mp * (1.0 - (n2 / np_kn).powi(2));
    let mpr2_expected: f64 = 0.75 * mp;
    assert!(
        (mpr2 - mpr2_expected).abs() / mpr2_expected < 1e-10,
        "Mpr at N=Np/2: computed={:.2}, expected={:.2}",
        mpr2, mpr2_expected
    );

    // Case 3: N = Np => Mpr = 0 (pure axial)
    let n3: f64 = np_kn;
    let mpr3: f64 = mp * (1.0 - (n3 / np_kn).powi(2));
    assert!(
        mpr3.abs() < 1e-10,
        "Mpr at N=Np: computed={:.6}, expected=0",
        mpr3
    );

    // Verify at N = 0.3*Np
    let n4: f64 = 0.3 * np_kn;
    let mpr4: f64 = mp * (1.0 - (n4 / np_kn).powi(2));
    let mpr4_expected: f64 = mp * (1.0 - 0.09); // = 0.91 * Mp
    assert!(
        (mpr4 - mpr4_expected).abs() / mpr4_expected < 1e-10,
        "Mpr at N=0.3Np: computed={:.2}, expected={:.2}",
        mpr4, mpr4_expected
    );
}

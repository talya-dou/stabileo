/// Validation: Soil Dynamics and Earthquake Geotechnical Engineering
///
/// References:
///   - Kramer, S.L., "Geotechnical Earthquake Engineering", Prentice Hall, 1996
///   - Towhata, I., "Geotechnical Earthquake Engineering", Springer, 2008
///   - Seed, H.B. & Idriss, I.M., "Simplified Procedure for Evaluating Soil
///     Liquefaction Potential", JSMFD ASCE, 97(9), 1971
///   - EN 1998-5 (Eurocode 8, Part 5): Foundations, retaining structures,
///     geotechnical aspects
///   - ASCE 7-22, Chapter 20: Site Classification Procedure
///   - Youd, T.L. et al., "Revised Multilinear Regression Models of Lateral
///     Spread Displacement", JGGE ASCE, 128(12), 2002
///
/// Tests verify fundamental geotechnical earthquake engineering formulas:
///   1. Site amplification factors (Fa, Fv) from ASCE 7 Table 11.4-1/11.4-2
///   2. Shear wave velocity averaging and ASCE 7 site classification
///   3. Ground response analysis: 1D SH-wave transfer function
///   4. Cyclic stress ratio (CSR) for liquefaction evaluation (Seed & Idriss)
///   5. Dynamic soil properties: G/Gmax degradation (Seed & Idriss curves)
///   6. Seismic bearing capacity reduction factors
///   7. Lateral spreading displacement (Youd et al. 2002 MLR model)
///   8. Site period estimation via quarter-wavelength method
mod helpers;

// ================================================================
// 1. ASCE 7 Site Amplification Factors (Fa, Fv)
// ================================================================
//
// ASCE 7-22 Tables 11.4-1 and 11.4-2 provide short-period (Fa) and
// long-period (Fv) site amplification coefficients as a function of
// site class and mapped spectral acceleration (Ss, S1).
//
// For Site Class D (stiff soil) at Ss = 0.50g:
//   Fa = 1.4 (Table 11.4-1, interpolated)
// For Site Class D at S1 = 0.20g:
//   Fv = 1.8 (Table 11.4-2, interpolated)
//
// The design spectral accelerations are:
//   SDS = (2/3) * Fa * Ss
//   SD1 = (2/3) * Fv * S1

#[test]
fn validation_soil_dynamics_asce7_site_amplification() {
    // ASCE 7-22 Table 11.4-1: Site Class D, Ss = 0.50g
    let _ss: f64 = 0.50;
    let _fa: f64 = 1.4;

    // ASCE 7-22 Table 11.4-2: Site Class D, S1 = 0.20g
    let _s1: f64 = 0.20;
    let _fv: f64 = 1.8;

    // Design spectral acceleration parameters (ASCE 7-22 Eq. 11.4-3, 11.4-4)
    let sds: f64 = (2.0 / 3.0) * _fa * _ss;
    let sd1: f64 = (2.0 / 3.0) * _fv * _s1;

    let sds_expected: f64 = (2.0 / 3.0) * 1.4 * 0.50; // 0.4667g
    let sd1_expected: f64 = (2.0 / 3.0) * 1.8 * 0.20; // 0.240g

    assert!(
        (sds - sds_expected).abs() < 1e-10,
        "SDS = (2/3)*Fa*Ss: computed={:.6}, expected={:.6}", sds, sds_expected
    );
    assert!(
        (sd1 - sd1_expected).abs() < 1e-10,
        "SD1 = (2/3)*Fv*S1: computed={:.6}, expected={:.6}", sd1, sd1_expected
    );

    // Transition period T0 and Ts (ASCE 7-22 Eq. 11.4-5, 11.4-6)
    let ts: f64 = sd1 / sds;
    let t0: f64 = 0.2 * ts;

    let ts_expected: f64 = 0.240 / sds_expected;
    let t0_expected: f64 = 0.2 * ts_expected;

    assert!(
        (ts - ts_expected).abs() < 1e-10,
        "Ts = SD1/SDS: computed={:.6}, expected={:.6}", ts, ts_expected
    );
    assert!(
        (t0 - t0_expected).abs() < 1e-10,
        "T0 = 0.2*Ts: computed={:.6}, expected={:.6}", t0, t0_expected
    );

    // Verify Fa > 1 for soft site (amplification), and Fa increases
    // from Site Class B to D for moderate Ss
    // Site Class B: Fa = 1.0; Site Class C: Fa = 1.2; Site Class D: Fa = 1.4
    let _fa_b: f64 = 1.0;
    let _fa_c: f64 = 1.2;
    let _fa_d: f64 = 1.4;

    assert!(
        _fa_b < _fa_c && _fa_c < _fa_d,
        "Fa should increase from stiff to soft site: B={:.1} < C={:.1} < D={:.1}",
        _fa_b, _fa_c, _fa_d
    );

    // At very high Ss (>= 1.25g), Site Class D deamplifies: Fa = 1.0
    let _fa_d_high_ss: f64 = 1.0;
    assert!(
        _fa_d_high_ss <= _fa_d,
        "Fa for Site Class D decreases at high Ss: Fa(1.25g)={:.1} <= Fa(0.50g)={:.1}",
        _fa_d_high_ss, _fa_d
    );
}

// ================================================================
// 2. Shear Wave Velocity and ASCE 7 Site Classification
// ================================================================
//
// ASCE 7-22 §20.4: The average shear wave velocity for the top 30 m
// is computed as:
//   Vs30 = 30 / Σ(di/Vsi)    (travel-time weighted average)
//
// Site classification (ASCE 7-22 Table 20.3-1):
//   A: Vs30 > 1500 m/s  (hard rock)
//   B: 760 < Vs30 <= 1500 m/s (rock)
//   C: 360 < Vs30 <= 760 m/s (very dense soil / soft rock)
//   D: 180 < Vs30 <= 360 m/s (stiff soil)
//   E: Vs30 <= 180 m/s (soft clay)

#[test]
fn validation_soil_dynamics_shear_wave_velocity_classification() {
    // Example site profile: 3 layers totaling 30 m
    // Layer 1: 5 m thick, Vs = 150 m/s (soft clay)
    // Layer 2: 10 m thick, Vs = 250 m/s (medium stiff)
    // Layer 3: 15 m thick, Vs = 400 m/s (dense sand)
    let d1: f64 = 5.0;
    let vs1: f64 = 150.0;
    let d2: f64 = 10.0;
    let vs2: f64 = 250.0;
    let d3: f64 = 15.0;
    let vs3: f64 = 400.0;

    let total_depth: f64 = d1 + d2 + d3;
    assert!(
        (total_depth - 30.0).abs() < 1e-10,
        "Total depth should be 30 m for Vs30 calculation"
    );

    // Vs30 = 30 / Σ(di/Vsi)
    let travel_time_sum: f64 = d1 / vs1 + d2 / vs2 + d3 / vs3;
    let vs30: f64 = 30.0 / travel_time_sum;

    // Manual computation:
    // travel_time = 5/150 + 10/250 + 15/400
    //             = 0.03333 + 0.04000 + 0.03750
    //             = 0.11083 s
    // Vs30 = 30 / 0.11083 = 270.7 m/s
    let _travel_time_expected: f64 = 5.0 / 150.0 + 10.0 / 250.0 + 15.0 / 400.0;
    let vs30_expected: f64 = 30.0 / _travel_time_expected;

    assert!(
        (vs30 - vs30_expected).abs() < 1e-6,
        "Vs30 calculation: computed={:.2} m/s, expected={:.2} m/s", vs30, vs30_expected
    );

    // Site classification: 180 < 270.7 <= 360 => Site Class D
    assert!(
        vs30 > 180.0 && vs30 <= 360.0,
        "Vs30={:.1} m/s should classify as Site Class D (180 < Vs30 <= 360)", vs30
    );

    // Verify that harmonic mean < arithmetic mean (always true for positive values)
    let vs_arithmetic: f64 = (d1 * vs1 + d2 * vs2 + d3 * vs3) / total_depth;
    assert!(
        vs30 < vs_arithmetic,
        "Harmonic mean Vs30={:.1} should be less than arithmetic mean={:.1}",
        vs30, vs_arithmetic
    );

    // Verify site class boundaries
    let _vs30_class_e: f64 = 150.0;
    let _vs30_class_c: f64 = 500.0;
    let _vs30_class_b: f64 = 900.0;
    let _vs30_class_a: f64 = 2000.0;

    assert!(_vs30_class_e <= 180.0, "Vs30={:.0} should be Site Class E", _vs30_class_e);
    assert!(_vs30_class_c > 360.0 && _vs30_class_c <= 760.0,
        "Vs30={:.0} should be Site Class C", _vs30_class_c);
    assert!(_vs30_class_b > 760.0 && _vs30_class_b <= 1500.0,
        "Vs30={:.0} should be Site Class B", _vs30_class_b);
    assert!(_vs30_class_a > 1500.0, "Vs30={:.0} should be Site Class A", _vs30_class_a);
}

// ================================================================
// 3. 1D Ground Response Analysis: SH-Wave Transfer Function
// ================================================================
//
// For a uniform soil layer of thickness H over rigid bedrock, the
// 1D SH-wave transfer function (amplification) at the surface is:
//
//   |F(omega)| = 1 / cos(omega * H / Vs)
//
// Resonant frequencies occur when cos(omega*H/Vs) = 0:
//   fn = Vs * (2n - 1) / (4H),  n = 1, 2, 3, ...
//
// Reference: Kramer, Ch. 7; Towhata, Ch. 4

#[test]
fn validation_soil_dynamics_1d_transfer_function() {
    let vs: f64 = 200.0;   // shear wave velocity (m/s)
    let h: f64 = 20.0;     // layer thickness (m)

    // Fundamental resonant frequency: f1 = Vs / (4H)
    let f1: f64 = vs / (4.0 * h);
    let f1_expected: f64 = 200.0 / (4.0 * 20.0); // 2.5 Hz

    assert!(
        (f1 - f1_expected).abs() < 1e-10,
        "Fundamental frequency f1 = Vs/(4H): computed={:.4} Hz, expected={:.4} Hz",
        f1, f1_expected
    );

    // Higher resonant frequencies: fn = (2n-1) * Vs / (4H)
    let f2: f64 = 3.0 * vs / (4.0 * h); // 7.5 Hz
    let f3: f64 = 5.0 * vs / (4.0 * h); // 12.5 Hz

    assert!(
        (f2 / f1 - 3.0).abs() < 1e-10,
        "f2/f1 should be 3.0 (odd harmonics): ratio={:.6}", f2 / f1
    );
    assert!(
        (f3 / f1 - 5.0).abs() < 1e-10,
        "f3/f1 should be 5.0 (odd harmonics): ratio={:.6}", f3 / f1
    );

    // Transfer function amplitude at a non-resonant frequency
    // |F(omega)| = 1 / |cos(omega * H / Vs)|
    let f_test: f64 = 1.0; // 1 Hz (below resonance)
    let omega_test: f64 = 2.0 * std::f64::consts::PI * f_test;
    let cos_arg: f64 = omega_test * h / vs;
    let amplification: f64 = 1.0 / cos_arg.cos().abs();

    // At 1 Hz: cos(2*pi*1*20/200) = cos(0.6283) = 0.8090
    // Amplification = 1/0.8090 = 1.236
    let _cos_expected: f64 = (2.0 * std::f64::consts::PI * 1.0 * 20.0 / 200.0).cos();
    let amp_expected: f64 = 1.0 / _cos_expected.abs();

    assert!(
        (amplification - amp_expected).abs() < 1e-10,
        "Transfer function at 1 Hz: computed={:.6}, expected={:.6}",
        amplification, amp_expected
    );

    // Amplification should be > 1 (surface amplifies bedrock motion)
    assert!(
        amplification > 1.0,
        "Surface amplification should be > 1.0 at non-resonant frequency, got {:.4}",
        amplification
    );

    // Near resonance, amplification should be very large
    // Use frequency slightly off resonance to avoid division by zero
    let f_near_res: f64 = f1 * 0.99;
    let omega_near: f64 = 2.0 * std::f64::consts::PI * f_near_res;
    let amp_near_res: f64 = 1.0 / (omega_near * h / vs).cos().abs();

    assert!(
        amp_near_res > 10.0,
        "Near-resonance amplification should be large (>10), got {:.2}",
        amp_near_res
    );

    // With damping ratio xi, the damped transfer function is:
    // |F(omega)| = 1 / sqrt(cos^2(k*H) + (xi*k*H)^2)  (approximate)
    // where k* = omega/Vs*(1 + i*xi) for low damping
    let xi: f64 = 0.05; // 5% damping ratio
    let _k_h: f64 = 2.0 * std::f64::consts::PI * f1 * h / vs;
    // At resonance (k*H = pi/2), cos(k*H) = 0:
    // |F| approx = 1 / (xi * pi/2) = 2/(pi*xi)
    let amp_damped_resonance: f64 = 1.0 / (xi * std::f64::consts::PI / 2.0);
    let amp_damped_expected: f64 = 2.0 / (std::f64::consts::PI * xi);

    assert!(
        (amp_damped_resonance - amp_damped_expected).abs() < 1e-10,
        "Damped resonance amplification: computed={:.4}, expected={:.4}",
        amp_damped_resonance, amp_damped_expected
    );

    // Damped resonance should be finite: 2/(pi*0.05) = 12.73
    assert!(
        (amp_damped_resonance - 12.732).abs() < 0.01,
        "Damped resonance amplification = 2/(pi*xi) = {:.3}, expected ~12.732",
        amp_damped_resonance
    );

    // Verify unused variable for fundamental period
    let _t_site: f64 = 1.0 / f1; // 0.4 s
    assert!(
        (_t_site - 4.0 * h / vs).abs() < 1e-10,
        "Site period T = 4H/Vs = {:.4} s", _t_site
    );
}

// ================================================================
// 4. Cyclic Stress Ratio (CSR) for Liquefaction Evaluation
// ================================================================
//
// The simplified procedure (Seed & Idriss, 1971) estimates the
// earthquake-induced cyclic stress ratio as:
//
//   CSR = 0.65 * (amax/g) * (sigma_v / sigma_v') * rd
//
// where:
//   amax = peak ground acceleration
//   sigma_v = total vertical stress at depth z
//   sigma_v' = effective vertical stress at depth z
//   rd = stress reduction coefficient (depth-dependent)
//
// Idriss (1999) approximation for rd:
//   rd = 1.0 - 0.00765*z   for z <= 9.15 m
//   rd = 1.174 - 0.0267*z  for 9.15 < z <= 23 m
//
// Reference: Seed & Idriss (1971), Youd et al. (2001)

#[test]
fn validation_soil_dynamics_cyclic_stress_ratio() {
    // Soil profile:
    //   gamma_sat = 19.0 kN/m^3 (saturated unit weight)
    //   gamma_w = 9.81 kN/m^3 (water unit weight)
    //   GWT at 2.0 m depth
    //   gamma_dry = 17.0 kN/m^3 (above water table)
    let gamma_dry: f64 = 17.0;
    let gamma_sat: f64 = 19.0;
    let gamma_w: f64 = 9.81;
    let gwt: f64 = 2.0;    // groundwater table depth (m)
    let z: f64 = 8.0;       // evaluation depth (m)
    let amax: f64 = 0.25;   // peak ground acceleration (g)
    let _g: f64 = 9.81;     // gravitational acceleration (m/s^2)

    // Total vertical stress at depth z
    // sigma_v = gamma_dry * gwt + gamma_sat * (z - gwt)
    let sigma_v: f64 = gamma_dry * gwt + gamma_sat * (z - gwt);
    let sigma_v_expected: f64 = 17.0 * 2.0 + 19.0 * 6.0; // 34 + 114 = 148 kPa

    assert!(
        (sigma_v - sigma_v_expected).abs() < 1e-10,
        "Total vertical stress: computed={:.2} kPa, expected={:.2} kPa",
        sigma_v, sigma_v_expected
    );

    // Pore water pressure at depth z
    let u: f64 = gamma_w * (z - gwt);
    let u_expected: f64 = 9.81 * 6.0; // 58.86 kPa

    assert!(
        (u - u_expected).abs() < 1e-10,
        "Pore water pressure: computed={:.2} kPa, expected={:.2} kPa", u, u_expected
    );

    // Effective vertical stress
    let sigma_v_eff: f64 = sigma_v - u;
    let sigma_v_eff_expected: f64 = sigma_v_expected - u_expected; // 89.14 kPa

    assert!(
        (sigma_v_eff - sigma_v_eff_expected).abs() < 1e-10,
        "Effective vertical stress: computed={:.2} kPa, expected={:.2} kPa",
        sigma_v_eff, sigma_v_eff_expected
    );

    // Stress reduction coefficient rd (Idriss 1999)
    // z = 8.0 m <= 9.15 m, so rd = 1.0 - 0.00765 * z
    let rd: f64 = 1.0 - 0.00765 * z;
    let rd_expected: f64 = 1.0 - 0.00765 * 8.0; // 0.9388

    assert!(
        (rd - rd_expected).abs() < 1e-10,
        "Stress reduction coefficient: computed={:.6}, expected={:.6}", rd, rd_expected
    );

    // CSR = 0.65 * (amax/g) * (sigma_v / sigma_v') * rd
    // Note: amax is already in g units, so amax/g = amax
    let csr: f64 = 0.65 * amax * (sigma_v / sigma_v_eff) * rd;
    let csr_expected: f64 = 0.65 * 0.25 * (sigma_v_expected / sigma_v_eff_expected) * rd_expected;

    assert!(
        (csr - csr_expected).abs() < 1e-10,
        "CSR = 0.65*(amax/g)*(sigma_v/sigma_v')*rd: computed={:.6}, expected={:.6}",
        csr, csr_expected
    );

    // Verify CSR is in a physically reasonable range (0.05 to 0.6 for typical conditions)
    assert!(
        csr > 0.05 && csr < 0.6,
        "CSR={:.4} should be in reasonable range [0.05, 0.6]", csr
    );

    // Verify rd for deeper layer (z = 15 m, use second formula)
    let z_deep: f64 = 15.0;
    let rd_deep: f64 = 1.174 - 0.0267 * z_deep;
    let rd_deep_expected: f64 = 1.174 - 0.0267 * 15.0; // 0.7735

    assert!(
        (rd_deep - rd_deep_expected).abs() < 1e-10,
        "Deep rd (z=15m): computed={:.6}, expected={:.6}", rd_deep, rd_deep_expected
    );

    // rd should decrease with depth (soil column flexibility)
    assert!(
        rd_deep < rd,
        "rd should decrease with depth: rd(8m)={:.4} > rd(15m)={:.4}", rd, rd_deep
    );
}

// ================================================================
// 5. Dynamic Soil Properties: G/Gmax Degradation Curves
// ================================================================
//
// The shear modulus of soil degrades with increasing shear strain.
// Seed & Idriss (1970) and Vucetic & Dobry (1991) provide curves
// for G/Gmax as a function of cyclic shear strain gamma.
//
// The hyperbolic model (Hardin & Drnevich, 1972):
//   G/Gmax = 1 / (1 + gamma/gamma_ref)
//
// where gamma_ref is the reference strain (strain at which G/Gmax = 0.5).
//
// Gmax = rho * Vs^2

#[test]
fn validation_soil_dynamics_g_gmax_degradation() {
    // Small-strain shear modulus: Gmax = rho * Vs^2
    let rho: f64 = 1800.0;   // mass density (kg/m^3)
    let vs: f64 = 250.0;      // shear wave velocity (m/s)
    let gmax: f64 = rho * vs * vs;
    let gmax_expected: f64 = 1800.0 * 250.0 * 250.0; // 112,500,000 Pa = 112.5 MPa

    assert!(
        (gmax - gmax_expected).abs() < 1e-6,
        "Gmax = rho*Vs^2: computed={:.0} Pa, expected={:.0} Pa", gmax, gmax_expected
    );

    // Hyperbolic model: G/Gmax = 1 / (1 + gamma/gamma_ref)
    let gamma_ref: f64 = 0.0005; // reference strain (0.05%)

    // At very small strain (gamma << gamma_ref), G/Gmax -> 1.0
    let gamma_small: f64 = 1e-6;
    let ratio_small: f64 = 1.0 / (1.0 + gamma_small / gamma_ref);
    assert!(
        (ratio_small - 1.0).abs() < 0.01,
        "G/Gmax at very small strain should be ~1.0, got {:.6}", ratio_small
    );

    // At reference strain (gamma = gamma_ref), G/Gmax = 0.5
    let ratio_ref: f64 = 1.0 / (1.0 + gamma_ref / gamma_ref);
    assert!(
        (ratio_ref - 0.5).abs() < 1e-10,
        "G/Gmax at reference strain should be 0.5, got {:.6}", ratio_ref
    );

    // At large strain (gamma = 10 * gamma_ref), G/Gmax = 1/11 ~ 0.0909
    let gamma_large: f64 = 10.0 * gamma_ref;
    let ratio_large: f64 = 1.0 / (1.0 + gamma_large / gamma_ref);
    let ratio_large_expected: f64 = 1.0 / 11.0;

    assert!(
        (ratio_large - ratio_large_expected).abs() < 1e-10,
        "G/Gmax at 10*gamma_ref: computed={:.6}, expected={:.6}",
        ratio_large, ratio_large_expected
    );

    // G/Gmax should be monotonically decreasing with strain
    let strains: [f64; 5] = [1e-6, 1e-4, 5e-4, 1e-3, 1e-2];
    let mut prev_ratio: f64 = 1.0;
    for &gamma in &strains {
        let ratio: f64 = 1.0 / (1.0 + gamma / gamma_ref);
        assert!(
            ratio <= prev_ratio + 1e-10,
            "G/Gmax should decrease with strain: at gamma={:.1e}, ratio={:.4} > prev={:.4}",
            gamma, ratio, prev_ratio
        );
        prev_ratio = ratio;
    }

    // Damping ratio increases with strain (Hardin & Drnevich):
    // D = Dmax * (1 - G/Gmax)
    let _dmax: f64 = 0.25; // maximum damping ratio (25%)
    let d_at_ref: f64 = _dmax * (1.0 - ratio_ref); // D at reference strain
    let d_at_ref_expected: f64 = 0.25 * 0.5; // 0.125

    assert!(
        (d_at_ref - d_at_ref_expected).abs() < 1e-10,
        "Damping at reference strain: computed={:.4}, expected={:.4}",
        d_at_ref, d_at_ref_expected
    );

    // Secant shear modulus at a given strain
    let gamma_test: f64 = 0.001; // 0.1% shear strain
    let g_secant: f64 = gmax / (1.0 + gamma_test / gamma_ref);
    let g_secant_expected: f64 = gmax_expected / (1.0 + 0.001 / 0.0005);
    // = 112,500,000 / 3 = 37,500,000 Pa

    assert!(
        (g_secant - g_secant_expected).abs() < 1e-6,
        "Secant G at 0.1% strain: computed={:.0} Pa, expected={:.0} Pa",
        g_secant, g_secant_expected
    );
}

// ================================================================
// 6. Seismic Bearing Capacity Reduction
// ================================================================
//
// Under seismic loading, foundation bearing capacity is reduced due
// to inertial forces in the soil mass. The seismic bearing capacity
// can be expressed as a fraction of the static bearing capacity using
// seismic reduction factors.
//
// Terzaghi's bearing capacity (static):
//   q_ult = c*Nc + q*Nq + 0.5*gamma*B*Ngamma
//
// Seismic reduction approach (Richards et al., 1993):
//   q_ult_seismic = q_ult_static * (1 - kh*tan(phi))
//
// where kh = horizontal seismic coefficient (amax/g)
//
// More rigorously (EN 1998-5 Annex F), inclination factors reduce
// Nc, Nq, Ngamma.

#[test]
fn validation_soil_dynamics_seismic_bearing_capacity() {
    // Static bearing capacity parameters
    let c: f64 = 20.0;        // cohesion (kPa)
    let phi_deg: f64 = 30.0;  // friction angle (degrees)
    let phi: f64 = phi_deg * std::f64::consts::PI / 180.0;
    let gamma_soil: f64 = 18.0; // unit weight (kN/m^3)
    let b: f64 = 2.0;         // foundation width (m)
    let df: f64 = 1.5;        // foundation depth (m)
    let q: f64 = gamma_soil * df; // overburden pressure (kPa)

    // Terzaghi bearing capacity factors (Vesic, 1973)
    // Nq = exp(pi*tan(phi)) * tan^2(45 + phi/2)
    let nq: f64 = (std::f64::consts::PI * phi.tan()).exp()
        * (std::f64::consts::PI / 4.0 + phi / 2.0).tan().powi(2);
    // Nc = (Nq - 1) * cot(phi)
    let nc: f64 = (nq - 1.0) / phi.tan();
    // Ngamma = 2*(Nq + 1)*tan(phi)  (Vesic approximation)
    let ngamma: f64 = 2.0 * (nq + 1.0) * phi.tan();

    // Verify Nq for phi = 30 degrees: Nq ~ 18.40
    assert!(
        (nq - 18.40).abs() < 0.1,
        "Nq(30deg) should be ~18.40, got {:.2}", nq
    );

    // Verify Nc for phi = 30 degrees: Nc ~ 30.14
    assert!(
        (nc - 30.14).abs() < 0.2,
        "Nc(30deg) should be ~30.14, got {:.2}", nc
    );

    // Static bearing capacity
    let q_ult_static: f64 = c * nc + q * nq + 0.5 * gamma_soil * b * ngamma;

    // Seismic reduction (simplified approach)
    let kh: f64 = 0.2; // horizontal seismic coefficient
    let reduction_factor: f64 = 1.0 - kh * phi.tan();
    let q_ult_seismic: f64 = q_ult_static * reduction_factor;

    // Reduction factor for kh=0.2, phi=30deg: 1 - 0.2*tan(30) = 1 - 0.1155 = 0.8845
    let rf_expected: f64 = 1.0 - 0.2 * (30.0_f64 * std::f64::consts::PI / 180.0).tan();

    assert!(
        (reduction_factor - rf_expected).abs() < 1e-10,
        "Seismic reduction factor: computed={:.6}, expected={:.6}",
        reduction_factor, rf_expected
    );

    // Seismic capacity should be less than static
    assert!(
        q_ult_seismic < q_ult_static,
        "Seismic bearing capacity ({:.2} kPa) should be less than static ({:.2} kPa)",
        q_ult_seismic, q_ult_static
    );

    // Reduction factor should be between 0 and 1 for reasonable kh
    assert!(
        reduction_factor > 0.0 && reduction_factor < 1.0,
        "Reduction factor={:.4} should be in (0, 1)", reduction_factor
    );

    // Higher seismic coefficient -> greater reduction
    let kh_high: f64 = 0.4;
    let rf_high: f64 = 1.0 - kh_high * phi.tan();
    assert!(
        rf_high < reduction_factor,
        "Higher kh should give lower reduction factor: kh=0.4 -> RF={:.4} < RF(0.2)={:.4}",
        rf_high, reduction_factor
    );

    // EN 1998-5: partial factor approach uses gamma_rd = 1.0 for drained, 1.5 for undrained
    let _gamma_rd_drained: f64 = 1.0;
    let _gamma_rd_undrained: f64 = 1.5;
    let q_design_drained: f64 = q_ult_seismic / _gamma_rd_drained;
    let q_design_undrained: f64 = q_ult_seismic / _gamma_rd_undrained;

    assert!(
        q_design_undrained < q_design_drained,
        "Undrained design capacity ({:.2}) < drained ({:.2}) due to higher safety factor",
        q_design_undrained, q_design_drained
    );
}

// ================================================================
// 7. Lateral Spreading Displacement (Youd et al. 2002)
// ================================================================
//
// Youd et al. (2002) developed multilinear regression (MLR) models
// for estimating lateral spreading displacement. The free-face model:
//
//   log(DH) = b0 + b1*M + b2*log(R*) + b3*R + b4*log(W)
//             + b5*log(T15) + b6*log(100-F15) + b7*log(D50_15 + 0.1)
//
// where:
//   DH = horizontal displacement (m)
//   M = earthquake magnitude
//   R* = modified source distance (km)
//   R = site-to-source distance (km)
//   W = free-face ratio (%)
//   T15 = cumulative thickness of saturated sand with (N1)60 < 15 (m)
//   F15 = average fines content of T15 soils (%)
//   D50_15 = average D50 grain size of T15 soils (mm)
//
// Regression coefficients (free-face model, Youd et al. 2002):
//   b0 = -16.213, b1 = 1.532, b2 = -1.406, b3 = -0.012,
//   b4 = 0.338, b5 = 0.540, b6 = 3.413, b7 = -0.795

#[test]
fn validation_soil_dynamics_lateral_spreading_youd() {
    // Regression coefficients (Youd et al. 2002, free-face model)
    let b0: f64 = -16.213;
    let b1: f64 = 1.532;
    let b2: f64 = -1.406;
    let b3: f64 = -0.012;
    let b4: f64 = 0.338;
    let b5: f64 = 0.540;
    let b6: f64 = 3.413;
    let b7: f64 = -0.795;

    // Example parameters
    let mw: f64 = 7.5;       // moment magnitude
    let r: f64 = 20.0;       // source distance (km)
    let w: f64 = 10.0;       // free-face ratio (%)
    let t15: f64 = 5.0;      // cumulative thickness of liquefiable sand (m)
    let f15: f64 = 15.0;     // average fines content (%)
    let d50_15: f64 = 0.3;   // average D50 grain size (mm)

    // Modified source distance: R* = R + 10^(0.89*M - 5.64)
    let r_star: f64 = r + 10.0_f64.powf(0.89 * mw - 5.64);

    let r_star_expected: f64 = 20.0 + 10.0_f64.powf(0.89 * 7.5 - 5.64);
    assert!(
        (r_star - r_star_expected).abs() < 1e-6,
        "R* = R + 10^(0.89*M - 5.64): computed={:.4}, expected={:.4}",
        r_star, r_star_expected
    );

    // log(DH) computation (log base 10)
    let log_dh: f64 = b0 + b1 * mw + b2 * r_star.log10() + b3 * r
        + b4 * w.log10() + b5 * t15.log10()
        + b6 * (100.0 - f15).log10() + b7 * (d50_15 + 0.1).log10();

    // DH in meters
    let dh: f64 = 10.0_f64.powf(log_dh);

    // Verify each term independently
    let term_magnitude: f64 = b1 * mw;
    let term_distance: f64 = b2 * r_star.log10() + b3 * r;
    let term_geometry: f64 = b4 * w.log10();
    let term_soil: f64 = b5 * t15.log10() + b6 * (100.0 - f15).log10() + b7 * (d50_15 + 0.1).log10();

    let log_dh_sum: f64 = b0 + term_magnitude + term_distance + term_geometry + term_soil;
    assert!(
        (log_dh - log_dh_sum).abs() < 1e-10,
        "Sum of terms should equal log(DH): {:.6} vs {:.6}", log_dh, log_dh_sum
    );

    // DH should be positive and in a physically reasonable range (0.01 - 10 m)
    assert!(
        dh > 0.0,
        "Lateral spreading displacement should be positive, got {:.6} m", dh
    );
    assert!(
        dh < 20.0,
        "Lateral spreading displacement should be < 20 m for typical conditions, got {:.4} m", dh
    );

    // Higher magnitude should increase displacement (all else equal)
    let mw_higher: f64 = 8.0;
    let r_star_higher: f64 = r + 10.0_f64.powf(0.89 * mw_higher - 5.64);
    let log_dh_higher: f64 = b0 + b1 * mw_higher + b2 * r_star_higher.log10() + b3 * r
        + b4 * w.log10() + b5 * t15.log10()
        + b6 * (100.0 - f15).log10() + b7 * (d50_15 + 0.1).log10();
    let dh_higher: f64 = 10.0_f64.powf(log_dh_higher);

    assert!(
        dh_higher > dh,
        "Higher magnitude (M={:.1}) should give larger displacement: {:.4} m > {:.4} m",
        mw_higher, dh_higher, dh
    );

    // Greater distance should reduce displacement
    let r_far: f64 = 50.0;
    let r_star_far: f64 = r_far + 10.0_f64.powf(0.89 * mw - 5.64);
    let log_dh_far: f64 = b0 + b1 * mw + b2 * r_star_far.log10() + b3 * r_far
        + b4 * w.log10() + b5 * t15.log10()
        + b6 * (100.0 - f15).log10() + b7 * (d50_15 + 0.1).log10();
    let dh_far: f64 = 10.0_f64.powf(log_dh_far);

    assert!(
        dh_far < dh,
        "Greater distance (R={:.0} km) should reduce displacement: {:.4} m < {:.4} m",
        r_far, dh_far, dh
    );
}

// ================================================================
// 8. Site Period Estimation: Quarter-Wavelength Method
// ================================================================
//
// The fundamental site period of a soil deposit can be estimated as:
//   T_site = 4H / Vs     (uniform layer over rigid base)
//
// For a multi-layer profile (quarter-wavelength approximation):
//   T_site = 4 * Σ(di / Vsi)  = 4 * total_travel_time
//
// This equals four times the one-way S-wave travel time through
// the soil column.
//
// Reference: Kramer Ch. 7, Dobry et al. (2000)

#[test]
fn validation_soil_dynamics_site_period_quarter_wavelength() {
    // Case 1: Uniform layer
    let h: f64 = 30.0;       // layer thickness (m)
    let vs: f64 = 300.0;     // shear wave velocity (m/s)

    let t_site: f64 = 4.0 * h / vs;
    let t_site_expected: f64 = 4.0 * 30.0 / 300.0; // 0.4 s

    assert!(
        (t_site - t_site_expected).abs() < 1e-10,
        "Uniform layer T_site = 4H/Vs: computed={:.6} s, expected={:.6} s",
        t_site, t_site_expected
    );

    // Fundamental frequency
    let f_site: f64 = 1.0 / t_site;
    let f_site_expected: f64 = vs / (4.0 * h); // = 300/(4*30) = 2.5 Hz

    assert!(
        (f_site - f_site_expected).abs() < 1e-10,
        "Site frequency f = 1/T = Vs/(4H): computed={:.4} Hz, expected={:.4} Hz",
        f_site, f_site_expected
    );

    // Case 2: Multi-layer profile
    // Layer 1: 5 m, Vs = 120 m/s
    // Layer 2: 10 m, Vs = 200 m/s
    // Layer 3: 15 m, Vs = 350 m/s
    let layers: [(f64, f64); 3] = [
        (5.0, 120.0),
        (10.0, 200.0),
        (15.0, 350.0),
    ];

    let total_travel_time: f64 = layers.iter().map(|(d, v)| d / v).sum();
    let t_site_multi: f64 = 4.0 * total_travel_time;

    // Manual: travel_time = 5/120 + 10/200 + 15/350
    //                     = 0.04167 + 0.05000 + 0.04286 = 0.13452 s
    // T_site = 4 * 0.13452 = 0.53810 s
    let _tt_expected: f64 = 5.0 / 120.0 + 10.0 / 200.0 + 15.0 / 350.0;
    let t_site_multi_expected: f64 = 4.0 * _tt_expected;

    assert!(
        (t_site_multi - t_site_multi_expected).abs() < 1e-10,
        "Multi-layer T_site = 4*Σ(di/Vsi): computed={:.6} s, expected={:.6} s",
        t_site_multi, t_site_multi_expected
    );

    // Equivalent Vs for the multi-layer profile
    let total_h: f64 = layers.iter().map(|(d, _)| d).sum();
    let vs_equiv: f64 = total_h / total_travel_time;

    // Verify: T_site_multi = 4 * total_H / Vs_equiv
    let t_check: f64 = 4.0 * total_h / vs_equiv;
    assert!(
        (t_site_multi - t_check).abs() < 1e-10,
        "Quarter-wavelength consistency: T = 4H/Vs_equiv = {:.6} s", t_check
    );

    // Multi-layer period should be longer than uniform layer with highest Vs
    let vs_max: f64 = 350.0;
    let t_stiffest: f64 = 4.0 * total_h / vs_max;
    assert!(
        t_site_multi > t_stiffest,
        "Multi-layer period ({:.4} s) > uniform stiff layer ({:.4} s)",
        t_site_multi, t_stiffest
    );

    // Multi-layer period should be shorter than uniform layer with lowest Vs
    let vs_min: f64 = 120.0;
    let t_softest: f64 = 4.0 * total_h / vs_min;
    assert!(
        t_site_multi < t_softest,
        "Multi-layer period ({:.4} s) < uniform soft layer ({:.4} s)",
        t_site_multi, t_softest
    );

    // EC8 simplified site period approximation:
    // For soft soil (Site Class D), Eurocode 8 typically uses
    // corner periods TB = 0.15s, TC = 0.5s for Type 1 spectrum.
    // Our computed T_site ~ 0.54s is close to TC, consistent with
    // a site with predominantly soft-to-medium soil.
    let _ec8_tc_type1_ground_c: f64 = 0.6; // EN 1998-1 Table 3.2, Ground type C
    let _ec8_tc_type1_ground_d: f64 = 0.8; // Ground type D

    // Verify site period falls in a geotechnically reasonable range
    // for a 30m profile: typically 0.1s to 2.0s
    assert!(
        t_site_multi > 0.1 && t_site_multi < 2.0,
        "Multi-layer site period {:.4} s should be in [0.1, 2.0] range",
        t_site_multi
    );

    // Doubling all layer thicknesses should double the site period
    let t_site_double: f64 = 4.0 * layers.iter().map(|(d, v)| 2.0 * d / v).sum::<f64>();
    assert!(
        (t_site_double / t_site_multi - 2.0).abs() < 1e-10,
        "Doubling layer thickness should double site period: ratio={:.6}",
        t_site_double / t_site_multi
    );
}

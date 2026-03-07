/// Validation: Advanced Structural Dynamics Formulas (Pure Formula Verification)
///
/// References:
///   - Chopra, "Dynamics of Structures", 5th Ed.
///   - Clough & Penzien, "Dynamics of Structures", 3rd Ed.
///   - Newmark, "A Method of Computation for Structural Dynamics" (1959)
///   - Den Hartog, "Mechanical Vibrations", 4th Ed.
///
/// Tests verify advanced dynamics formulas without calling the solver.
///   1. Newmark-beta method: single step verification
///   2. Wilson-theta method: single step verification
///   3. Duhamel integral: impulse response of SDOF
///   4. Transmissibility ratio
///   5. Tuned mass damper: optimal parameters
///   6. Base isolation: period shift
///   7. Floor response spectrum amplification factor
///   8. Rayleigh damping coefficients from two target frequencies

mod helpers;

use std::f64::consts::PI;

// ================================================================
// 1. Newmark-Beta Method: Single Step Verification
// ================================================================
//
// The Newmark-beta method with beta=1/4 and gamma=1/2 (average
// acceleration method, unconditionally stable).
//
// For SDOF system: m*a + c*v + k*u = f(t)
//
// Newmark update equations:
//   u_{n+1} = u_n + dt*v_n + dt^2*((0.5-beta)*a_n + beta*a_{n+1})
//   v_{n+1} = v_n + dt*((1-gamma)*a_n + gamma*a_{n+1})
//
// Effective stiffness: k_eff = k + gamma/(beta*dt)*c + 1/(beta*dt^2)*m
// Effective force: f_eff = f_{n+1} + m*(...) + c*(...)
//
// For a simple case: m=1, c=0, k=100, f=0 (free vibration),
// u_0=1, v_0=0, a_0=-k*u_0/m = -100
//
// Reference: Chopra, Ch. 5; Newmark (1959)

#[test]
fn validation_dynamics_newmark_beta_single_step() {
    let m: f64 = 1.0;     // kg
    let c: f64 = 0.0;     // N*s/m (undamped)
    let k: f64 = 100.0;   // N/m
    let dt: f64 = 0.1;    // s, time step

    let beta: f64 = 0.25;
    let gamma: f64 = 0.5;

    // Initial conditions
    let u0 = 1.0;    // m
    let v0 = 0.0;    // m/s
    let a0 = -k * u0 / m; // = -100 m/s^2

    // External force at t = dt
    let f1 = 0.0;

    // Effective stiffness
    let k_eff = k + gamma / (beta * dt) * c + 1.0 / (beta * dt * dt) * m;
    // = 100 + 0 + 1/(0.25*0.01) = 100 + 400 = 500
    let k_eff_expected = 500.0;
    assert!(
        (k_eff - k_eff_expected).abs() / k_eff_expected < 1e-10,
        "k_eff: computed={:.2}, expected={:.2}",
        k_eff, k_eff_expected
    );

    // Effective force at step n+1
    let f_eff = f1
        + m * (1.0 / (beta * dt * dt) * u0 + 1.0 / (beta * dt) * v0 + (1.0 / (2.0 * beta) - 1.0) * a0)
        + c * (gamma / (beta * dt) * u0 + (gamma / beta - 1.0) * v0 + dt * (gamma / (2.0 * beta) - 1.0) * a0);

    // For c=0: f_eff = 0 + m*(400*1 + 0 + (2-1)*(-100)) + 0
    //        = 1*(400 + 0 - 100) = 300
    let f_eff_expected = 300.0;
    assert!(
        (f_eff - f_eff_expected).abs() / f_eff_expected < 1e-10,
        "f_eff: computed={:.2}, expected={:.2}",
        f_eff, f_eff_expected
    );

    // Displacement at t = dt
    let u1 = f_eff / k_eff;
    let u1_expected = 300.0 / 500.0; // = 0.6
    assert!(
        (u1 - u1_expected).abs() / u1_expected < 1e-10,
        "u1: computed={:.6}, expected={:.6}",
        u1, u1_expected
    );

    // Acceleration at t = dt
    let a1 = 1.0 / (beta * dt * dt) * (u1 - u0) - 1.0 / (beta * dt) * v0 - (1.0 / (2.0 * beta) - 1.0) * a0;
    // = 400*(0.6-1) - 0 - 1*(-100) = 400*(-0.4) + 100 = -160 + 100 = -60
    let a1_expected = -60.0;
    assert!(
        (a1 - a1_expected).abs() / a1_expected.abs() < 1e-10,
        "a1: computed={:.6}, expected={:.6}",
        a1, a1_expected
    );

    // Velocity at t = dt
    let v1 = v0 + dt * ((1.0 - gamma) * a0 + gamma * a1);
    // = 0 + 0.1*(0.5*(-100) + 0.5*(-60)) = 0.1*(-50 - 30) = -8
    let v1_expected = -8.0;
    assert!(
        (v1 - v1_expected).abs() / v1_expected.abs() < 1e-10,
        "v1: computed={:.6}, expected={:.6}",
        v1, v1_expected
    );

    // Verify equilibrium at t = dt: m*a1 + k*u1 = f1
    let residual = m * a1 + k * u1 - f1;
    // = 1*(-60) + 100*0.6 - 0 = -60 + 60 = 0
    assert!(
        residual.abs() < 1e-10,
        "Equilibrium residual: {:.6e}",
        residual
    );
}

// ================================================================
// 2. Wilson-Theta Method: Single Step (theta = 1.4)
// ================================================================
//
// The Wilson-theta method extends the acceleration linearly over
// a time interval theta*dt (theta >= 1.0 for stability, typically 1.4).
//
// For SDOF: m*a + c*v + k*u = f(t)
//
// Effective stiffness: k_hat = k + 3c/(theta*dt) + 6m/(theta*dt)^2
// The method computes displacement at t + theta*dt, then
// interpolates back to t + dt.
//
// Reference: Wilson, Farhoomand & Bathe (1972); Chopra, Ch. 5

#[test]
fn validation_dynamics_wilson_theta_single_step() {
    let m: f64 = 2.0;      // kg
    let c: f64 = 0.0;      // undamped
    let k: f64 = 200.0;    // N/m
    let dt: f64 = 0.05;    // s
    let theta: f64 = 1.4;

    // Initial conditions
    let u0 = 0.5;
    let v0 = 0.0;
    let a0 = -k * u0 / m; // = -50 m/s^2

    // Force at t and t+theta*dt
    let _f0 = 0.0;
    let f_theta = 0.0; // free vibration

    // Extended time step
    let dt_theta = theta * dt; // 0.07 s

    // Effective stiffness for Wilson-theta
    let k_hat = k + 3.0 * c / dt_theta + 6.0 * m / (dt_theta * dt_theta);
    // = 200 + 0 + 6*2/(0.07^2) = 200 + 12/0.0049 = 200 + 2448.98
    let k_hat_expected = 200.0 + 6.0 * 2.0 / (0.07 * 0.07);
    assert!(
        (k_hat - k_hat_expected).abs() / k_hat_expected < 1e-10,
        "Wilson k_hat: computed={:.2}, expected={:.2}",
        k_hat, k_hat_expected
    );

    // Effective force at t + theta*dt
    let f_hat = f_theta
        + m * (6.0 / (dt_theta * dt_theta) * u0 + 6.0 / dt_theta * v0 + 2.0 * a0)
        + c * (3.0 / dt_theta * u0 + 2.0 * v0 + dt_theta / 2.0 * a0);
    // = 0 + 2*(6/(0.0049)*0.5 + 0 + 2*(-50)) + 0
    // = 2*(612.24*0.5 - 100) ... let me compute step by step
    // 6/dt_theta^2 * u0 = 6/0.0049 * 0.5 = 612.2449
    // 6/dt_theta * v0 = 0
    // 2 * a0 = -100
    // sum = 512.2449
    // m * sum = 2 * 512.2449 = 1024.4898
    let term1 = 6.0 / (dt_theta * dt_theta) * u0;
    let term2 = 6.0 / dt_theta * v0;
    let term3 = 2.0 * a0;
    let f_hat_expected = m * (term1 + term2 + term3);
    assert!(
        (f_hat - f_hat_expected).abs() / f_hat_expected.abs() < 1e-10,
        "Wilson f_hat: computed={:.4}, expected={:.4}",
        f_hat, f_hat_expected
    );

    // Displacement at t + theta*dt
    let u_theta = f_hat / k_hat;

    // Acceleration at t + theta*dt
    let a_theta = 6.0 / (dt_theta * dt_theta) * (u_theta - u0) - 6.0 / dt_theta * v0 - 2.0 * a0;

    // Interpolate acceleration back to t + dt:
    // a1 = a0 + (a_theta - a0) / theta
    let a1 = a0 + (a_theta - a0) / theta;

    // Update velocity and displacement at t + dt
    let _v1 = v0 + dt / 2.0 * (a0 + a1);
    let u1 = u0 + dt * v0 + dt * dt / 6.0 * (2.0 * a0 + a1);

    // Verify that displacement is physically reasonable
    // Free vibration from u0=0.5, period T = 2*pi*sqrt(m/k) = 2*pi*sqrt(0.01) = 0.628 s
    let period = 2.0 * PI * (m / k).sqrt();
    let period_expected = 2.0 * PI * 0.1; // ~0.6283 s
    assert!(
        (period - period_expected).abs() / period_expected < 1e-10,
        "Period: computed={:.4}, expected={:.4}",
        period, period_expected
    );

    // At t=dt=0.05, the displacement should still be positive
    // (we haven't reached quarter period yet, dt/T ~ 0.08)
    assert!(
        u1 > 0.0 && u1 < u0,
        "Wilson u1={:.6}: should be between 0 and u0={:.6}",
        u1, u0
    );

    // Verify equilibrium at t+dt (approximately)
    let eq_residual = m * a1 + k * u1;
    // Should be close to zero for free vibration
    assert!(
        eq_residual.abs() < k * u0 * 0.05,
        "Wilson equilibrium residual: {:.6} (tol: {:.6})",
        eq_residual, k * u0 * 0.05
    );
}

// ================================================================
// 3. Duhamel Integral: Impulse Response of SDOF
// ================================================================
//
// For an undamped SDOF system subjected to an impulse I at t=0:
//   u(t) = I / (m * omega_n) * sin(omega_n * t)
// where omega_n = sqrt(k/m) is the natural frequency.
//
// For a damped SDOF system:
//   u(t) = I / (m * omega_d) * exp(-zeta*omega_n*t) * sin(omega_d*t)
// where omega_d = omega_n * sqrt(1 - zeta^2) is the damped frequency.
//
// Reference: Chopra, Ch. 4; Clough & Penzien, Ch. 3

#[test]
fn validation_dynamics_duhamel_impulse_response() {
    let m: f64 = 5.0;      // kg
    let k: f64 = 500.0;    // N/m
    let impulse: f64 = 10.0; // N*s

    // Natural frequency
    let omega_n = (k / m).sqrt();
    let omega_n_expected = 10.0; // sqrt(100) = 10 rad/s
    assert!(
        (omega_n - omega_n_expected).abs() / omega_n_expected < 1e-10,
        "omega_n: computed={:.4}, expected={:.4}",
        omega_n, omega_n_expected
    );

    // Undamped response at t = pi/(2*omega_n) (quarter period: maximum)
    let t_quarter = PI / (2.0 * omega_n);
    let u_max_undamped = impulse / (m * omega_n) * (omega_n * t_quarter).sin();
    // = 10/(5*10) * sin(pi/2) = 0.2 * 1.0 = 0.2
    let u_max_expected = 0.2;
    assert!(
        (u_max_undamped - u_max_expected).abs() / u_max_expected < 1e-10,
        "Undamped max response: computed={:.6}, expected={:.6}",
        u_max_undamped, u_max_expected
    );

    // Undamped response at t = pi/omega_n (half period: zero crossing)
    let t_half = PI / omega_n;
    let u_half_undamped = impulse / (m * omega_n) * (omega_n * t_half).sin();
    assert!(
        u_half_undamped.abs() < 1e-10,
        "Response at half period: {:.6e} (should be ~0)",
        u_half_undamped
    );

    // Damped response
    let zeta: f64 = 0.05; // 5% damping
    let omega_d = omega_n * (1.0 - zeta * zeta).sqrt();
    let omega_d_expected = 10.0 * (1.0 - 0.0025_f64).sqrt();
    assert!(
        (omega_d - omega_d_expected).abs() / omega_d_expected < 1e-10,
        "omega_d: computed={:.6}, expected={:.6}",
        omega_d, omega_d_expected
    );

    // Damped response at t = pi/(2*omega_d)
    let t_d_quarter = PI / (2.0 * omega_d);
    let u_damped = impulse / (m * omega_d)
        * (-zeta * omega_n * t_d_quarter).exp()
        * (omega_d * t_d_quarter).sin();

    // The damped response should be less than undamped
    assert!(
        u_damped < u_max_undamped,
        "Damped ({:.6}) < Undamped ({:.6})",
        u_damped, u_max_undamped
    );
    assert!(
        u_damped > 0.0,
        "Damped response should be positive at quarter period: {:.6}",
        u_damped
    );

    // Peak amplitude decays exponentially: ratio of successive peaks
    // u_n/u_{n+1} = exp(2*pi*zeta/sqrt(1-zeta^2))
    let log_dec = 2.0 * PI * zeta / (1.0 - zeta * zeta).sqrt();
    let peak_ratio = log_dec.exp();
    let peak_ratio_expected = (2.0 * PI * 0.05 / (0.9975_f64).sqrt()).exp();
    assert!(
        (peak_ratio - peak_ratio_expected).abs() / peak_ratio_expected < 1e-6,
        "Peak ratio: computed={:.6}, expected={:.6}",
        peak_ratio, peak_ratio_expected
    );
}

// ================================================================
// 4. Transmissibility: Force or Displacement Transmission
// ================================================================
//
// Transmissibility T_r is the ratio of transmitted force (or displacement)
// to applied force (or base displacement):
//
//   T_r = sqrt((1 + (2*zeta*r)^2) / ((1-r^2)^2 + (2*zeta*r)^2))
//
// where r = omega/omega_n is the frequency ratio.
//
// Key properties:
//   - At r=0: T_r = 1 (quasi-static)
//   - At r=1: T_r = sqrt(1+4*zeta^2)/(2*zeta) (resonance)
//   - At r=sqrt(2): T_r = 1 (crossover point, independent of damping)
//   - For r >> 1: T_r → 2*zeta*r / r^2 → 0 (isolation region)
//
// Reference: Den Hartog, Ch. 2; Chopra, Ch. 3

#[test]
fn validation_dynamics_transmissibility() {
    let zeta: f64 = 0.05; // 5% damping

    // Transmissibility function
    let transmissibility = |r: f64, z: f64| -> f64 {
        let num = 1.0 + (2.0 * z * r).powi(2);
        let den = (1.0 - r * r).powi(2) + (2.0 * z * r).powi(2);
        (num / den).sqrt()
    };

    // At r = 0: T = 1 (quasi-static)
    let t_0 = transmissibility(0.0, zeta);
    assert!(
        (t_0 - 1.0).abs() < 1e-10,
        "T(r=0) = {:.6}, expected 1.0",
        t_0
    );

    // At r = 1 (resonance): T = sqrt(1 + 4*zeta^2) / (2*zeta)
    let t_resonance = transmissibility(1.0, zeta);
    let t_resonance_expected = (1.0 + 4.0 * zeta * zeta).sqrt() / (2.0 * zeta);
    assert!(
        (t_resonance - t_resonance_expected).abs() / t_resonance_expected < 1e-10,
        "T(r=1): computed={:.4}, expected={:.4}",
        t_resonance, t_resonance_expected
    );
    // For zeta=0.05: T ~ sqrt(1.01)/0.1 ~ 10.05
    assert!(
        t_resonance > 9.0 && t_resonance < 11.0,
        "T at resonance should be ~10 for zeta=0.05: T={:.4}",
        t_resonance
    );

    // At r = sqrt(2): T = 1 (crossover point, independent of damping)
    let r_cross = 2.0_f64.sqrt();
    let t_cross = transmissibility(r_cross, zeta);
    // T = sqrt((1+8*zeta^2) / ((1-2)^2 + 8*zeta^2))
    //   = sqrt((1+8z^2) / (1+8z^2)) = 1
    let t_cross_expected = 1.0;
    assert!(
        (t_cross - t_cross_expected).abs() < 1e-10,
        "T(r=sqrt(2)): computed={:.6}, expected={:.6}",
        t_cross, t_cross_expected
    );

    // Verify crossover is independent of damping
    let t_cross_z10 = transmissibility(r_cross, 0.10);
    let t_cross_z20 = transmissibility(r_cross, 0.20);
    assert!(
        (t_cross_z10 - 1.0).abs() < 1e-10,
        "T(r=sqrt(2), zeta=0.10): {:.6}",
        t_cross_z10
    );
    assert!(
        (t_cross_z20 - 1.0).abs() < 1e-10,
        "T(r=sqrt(2), zeta=0.20): {:.6}",
        t_cross_z20
    );

    // For r >> 1 (isolation region): T < 1
    let t_high = transmissibility(3.0, zeta);
    assert!(
        t_high < 1.0,
        "T(r=3) = {:.4} should be < 1 (isolation)",
        t_high
    );
}

// ================================================================
// 5. Tuned Mass Damper: Optimal Parameters
// ================================================================
//
// For a TMD (tuned mass damper) with mass ratio mu = m_d / m_s:
//
// Optimal tuning (Den Hartog formulas):
//   f_opt = 1 / (1 + mu)          (frequency ratio)
//   zeta_opt = sqrt(3*mu / (8*(1+mu)^3))  (optimal TMD damping)
//
// The TMD reduces the peak dynamic amplification of the main structure.
//
// Reference: Den Hartog, Ch. 3; Chopra, Sec. 14.4

#[test]
fn validation_dynamics_tuned_mass_damper() {
    // Mass ratio mu = m_TMD / m_structure
    let mu: f64 = 0.02; // 2% mass ratio (typical)

    // Optimal frequency ratio
    let f_opt = 1.0 / (1.0 + mu);
    let f_opt_expected = 1.0 / 1.02;
    assert!(
        (f_opt - f_opt_expected).abs() / f_opt_expected < 1e-10,
        "f_opt: computed={:.6}, expected={:.6}",
        f_opt, f_opt_expected
    );

    // Optimal damping ratio for TMD
    let zeta_opt = (3.0 * mu / (8.0 * (1.0 + mu).powi(3))).sqrt();
    let zeta_opt_expected = (3.0 * 0.02 / (8.0 * 1.02_f64.powi(3))).sqrt();
    assert!(
        (zeta_opt - zeta_opt_expected).abs() / zeta_opt_expected < 1e-10,
        "zeta_opt: computed={:.6}, expected={:.6}",
        zeta_opt, zeta_opt_expected
    );

    // For mu = 0.02: zeta_opt ~ sqrt(0.06/8.489) ~ sqrt(0.00707) ~ 0.0841
    assert!(
        zeta_opt > 0.05 && zeta_opt < 0.15,
        "TMD optimal damping should be ~8% for mu=2%: zeta={:.4}",
        zeta_opt
    );

    // Verify for a larger mass ratio: mu = 0.05
    let mu2: f64 = 0.05;
    let f_opt2 = 1.0 / (1.0 + mu2);
    let zeta_opt2 = (3.0 * mu2 / (8.0 * (1.0 + mu2).powi(3))).sqrt();

    // Larger mu → more detuning (lower f_opt) and higher damping
    assert!(
        f_opt2 < f_opt,
        "Larger mu → lower f_opt: {:.6} < {:.6}",
        f_opt2, f_opt
    );
    assert!(
        zeta_opt2 > zeta_opt,
        "Larger mu → higher zeta_opt: {:.6} > {:.6}",
        zeta_opt2, zeta_opt
    );

    // Peak amplification with optimal TMD: DAF ~ sqrt(2/mu) (approximate)
    let daf_approx = (2.0 / mu).sqrt();
    // For mu = 0.02: DAF ~ sqrt(100) = 10
    let daf_expected = 10.0;
    assert!(
        (daf_approx - daf_expected).abs() / daf_expected < 1e-10,
        "DAF with TMD: computed={:.4}, expected={:.4}",
        daf_approx, daf_expected
    );
}

// ================================================================
// 6. Base Isolation: Period Shift
// ================================================================
//
// The isolated period of a structure on base isolators:
//   T_iso = 2*pi * sqrt(m / k_iso)
//
// where k_iso is the total isolation system stiffness.
//
// The effective period is related to the fixed-base period:
//   T_iso / T_fixed = sqrt(1 + k_fixed/k_iso)
//
// For effective isolation: T_iso / T_fixed > 3 (typically)
//
// Reference: Chopra, Ch. 21; Naeim & Kelly, "Design of Seismic
//   Isolated Structures"

#[test]
fn validation_dynamics_base_isolation_period_shift() {
    let m: f64 = 500_000.0;    // kg (500 tonnes, typical building)
    let k_fixed: f64 = 2e8;    // N/m (fixed-base stiffness)
    let k_iso: f64 = 5e6;      // N/m (isolation stiffness, much softer)

    // Fixed-base period
    let t_fixed = 2.0 * PI * (m / k_fixed).sqrt();
    // = 2*pi*sqrt(500000/2e8) = 2*pi*sqrt(0.0025) = 2*pi*0.05 = 0.3142 s
    let t_fixed_expected = 2.0 * PI * 0.05;
    assert!(
        (t_fixed - t_fixed_expected).abs() / t_fixed_expected < 1e-10,
        "T_fixed: computed={:.4}, expected={:.4}",
        t_fixed, t_fixed_expected
    );

    // Isolated period
    let t_iso = 2.0 * PI * (m / k_iso).sqrt();
    // = 2*pi*sqrt(500000/5e6) = 2*pi*sqrt(0.1) = 2*pi*0.3162 = 1.987 s
    let t_iso_expected = 2.0 * PI * (0.1_f64).sqrt();
    assert!(
        (t_iso - t_iso_expected).abs() / t_iso_expected < 1e-10,
        "T_iso: computed={:.4}, expected={:.4}",
        t_iso, t_iso_expected
    );

    // Period ratio: T_iso / T_fixed should be > 3 for effective isolation
    let period_ratio = t_iso / t_fixed;
    assert!(
        period_ratio > 3.0,
        "Period ratio {:.2} should be > 3 for effective isolation",
        period_ratio
    );

    // Verify the relationship: T_iso/T_fixed = sqrt(k_fixed/k_iso)
    // (when the structure is rigid relative to the isolator)
    let ratio_from_stiffness = (k_fixed / k_iso).sqrt();
    assert!(
        (period_ratio - ratio_from_stiffness).abs() / ratio_from_stiffness < 1e-10,
        "Period ratio from stiffness: computed={:.4}, expected={:.4}",
        period_ratio, ratio_from_stiffness
    );

    // Force reduction: the base shear is reduced by approximately
    // (T_fixed/T_iso)^2 in the displacement-sensitive range of the spectrum
    let force_reduction = (t_fixed / t_iso).powi(2);
    // = (0.3142/1.987)^2 = (0.1581)^2 = 0.025
    assert!(
        force_reduction < 0.1,
        "Force reduction: {:.4} should be < 0.1 (significant reduction)",
        force_reduction
    );
}

// ================================================================
// 7. Floor Response Spectrum Amplification Factor
// ================================================================
//
// The amplification of floor acceleration relative to ground
// acceleration depends on the building's natural period and
// the floor level.
//
// For a building with fundamental period T_b and a piece of
// equipment with period T_e at floor level:
//
// Maximum amplification at a floor (simplified formula):
//   A_floor = 1 + 2 * (z/H)     (linear approximation)
// where z = height of floor, H = total height.
//
// This gives A_floor = 1 at ground (z=0) and A_floor = 3 at roof (z=H).
//
// More refined (ASCE 7):
//   a_p * S_DS * (1 + 2*z/h) / (Rp/Ip)
//
// The floor amplification factor A_f:
//   A_f = S_a(T_e, zeta_e) / S_a(T_b, zeta_b) * phi(z)
// where phi(z) is the mode shape at height z.
//
// Reference: ASCE 7-22 Chapter 13; Chopra, Ch. 13

#[test]
fn validation_dynamics_floor_response_amplification() {
    // Building parameters
    let h_total: f64 = 30.0; // m, total building height
    let n_floors: i32 = 10;

    // Floor amplification factor (ASCE 7 simplified)
    // A_f(z) = 1 + 2 * z / H
    let amplification = |z: f64| -> f64 {
        1.0 + 2.0 * z / h_total
    };

    // Ground floor (z = 0): A = 1.0
    let a_ground = amplification(0.0);
    assert!(
        (a_ground - 1.0).abs() < 1e-10,
        "Ground amplification: {:.4}, expected 1.0",
        a_ground
    );

    // Roof (z = H): A = 3.0
    let a_roof = amplification(h_total);
    assert!(
        (a_roof - 3.0).abs() < 1e-10,
        "Roof amplification: {:.4}, expected 3.0",
        a_roof
    );

    // Mid-height (z = H/2): A = 2.0
    let a_mid = amplification(h_total / 2.0);
    assert!(
        (a_mid - 2.0).abs() < 1e-10,
        "Mid-height amplification: {:.4}, expected 2.0",
        a_mid
    );

    // Verify amplification increases monotonically with height
    let floor_height = h_total / n_floors as f64;
    let mut prev_a = 0.0;
    for i in 0..=n_floors {
        let z = i as f64 * floor_height;
        let a = amplification(z);
        assert!(
            a >= prev_a,
            "Amplification at z={:.1} ({:.4}) should be >= prev ({:.4})",
            z, a, prev_a
        );
        prev_a = a;
    }

    // Component amplification factor a_p (ASCE 7 Table 13.5-1)
    // For flexible equipment: a_p = 2.5
    // For rigid equipment: a_p = 1.0
    let a_p_flexible: f64 = 2.5;
    let a_p_rigid: f64 = 1.0;
    let s_ds: f64 = 1.0;  // Short-period spectral acceleration (g)
    let rp: f64 = 2.5;    // Response modification factor
    let ip: f64 = 1.0;    // Importance factor

    // Design force at roof for flexible equipment
    let fp_roof = a_p_flexible * s_ds * (1.0 + 2.0 * h_total / h_total) / (rp / ip);
    // = 2.5 * 1.0 * 3.0 / 2.5 = 3.0
    let fp_roof_expected = 3.0;
    assert!(
        (fp_roof - fp_roof_expected).abs() / fp_roof_expected < 1e-10,
        "Fp at roof (flexible): computed={:.4}, expected={:.4}",
        fp_roof, fp_roof_expected
    );

    // Design force at ground for rigid equipment
    let fp_ground = a_p_rigid * s_ds * (1.0 + 2.0 * 0.0 / h_total) / (rp / ip);
    // = 1.0 * 1.0 * 1.0 / 2.5 = 0.4
    let fp_ground_expected = 0.4;
    assert!(
        (fp_ground - fp_ground_expected).abs() / fp_ground_expected < 1e-10,
        "Fp at ground (rigid): computed={:.4}, expected={:.4}",
        fp_ground, fp_ground_expected
    );
}

// ================================================================
// 8. Rayleigh Damping Coefficients from Two Target Frequencies
// ================================================================
//
// Rayleigh (proportional) damping: C = a0*M + a1*K
//
// The damping ratio at frequency omega_i:
//   zeta_i = a0/(2*omega_i) + a1*omega_i/2
//
// Given target damping ratio zeta at two frequencies omega_1 and omega_2:
//   a0 = 2*zeta*omega_1*omega_2 / (omega_1 + omega_2)
//   a1 = 2*zeta / (omega_1 + omega_2)
//
// Reference: Chopra, Sec. 11.4; Clough & Penzien, Sec. 12.4

#[test]
fn validation_dynamics_rayleigh_damping_coefficients() {
    let zeta: f64 = 0.05; // 5% target damping ratio

    // Two target frequencies
    let f1: f64 = 2.0;   // Hz (first mode)
    let f2: f64 = 10.0;  // Hz (higher mode)
    let omega_1 = 2.0 * PI * f1;
    let omega_2 = 2.0 * PI * f2;

    // Rayleigh damping coefficients
    let a0 = 2.0 * zeta * omega_1 * omega_2 / (omega_1 + omega_2);
    let a1 = 2.0 * zeta / (omega_1 + omega_2);

    // Verify: a0 = 2*0.05*(4*pi)*(20*pi)/(24*pi) = 2*0.05*80*pi^2/(24*pi)
    //        = 8*pi^2/(24*pi) * 0.1 = pi/3 * 0.1 ... let me just compute numerically
    let a0_expected = 2.0 * 0.05 * (2.0 * PI * 2.0) * (2.0 * PI * 10.0)
        / (2.0 * PI * 2.0 + 2.0 * PI * 10.0);
    assert!(
        (a0 - a0_expected).abs() / a0_expected < 1e-10,
        "a0: computed={:.6}, expected={:.6}",
        a0, a0_expected
    );

    let a1_expected = 2.0 * 0.05 / (2.0 * PI * 2.0 + 2.0 * PI * 10.0);
    assert!(
        (a1 - a1_expected).abs() / a1_expected < 1e-10,
        "a1: computed={:.6e}, expected={:.6e}",
        a1, a1_expected
    );

    // Verify that damping is exactly zeta at both target frequencies
    let zeta_at_f1 = a0 / (2.0 * omega_1) + a1 * omega_1 / 2.0;
    assert!(
        (zeta_at_f1 - zeta).abs() / zeta < 1e-10,
        "zeta at f1: computed={:.6}, expected={:.6}",
        zeta_at_f1, zeta
    );

    let zeta_at_f2 = a0 / (2.0 * omega_2) + a1 * omega_2 / 2.0;
    assert!(
        (zeta_at_f2 - zeta).abs() / zeta < 1e-10,
        "zeta at f2: computed={:.6}, expected={:.6}",
        zeta_at_f2, zeta
    );

    // Damping at intermediate frequency (should be close to but less than zeta)
    let f_mid = (f1 * f2).sqrt(); // geometric mean ~ 4.47 Hz
    let omega_mid = 2.0 * PI * f_mid;
    let zeta_mid = a0 / (2.0 * omega_mid) + a1 * omega_mid / 2.0;

    // At the geometric mean of the two target frequencies,
    // the Rayleigh damping is slightly less than the target
    // (the damping curve has a minimum between the two target frequencies)
    assert!(
        zeta_mid < zeta,
        "Damping at midpoint ({:.6}) should be less than target ({:.6})",
        zeta_mid, zeta
    );
    assert!(
        zeta_mid > 0.0,
        "Damping should be positive: {:.6}",
        zeta_mid
    );

    // Damping outside the range should exceed zeta
    let f_high = 20.0; // Hz, above the target range
    let omega_high = 2.0 * PI * f_high;
    let zeta_high = a0 / (2.0 * omega_high) + a1 * omega_high / 2.0;
    assert!(
        zeta_high > zeta,
        "Damping at f={} Hz ({:.6}) should exceed target ({:.6})",
        f_high, zeta_high, zeta
    );
}

/// Compute Rayleigh damping coefficients from two target frequencies and damping ratio.
/// omega1, omega2: circular frequencies (rad/s)
/// xi: damping ratio (typically 0.02-0.05)
/// Returns (a0, a1) such that C = a0*M + a1*K
pub fn rayleigh_coefficients(omega1: f64, omega2: f64, xi: f64) -> (f64, f64) {
    let a0 = 2.0 * xi * omega1 * omega2 / (omega1 + omega2);
    let a1 = 2.0 * xi / (omega1 + omega2);
    (a0, a1)
}

/// Assemble Rayleigh damping matrix C = a0*M + a1*K.
/// m, k: mass and stiffness matrices (n x n, row-major dense)
pub fn rayleigh_damping_matrix(m: &[f64], k: &[f64], n: usize, a0: f64, a1: f64) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n * n {
        c[i] = a0 * m[i] + a1 * k[i];
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rayleigh_coefficients_symmetry() {
        let xi = 0.05;
        let omega1 = 10.0;
        let omega2 = 30.0;
        let (a0, a1) = rayleigh_coefficients(omega1, omega2, xi);

        // At omega1: xi = a0/(2*omega1) + a1*omega1/2
        let xi_at_1 = a0 / (2.0 * omega1) + a1 * omega1 / 2.0;
        assert!((xi_at_1 - xi).abs() < 1e-12, "xi at omega1 = {}", xi_at_1);

        // At omega2: xi = a0/(2*omega2) + a1*omega2/2
        let xi_at_2 = a0 / (2.0 * omega2) + a1 * omega2 / 2.0;
        assert!((xi_at_2 - xi).abs() < 1e-12, "xi at omega2 = {}", xi_at_2);
    }

    #[test]
    fn test_rayleigh_damping_matrix() {
        let m = vec![2.0, 0.0, 0.0, 3.0];
        let k = vec![10.0, -5.0, -5.0, 10.0];
        let c = rayleigh_damping_matrix(&m, &k, 2, 0.5, 0.1);
        // c[0] = 0.5*2 + 0.1*10 = 2.0
        assert!((c[0] - 2.0).abs() < 1e-12);
        // c[1] = 0.5*0 + 0.1*(-5) = -0.5
        assert!((c[1] - (-0.5)).abs() < 1e-12);
        // c[3] = 0.5*3 + 0.1*10 = 2.5
        assert!((c[3] - 2.5).abs() < 1e-12);
    }
}

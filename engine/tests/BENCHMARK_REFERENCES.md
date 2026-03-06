# Benchmark References for Proper Test Reproductions

## Curved Beams

### VM18 — Out-of-Plane Bending of a Curved Bar
- **Source**: ANSYS Verification Manual, VM18
- **URL**: https://examples.mapdl.docs.pyansys.com/verif-manual/vm-018.html
- **Problem**: Quarter-circle cantilever bar with out-of-plane load at free end
- **Parameters**: R = 100 in, cross-section 2×2 in (square), E = 30×10⁶ psi, ν = 0.0, F = 50 lb (out-of-plane)
- **Reference Answer**: δ_B = -2.648 in (displacement at free end in load direction)
- **Analytical Source**: Timoshenko, *Strength of Materials*, Part II, 3rd Ed., p. 412, eq. (c)
- **Formula**: δ = F·R³/(4EI) · (π + 2π·EI/(GJ))  where GJ = GJ_torsion
- **Note**: This is a 3D problem with out-of-plane loading. Our solver handles this via curved beam segment subdivision + 3D frame elements.

### Roark Ring Formula — Circular Ring Under Diametrically Opposite Loads
- **Source**: Roark, *Formulas for Stress and Strain*, 8th Ed., Table 9.2, Case 1
- **Problem**: Complete circular ring loaded by two equal and opposite forces P along a diameter
- **Formula**: δ_vertical = P·R³/(EI) · (π/4 - 2/π) ≈ 0.1488 · P·R³/(EI)
- **Note**: Ring must be modeled as full circle (not half/quarter). Two curved beams or one 360° arc.

### Semicircular Arch Under Crown Load (Pinned-Pinned)
- **Source**: Timoshenko & Goodier, *Theory of Elasticity*; Roark Table 9.2
- **Problem**: Semicircular arch, pinned at both ends, vertical point load at crown
- **Horizontal thrust**: H = P/π
- **Crown vertical deflection**: δ_v = P·R³/(EI) · (π/4 - 2/π) (for pin-pin semicircular arch)

---

## Plates / Shells

### Navier Solution — Simply Supported Square Plate Under Uniform Pressure
- **Source**: Timoshenko & Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd Ed., Table 8
- **Also**: Cook et al., *Concepts and Applications of FEA*, 4th Ed.
- **Problem**: Square plate, side a, thickness t, all edges simply supported, uniform pressure p
- **Formula**: w_max = α · p · a⁴ / D, where D = E·t³/(12(1-ν²))
- **Key Coefficient**: α = 0.00406 (for ν = 0.3, square plate, Table 8)
- **Also known as**: α = 0.00416 for ν = 0.0 (important to match ν!)
- **Convergence**: Need 8×8 mesh or finer for <5% error with DKT elements
- **Reference Answer Example**: a=1m, t=0.01m, E=200GPa, ν=0.3, p=1kPa → D=18315.02, w=0.00406·1·1/(18315.02)=2.216×10⁻⁷ m

### Scordelis-Lo Roof (NAFEMS)
- **Source**: NAFEMS standard test, MacNeal & Harder (1985)
- **Problem**: Cylindrical roof shell, gravity loading
- **Note**: Requires shell elements (not just flat plates). Future benchmark.

---

## Nonlinear Material / Plastic Collapse

### Plastic Collapse Load Factors (Textbook)
- **Source**: Neal, *Plastic Methods of Structural Analysis*; Chen & Sohal, *Plastic Design and Second-Order Analysis*
- **Fixed-Fixed Beam, Central Point Load**: Pc = 8·Mp/L (first hinge at midspan, then ends)
- **Fixed-Fixed Beam, UDL**: wc = 16·Mp/L² (hinges at ends then midspan)
- **Propped Cantilever, Central Point Load**: Pc = 6·Mp/L
- **Simply Supported Beam**: Pc = 4·Mp/L (single hinge at midspan)

### VM15-style — Fixed-Fixed Beam Plastic Collapse
- **Note**: ANSYS VM15 is actually a different problem (circular plate). We use the textbook plastic collapse benchmark instead.
- **Problem**: Fixed-fixed beam, length L, central point load P
- **Parameters**: L=4m, E=200GPa, b=0.15m, h=0.30m → A=0.045m², Iz=3.375×10⁻⁴ m⁴, Zp=b·h²/4=3.375×10⁻³ m³, Fy=250MPa → Mp=Fy·Zp=843.75 kN·m
- **Collapse Load**: Pc = 8·Mp/L = 8×843.75/4 = 1687.5 kN
- **Tolerance**: <5% on collapse load factor

---

## Co-rotational / Large Displacement

### ANSYS VM14 — Eccentric Compression of a Column (Large Deflection)
- **Source**: ANSYS Verification Manual, VM14
- **URL**: https://examples.mapdl.docs.pyansys.com/verif-manual/vm-014.html
- **Problem**: Simply supported column, eccentric axial load P, length L
- **Parameters**: L=120 in, rectangular 3×5 in cross section, E=30×10⁶ psi, P=4000 lb, eccentricity e=0.3 in
- **Reference Answer**: δ_mid = 0.1086 in (midspan lateral deflection, secant formula)
- **Analytical**: δ = e · [sec(L/2 · √(P/EI)) - 1]
- **Note**: This is a moderate-displacement problem suitable for co-rotational solver

### Mattiasson Elastica — Cantilever Large Deflection
- **Source**: Mattiasson (1981), "Numerical results from large deflection beam and frame problems analysed by means of elliptic integrals", Int. J. Numer. Meth. Engng., 17, 145-153
- **Problem**: Cantilever beam, length L, tip load P, dimensionless load P·L²/(EI)
- **Reference Table** (P·L²/(EI) → u_tip/L, v_tip/L):
  - 1.0 → u/L = 0.0566, v/L = 0.3015 (from Mattiasson Table 1)
  - 2.0 → u/L = 0.1632, v/L = 0.4929
  - 5.0 → u/L = 0.4054, v/L = 0.6647
  - 10.0 → u/L = 0.5584, v/L = 0.6854
- **Note**: Need to verify exact table values from the original paper

### Williams Toggle Frame
- **Source**: Williams (1964), "An approach to the non-linear behaviour of members of a rigid jointed plane framework"
- **Problem**: Symmetric two-bar toggle, apex load, snap-through behavior
- **Geometry**: Two bars meeting at shallow angle, pinned at base, load at apex
- **Critical Load**: Function of geometry and stiffness
- **Note**: Primarily tests convergence through snap-through, not a precise deflection benchmark

---

## Time Integration / Dynamics

### SDOF Free Vibration
- **Source**: Chopra, *Dynamics of Structures*, 5th Ed., Chapter 2
- **Problem**: SDOF system, free vibration after initial displacement
- **Key Result**: T = 2π√(m/k), amplitude constant (undamped) or decaying exponentially (viscous damping)

### SDOF Step Load — Dynamic Amplification Factor
- **Source**: Chopra, *Dynamics of Structures*, 5th Ed., Chapter 4
- **Problem**: Undamped SDOF under suddenly applied constant force F₀
- **Key Result**: DAF = u_max / u_static = 2.0 exactly (for undamped system)
- **Formula**: u(t) = (F₀/k)(1 - cos ωt), max when cos ωt = -1 → u_max = 2·F₀/k
- **Tolerance**: Should match 2.0 within 5%

### Newmark Energy Conservation
- **Source**: Newmark (1959), "A method of computation for structural dynamics"
- **Key Property**: Average acceleration (β=1/4, γ=1/2) is unconditionally stable and conserves energy for linear undamped systems
- **Test**: Run undamped free vibration for many periods, verify amplitude doesn't grow

### HHT-alpha Numerical Dissipation
- **Source**: Hilber, Hughes & Taylor (1977), "Improved numerical dissipation for time integration algorithms in structural dynamics"
- **Key Property**: α ∈ [-1/3, 0], introduces high-frequency damping without affecting low-frequency response
- **Test**: Compare late-time vs early-time amplitudes; HHT should show decay in high-frequency content

### Bathe SDOF Benchmarks
- **Source**: Bathe, *Finite Element Procedures*, 2nd Ed., Chapter 9
- **Problem**: SDOF spring-mass under step/impulse/sinusoidal loading
- **Note**: Exact solutions available for comparison with numerical integration

---

## Benchmark Status Summary

| Benchmark | Category | Solver Feature | Status | Notes |
|-----------|----------|----------------|--------|-------|
| Navier SS plate | Plates | DKT+CST | CAPABILITY → need 8×8 mesh, <5% | α=0.00406 |
| Roark ring | Curved beams | 3D curved expansion | NEW | Full ring needed |
| VM18 quarter-circle | Curved beams | 3D curved expansion | CAPABILITY → tighten | R=100in, δ=-2.648 |
| Plastic collapse 8Mp/L | Nonlinear material | Bilinear N-R | CAPABILITY → tighten | Fixed-fixed, central P |
| VM14 eccentric column | Co-rotational | Large displacement | NEW | δ=0.1086in |
| Mattiasson elastica | Co-rotational | Large displacement | NEW | Need paper values |
| SDOF DAF=2.0 | Time history | Newmark | CAPABILITY → tighten | Chopra Ch.4 |
| Newmark energy | Time history | Newmark | CAPABILITY → tighten | β=1/4, γ=1/2 |
| NAFEMS LE5 Z-section | Warping torsion | 7th DOF | BLOCKED | Assembly not wired |

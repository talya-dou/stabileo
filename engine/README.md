# Dedaliano Engine

High-performance 2D/3D structural analysis engine in Rust, implementing the Direct Stiffness Method from scratch with no external linear algebra dependencies.

## Analysis Types

- **Linear static** (2D & 3D): direct stiffness method, sparse Cholesky solver
- **P-Delta** (2D): second-order geometric nonlinearity with iterative convergence
- **Buckling** (2D): linearized eigenvalue buckling (Lanczos eigensolver)
- **Modal** (2D): natural frequencies and mode shapes via consistent mass matrix
- **Spectral** (2D): response spectrum analysis with SRSS/CQC combination
- **Moving loads** (2D): load envelope by stepping axle groups across the structure
- **Influence lines** (2D): Muller-Breslau virtual unit load method
- **Plastic collapse** (2D): incremental hinge formation to mechanism
- **Kinematic analysis**: mechanism detection, degree of indeterminacy, rank check

## Running Tests

```bash
cd engine && cargo test              # full suite (466 tests)
cd engine && cargo test validation_  # validation tests only (226 tests)
cd engine && cargo test diff_fuzz    # differential fuzz tests (90 tests)
```

## Validation Test Suite

226 validation tests across 25 files, verified against published analytical solutions, industry codes, and commercial software results. See [`tests/BENCHMARK_TRACKING.md`](tests/BENCHMARK_TRACKING.md) for detailed status of each benchmark.

### Tier 1: Industry Standards

| Source | File | Tests | What it covers |
|--------|------|-------|----------------|
| AISC 360-22 Ch. C | `validation_aisc_stability.rs` | 6 | B1/B2 amplification factors, P-delta convergence, equilibrium |
| Eurocode 3 EN 1993-1-1 | `validation_eurocode3_buckling.rs` | 6 | Elastic critical load factor alpha_cr, Horne's method, multi-story sway |

### Tier 2: Textbook Classics

| Source | File | Tests | What it covers |
|--------|------|-------|----------------|
| Timoshenko, Gere | `validation_beam_formulas.rs` | 14 | Euler-Bernoulli exact: SS, cantilever, fixed-fixed, propped cantilever |
| Timoshenko & Gere | `validation_pdelta_stability.rs` | 8 | Euler column buckling: 4 BCs x 4 mesh densities, P_cr convergence |
| Blevins | `validation_modal_frequencies.rs` | 16 | Beam natural frequencies: SS, cantilever, fixed-fixed x 4 modes |
| Przemieniecki | `validation_3d_analysis.rs` | 10 | 3D cantilever biaxial bending, pure torsion, axial, weak-axis, equilibrium |
| Przemieniecki | `validation_przemieniecki_extended.rs` | 6 | 12x12 stiffness symmetry, patch test, coordinate transform, hinge |
| Ziemian & Ziemian (2021) | `validation_mastan2_frames.rs` | 15 | 5 benchmark frames: alpha_cr, P-delta drift, B2 consistency, equilibrium |
| Kassimali | `validation_continuous_beams.rs` | 6 | 2-span/3-span continuous beams, settlement |
| Kassimali | `validation_frames.rs` | 7 | Portal frames, multi-story frames |
| Kassimali | `validation_kassimali_extended.rs` | 6 | Propped cantilever, partial UDL, 3-story 2-bay, truss joints, influence line |
| Biggs, Chopra | `validation_spectral_response.rs` | 8 | Spectral base shear, SRSS vs CQC, importance/reduction factors |
| Biggs, Chopra, EC8 | `validation_biggs_extended.rs` | 4 | Triangular spectrum, EC8 Type 1, shear building forces, overturning moment |

### Tier 3: Cross-Validation with Commercial Software

| Source | File | Tests | What it covers |
|--------|------|-------|----------------|
| ANSYS VM | `validation_ansys_vm.rs` | 7 | VM1 (3-bar truss), VM2 (overhangs), VM4 (V-truss), VM10 (eccentric load) |
| ANSYS VM | `validation_ansys_vm_extended.rs` | 18 | VM3/5/6/7/8/9/12/13/14/21/156: stepped beam, thermal, torsion, space truss |
| SAP2000/CSI | `validation_sap2000.rs` | 10 | Beam, continuous, portal, modal, leaning column, springs, P-delta |
| Code_Aster SSLL | `validation_code_aster.rs` | 9 | SSLL010/012/014/100/102/103/105/110/400: trusses, frames, buckling |

### Tier 4: Research Datasets

| Source | File | Tests | What it covers |
|--------|------|-------|----------------|
| Zubydan / Ziemian 22 | `validation_mastan2_frames.rs` | 15 | Same dataset as Tier 2 — 5 frames with alpha_cr + P-delta cross-validation |

### Additional Validation Files

| File | Tests | What it covers |
|------|-------|----------------|
| `validation_thermal_settlement.rs` | 10 | Thermal expansion, prescribed displacements |
| `validation_influence_lines.rs` | 8 | Muller-Breslau influence lines for beams and frames |
| `validation_section_stress.rs` | 8 | Navier bending, Jourawski shear, section properties |
| `validation_plastic_collapse.rs` | 8 | Incremental plastic hinge formation, collapse load factors |
| `validation_kinematic.rs` | 6 | Mechanism detection, isostatic/hyperstatic classification |
| `validation_moving_loads.rs` | 8 | Single axle, HL-93 truck, continuous beam envelopes |
| `validation_trusses.rs` | 6 | Planar trusses, method of joints verification |

### Blocked Benchmarks (Need New Solver Features)

| Benchmark | Blocked by |
|-----------|------------|
| NAFEMS LE5 (Z-section) | Warping torsion (Cw) |
| NAFEMS R0024 (3D beam) | Co-rotational large displacement |
| ANSYS VM11 (plate bending) | Plate/shell elements |
| ANSYS VM15 (nonlinear) | Elastoplastic material model |
| ANSYS VM18/VM44 (curved) | Curved beam elements |
| Code_Aster SSLL112 (arch) | Curved beam elements |
| DYNA3D benchmarks | Newmark/HHT time integration |

## Differential Fuzz Tests

90 tests comparing the Rust engine output against the TypeScript reference solver across random seeds, validating:
- Displacements, reactions, element forces
- Internal force diagrams at 9 interior points per element
- Support for all element types, load patterns, and boundary conditions

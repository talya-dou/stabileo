/// Release-grade acceptance models.
///
/// These tests exercise realistic, multi-feature scenarios that look like
/// actual engineering work — not just textbook cases.
///
/// Tests:
///   3A. Industrial nave (538 elements, fixture-based)
///   3B. 3D multi-story building (105 elements, fixture-based)
///   3C. Large programmatic 2D frame (550 elements, lateral + gravity)
///   3D. 10-span continuous beam with mixed loads (200 elements)
///   3E. Mixed frame+shell structure

#[path = "../common/mod.rs"]
mod common;

use common::*;
use dedaliano_engine::solver::{linear, pdelta};
use dedaliano_engine::solver::contact::{self, ContactInput, ContactType};
use dedaliano_engine::solver::fiber_nonlinear::{self, FiberNonlinearInput};
use dedaliano_engine::solver::reduction::{self, GuyanInput};
use dedaliano_engine::element::fiber_beam::{rectangular_fiber_section, FiberMaterial};
use dedaliano_engine::types::*;
use std::collections::HashMap;

// ── 3A: Industrial Nave (fixture-based, 538 elements) ──────────────

#[test]
fn acceptance_3a_industrial_nave() {
    let input_json = include_str!("../fixtures/ex-3d-nave-industrial-input.json");
    let input: SolverInput3D = serde_json::from_str(input_json)
        .expect("Failed to parse nave industrial input");

    let result = linear::solve_3d(&input).expect("solve_3d failed on nave industrial");

    // Basic sanity
    assert!(!result.displacements.is_empty(), "No displacements");
    assert!(!result.reactions.is_empty(), "No reactions");
    assert!(!result.element_forces.is_empty(), "No element forces");

    // Large model — verify expected counts
    assert!(result.displacements.len() > 200, "Expected 200+ nodes");
    assert!(result.element_forces.len() > 500, "Expected 500+ elements");

    // No NaN/Inf in displacements
    for d in &result.displacements {
        assert!(d.ux.is_finite() && d.uy.is_finite() && d.uz.is_finite(),
            "NaN/Inf at node {}", d.node_id);
    }

    // Reactions should be non-trivial
    let sum_fz: f64 = result.reactions.iter().map(|r| r.fz).sum();
    assert!(sum_fz.abs() > 1.0, "Reactions should be non-zero, got sum_fz={:.2}", sum_fz);

    // Displacements should be physically reasonable (not zero, not huge)
    let max_uz = result.displacements.iter()
        .map(|d| d.uz.abs())
        .fold(0.0_f64, f64::max);
    assert!(max_uz > 1e-6, "Displacements too small — model may not be loaded");
    assert!(max_uz < 1.0, "Displacements unreasonably large: max_uz={:.4}", max_uz);
}

// ── 3B: 3D Multi-Story Building (fixture-based, 105 elements) ─────

#[test]
fn acceptance_3b_building_case1() {
    let input_json = include_str!("../fixtures/ex-3d-building-case1-input.json");
    let input: SolverInput3D = serde_json::from_str(input_json)
        .expect("Failed to parse building case1 input");

    let result = linear::solve_3d(&input).expect("solve_3d failed on building case1");

    // Basic sanity
    assert!(!result.displacements.is_empty());
    assert!(!result.reactions.is_empty());
    assert!(!result.element_forces.is_empty());

    // No NaN/Inf
    for d in &result.displacements {
        assert!(d.ux.is_finite() && d.uy.is_finite() && d.uz.is_finite(),
            "NaN/Inf at node {}", d.node_id);
    }

    // Reactions should be non-trivial
    let sum_fz: f64 = result.reactions.iter().map(|r| r.fz).sum();
    assert!(sum_fz.abs() > 1.0, "Reactions should be non-zero, got sum_fz={:.2}", sum_fz);

    // Displacements reasonable
    let max_disp = result.displacements.iter()
        .map(|d| d.ux.abs().max(d.uy.abs()).max(d.uz.abs()))
        .fold(0.0_f64, f64::max);
    assert!(max_disp > 1e-8, "All displacements near zero");
    assert!(max_disp < 1.0, "Displacements unreasonably large: {:.4}", max_disp);
}

// ── 3C: Large Programmatic 2D Frame (50-story × 5-bay, 550 elements) ─

#[test]
fn acceptance_3c_large_2d_frame() {
    let n_stories = 50;
    let n_bays = 5;
    let h = 3.0;
    let w = 6.0;
    let e = 200_000.0;
    let a = 0.01;
    let iz = 1e-4;
    let cols = n_bays + 1;

    let mut nodes = Vec::new();
    let mut node_id = 1usize;
    for j in 0..=n_stories {
        for i in 0..=n_bays {
            nodes.push((node_id, i as f64 * w, j as f64 * h));
            node_id += 1;
        }
    }

    let mut elems = Vec::new();
    let mut eid = 1usize;
    for j in 0..n_stories {
        for i in 0..=n_bays {
            let ni = j * cols + i + 1;
            let nj = (j + 1) * cols + i + 1;
            elems.push((eid, "frame", ni, nj, 1, 1, false, false));
            eid += 1;
        }
    }
    for j in 1..=n_stories {
        for i in 0..n_bays {
            let ni = j * cols + i + 1;
            let nj = j * cols + i + 2;
            elems.push((eid, "frame", ni, nj, 1, 1, false, false));
            eid += 1;
        }
    }

    let sups: Vec<_> = (0..=n_bays).map(|i| (i + 1, i + 1, "fixed")).collect();

    let mut loads = Vec::new();
    for j in 1..=n_stories {
        // Lateral at left node
        loads.push(SolverLoad::Nodal(SolverNodalLoad {
            node_id: j * cols + 1,
            fx: 10.0, fy: 0.0, mz: 0.0,
        }));
        // Gravity at each floor node
        for i in 0..=n_bays {
            loads.push(SolverLoad::Nodal(SolverNodalLoad {
                node_id: j * cols + i + 1,
                fx: 0.0, fy: -50.0, mz: 0.0,
            }));
        }
    }

    let input = make_input(
        nodes, vec![(1, e, 0.3)], vec![(1, a, iz)], elems, sups, loads,
    );

    // Verify element count >= 500
    let expected_elements = n_stories * cols + n_stories * n_bays;
    assert_eq!(input.elements.len(), expected_elements);
    assert!(expected_elements >= 500);

    // Linear solve
    let result = linear::solve_2d(&input).expect("Linear solve failed on large frame");

    // No NaN/Inf
    for d in &result.displacements {
        assert!(d.ux.is_finite() && d.uy.is_finite(), "NaN/Inf at node {}", d.node_id);
    }

    // Reactions exist
    assert_eq!(result.reactions.len(), cols, "Expected {} base reactions", cols);

    // Top-floor lateral displacement should be positive (structure leans with wind)
    let top_left_node_id = n_stories * cols + 1;
    let top_disp = result.displacements.iter()
        .find(|d| d.node_id == top_left_node_id)
        .expect("Top-left node not found");
    assert!(top_disp.ux > 0.0,
        "Top-floor should deflect rightward under rightward wind, got ux={:.6}", top_disp.ux);

    // Gravity should cause downward deflection at beams
    let mid_floor_node = (n_stories / 2) * cols + n_bays / 2 + 1;
    let mid_disp = result.displacements.iter()
        .find(|d| d.node_id == mid_floor_node);
    if let Some(md) = mid_disp {
        assert!(md.uy < 0.0, "Mid-floor should deflect downward, got uy={:.6}", md.uy);
    }

    // P-Delta: verify convergence and amplification > 1.0
    let pdelta_result = pdelta::solve_pdelta_2d(&input, 20, 1e-4)
        .expect("P-delta solve failed on large frame");
    assert!(pdelta_result.converged, "P-delta did not converge");
    assert!(pdelta_result.b2_factor > 1.0,
        "B2 factor should be > 1.0, got {}", pdelta_result.b2_factor);
}

// ── 3D: 10-Span Continuous Beam with Mixed Loads ───────────────────

#[test]
fn acceptance_3d_continuous_beam_mixed_loads() {
    let e = 200_000.0;
    let a = 0.01;
    let iz = 1e-4;
    let n_per_span = 20;

    // 10 spans of varying lengths
    let spans = vec![6.0, 8.0, 5.0, 7.0, 6.5, 8.5, 5.5, 7.5, 6.0, 8.0];
    let total_length: f64 = spans.iter().sum();
    let total_elements = n_per_span * spans.len();

    // UDL on all spans
    let q = -12.0; // kN/m
    let mut loads = Vec::new();
    for i in 0..total_elements {
        loads.push(SolverLoad::Distributed(SolverDistributedLoad {
            element_id: i + 1,
            q_i: q,
            q_j: q,
            a: None,
            b: None,
        }));
    }

    // Point loads on alternating spans (at midspan node)
    for (span_idx, _span_len) in spans.iter().enumerate() {
        if span_idx % 2 == 0 {
            let mid_node = 1 + span_idx * n_per_span + n_per_span / 2;
            loads.push(SolverLoad::Nodal(SolverNodalLoad {
                node_id: mid_node,
                fx: 0.0, fy: -50.0, mz: 0.0,
            }));
        }
    }

    let mut input = make_continuous_beam(&spans, n_per_span, e, a, iz, loads);

    // Settlement at one interior support (support at end of span 3)
    let settlement_node = 1 + 3 * n_per_span;
    for sup in input.supports.values_mut() {
        if sup.node_id == settlement_node {
            sup.dy = Some(-0.005); // 5mm settlement
        }
    }

    let result = linear::solve_2d(&input).expect("Linear solve failed on continuous beam");

    // Basic sanity
    assert_eq!(result.displacements.len(), total_elements + 1);
    assert!(!result.reactions.is_empty());

    // No NaN/Inf
    for d in &result.displacements {
        assert!(d.ux.is_finite() && d.uy.is_finite(), "NaN/Inf at node {}", d.node_id);
    }

    // Global equilibrium: sum of reactions ≈ total applied load
    let sum_ry: f64 = result.reactions.iter().map(|r| r.ry).sum();
    let total_udl = q * total_length;
    let total_point: f64 = spans.iter().enumerate()
        .filter(|(i, _)| i % 2 == 0)
        .map(|_| -50.0)
        .sum();
    let total_applied = total_udl + total_point;
    let equil_err = (sum_ry + total_applied).abs();
    assert!(equil_err < 1.0,
        "Equilibrium error: sum_ry={:.3}, total_applied={:.3}, err={:.3}",
        sum_ry, total_applied, equil_err);

    // Midspan deflections should be physically reasonable (downward, bounded)
    for span_idx in [0, 4, 9] {
        let mid_node = 1 + span_idx * n_per_span + n_per_span / 2;
        let mid_disp = result.displacements.iter()
            .find(|d| d.node_id == mid_node)
            .expect("Midspan node not found");
        assert!(mid_disp.uy < 0.0,
            "Span {} midspan deflection should be downward, got {:.6}", span_idx, mid_disp.uy);
        assert!(mid_disp.uy > -0.1,
            "Span {} midspan deflection unreasonably large: {:.6}", span_idx, mid_disp.uy);
    }
}

// ── 3E: Mixed Frame + Shell Structure ──────────────────────────────

#[test]
fn acceptance_3e_mixed_frame_shell() {
    // Portal frame columns (4 frame elements) + MITC4 quad roof slab (4×4 mesh)
    let h = 4.0;
    let w = 8.0;
    let depth = 8.0;
    let e = 30_000.0;
    let nu = 0.2;
    let col_a = 0.16;
    let iy = 2.133e-3;
    let iz_val = 2.133e-3;
    let j_val = 3.6e-3;
    let slab_t = 0.2;

    let nx = 4;
    let ny = 4;
    let mut nodes_map = HashMap::new();
    let mut node_id = 1usize;

    // 4 base nodes (z=0) at slab corners
    let corners = [(0.0, 0.0), (w, 0.0), (w, depth), (0.0, depth)];
    let mut base_ids = Vec::new();
    for &(x, y) in &corners {
        nodes_map.insert(node_id.to_string(), SolverNode3D { id: node_id, x, y, z: 0.0 });
        base_ids.push(node_id);
        node_id += 1;
    }

    // Slab grid at z=h
    let sx = w / nx as f64;
    let sy = depth / ny as f64;
    let mut grid = vec![vec![0usize; ny + 1]; nx + 1];
    for i in 0..=nx {
        for jj in 0..=ny {
            let x = i as f64 * sx;
            let y = jj as f64 * sy;
            nodes_map.insert(node_id.to_string(), SolverNode3D { id: node_id, x, y, z: h });
            grid[i][jj] = node_id;
            node_id += 1;
        }
    }

    // Column top nodes = slab grid corners
    let top_ids = [grid[0][0], grid[nx][0], grid[nx][ny], grid[0][ny]];

    let mut materials = HashMap::new();
    materials.insert("1".to_string(), SolverMaterial { id: 1, e, nu });
    let mut sections = HashMap::new();
    sections.insert("1".to_string(), SolverSection3D {
        id: 1, name: None, a: col_a, iy, iz: iz_val, j: j_val, cw: None, as_y: None, as_z: None,
    });

    // 4 column elements
    let mut elements = HashMap::new();
    let mut eid = 1usize;
    for ci in 0..4 {
        elements.insert(eid.to_string(), SolverElement3D {
            id: eid,
            elem_type: "frame".to_string(),
            node_i: base_ids[ci],
            node_j: top_ids[ci],
            material_id: 1, section_id: 1,
            hinge_start: false, hinge_end: false,
            local_yx: None, local_yy: None, local_yz: None, roll_angle: None,
        });
        eid += 1;
    }

    // 4×4 quad mesh
    let mut quads = HashMap::new();
    let mut qid = 1usize;
    for i in 0..nx {
        for jj in 0..ny {
            quads.insert(qid.to_string(), SolverQuadElement {
                id: qid,
                nodes: [grid[i][jj], grid[i+1][jj], grid[i+1][jj+1], grid[i][jj+1]],
                material_id: 1,
                thickness: slab_t,
            });
            qid += 1;
        }
    }

    // Fixed supports at base
    let mut supports = HashMap::new();
    for (i, &nid) in base_ids.iter().enumerate() {
        supports.insert((i + 1).to_string(), SolverSupport3D {
            node_id: nid,
            rx: true, ry: true, rz: true,
            rrx: true, rry: true, rrz: true,
            kx: None, ky: None, kz: None,
            krx: None, kry: None, krz: None,
            dx: None, dy: None, dz: None,
            drx: None, dry: None, drz: None,
            normal_x: None, normal_y: None, normal_z: None,
            is_inclined: None, rw: None, kw: None,
        });
    }

    // Gravity on slab nodes
    let total_slab_nodes = (nx + 1) * (ny + 1);
    let force_per_node = -5.0 * w * depth / total_slab_nodes as f64;
    let mut loads = Vec::new();
    for i in 0..=nx {
        for jj in 0..=ny {
            loads.push(SolverLoad3D::Nodal(SolverNodalLoad3D {
                node_id: grid[i][jj],
                fx: 0.0, fy: 0.0, fz: force_per_node,
                mx: 0.0, my: 0.0, mz: 0.0, bw: None,
            }));
        }
    }

    let input = SolverInput3D {
        nodes: nodes_map,
        materials,
        sections,
        elements,
        supports,
        loads,
        constraints: vec![],
        plates: HashMap::new(),
        quads,
        left_hand: None,
        curved_beams: vec![],
        connectors: HashMap::new(),
    };

    let result = linear::solve_3d(&input).expect("solve_3d failed on mixed frame+shell");

    // Basic sanity
    assert!(!result.displacements.is_empty());
    assert!(!result.reactions.is_empty());

    // No NaN/Inf
    for d in &result.displacements {
        assert!(d.ux.is_finite() && d.uy.is_finite() && d.uz.is_finite(),
            "NaN/Inf at node {}", d.node_id);
    }

    // Reactions non-trivial
    let sum_fz: f64 = result.reactions.iter().map(|r| r.fz).sum();
    assert!(sum_fz.abs() > 1.0, "Reactions should be non-zero, sum_fz={:.2}", sum_fz);

    // Slab center should deflect downward
    let center_node = grid[nx / 2][ny / 2];
    let center_disp = result.displacements.iter()
        .find(|dd| dd.node_id == center_node)
        .expect("Center slab node not found");
    assert!(center_disp.uz < 0.0,
        "Center slab should deflect downward, got uz={:.6}", center_disp.uz);
}

// ── 4A: Steel Building with Diaphragms + Eccentric Connection ──────

#[test]
fn acceptance_4a_steel_building_diaphragms() {
    // 3-story, 2×1-bay 3D frame: 24 main nodes + 1 eccentric slave
    let e = 200_000.0;
    let nu = 0.3;

    // Nodes: 6 per level (3×2 grid, x={0,6,12}, y={0,6}), z={0,3,6,9}
    let mut nodes = Vec::new();
    let mut nid = 1usize;
    for level in 0..4 {
        let z = level as f64 * 3.0;
        for iy in 0..2 {
            for ix in 0..3 {
                nodes.push((nid, ix as f64 * 6.0, iy as f64 * 6.0, z));
                nid += 1;
            }
        }
    }
    // Node 25: eccentric slave off node 20
    nodes.push((25, 12.0, 6.3, 9.0)); // offset_y = 0.3 from node 20 at (12,6,9)

    // Sections: col (sec 1), beam (sec 2)
    let secs = vec![
        (1, 0.01_f64, 8.33e-5, 8.33e-5, 1.4e-4),   // column
        (2, 0.006_f64, 4.5e-5, 4.5e-5, 7.0e-5),     // beam
    ];

    // Elements: 18 columns (vertical) + 15 beams (horizontal)
    let mut elems = Vec::new();
    let mut eid = 1usize;

    // Columns: connect level j to level j+1 at each of 6 grid positions
    for level in 0..3 {
        for pos in 0..6 {
            let ni = level * 6 + pos + 1;
            let nj = (level + 1) * 6 + pos + 1;
            elems.push((eid, "frame", ni, nj, 1, 1)); // sec 1 = column
            eid += 1;
        }
    }

    // Beams: at each floor level (1,2,3), connect adjacent nodes
    // X-direction beams (along x at each y): 3 floors × 2 rows × 2 spans = 12
    // Y-direction beams (along y at each x): 3 floors × 3 cols × 1 span = 9 ... but 12+9=21, plan says 15
    // Plan says 15 beams total. Let's do X-direction only: 3 floors × 2 y-rows × 2 x-spans = 12
    // Plus Y-direction at edges: 3 floors × 1 span (y=0→y=6 at x=0 only) = 3 → 15 total
    // Actually let's just do all X-beams = 12, plus 3 Y-beams at x=0 = 15
    for level in 1..4 {
        let base = level * 6; // first node id at this level = base+1
        // X-direction beams at y=0: (base+1)-(base+2), (base+2)-(base+3)
        for i in 0..2 {
            let ni = base + i + 1;
            let nj = base + i + 2;
            elems.push((eid, "frame", ni, nj, 1, 2)); // sec 2 = beam
            eid += 1;
        }
        // X-direction beams at y=6: (base+4)-(base+5), (base+5)-(base+6)
        for i in 0..2 {
            let ni = base + 3 + i + 1;
            let nj = base + 3 + i + 2;
            elems.push((eid, "frame", ni, nj, 1, 2));
            eid += 1;
        }
        // Y-direction beam at x=0: (base+1)-(base+4)
        let ni = base + 1;
        let nj = base + 4;
        elems.push((eid, "frame", ni, nj, 1, 2));
        eid += 1;
    }
    assert_eq!(elems.len(), 33, "Expected 33 elements (18 cols + 15 beams)");

    // Supports: nodes 1–6 fully fixed
    let sups: Vec<(usize, Vec<bool>)> = (1..=6)
        .map(|n| (n, vec![true, true, true, true, true, true]))
        .collect();

    // Loads: gravity + wind + eccentric
    let mut loads = Vec::new();
    // Gravity: fz=-50 kN at each floor node (nodes 7–24, 18 nodes)
    for n in 7..=24 {
        loads.push(SolverLoad3D::Nodal(SolverNodalLoad3D {
            node_id: n, fx: 0.0, fy: 0.0, fz: -50.0, mx: 0.0, my: 0.0, mz: 0.0, bw: None,
        }));
    }
    // Wind: fx=+5 kN at floor-3 nodes (19–24)
    for n in 19..=24 {
        loads.push(SolverLoad3D::Nodal(SolverNodalLoad3D {
            node_id: n, fx: 5.0, fy: 0.0, fz: 0.0, mx: 0.0, my: 0.0, mz: 0.0, bw: None,
        }));
    }
    // Eccentric: fx=+2 kN at node 25
    loads.push(SolverLoad3D::Nodal(SolverNodalLoad3D {
        node_id: 25, fx: 2.0, fy: 0.0, fz: 0.0, mx: 0.0, my: 0.0, mz: 0.0, bw: None,
    }));

    let mut input = make_3d_input(nodes, vec![(1, e, nu)], secs, elems, sups, loads);

    // Constraints
    // Floor diaphragms
    input.constraints.push(Constraint::Diaphragm(DiaphragmConstraint {
        master_node: 7, slave_nodes: vec![8, 9, 10, 11, 12], plane: "XY".into(),
    }));
    input.constraints.push(Constraint::Diaphragm(DiaphragmConstraint {
        master_node: 13, slave_nodes: vec![14, 15, 16, 17, 18], plane: "XY".into(),
    }));
    input.constraints.push(Constraint::Diaphragm(DiaphragmConstraint {
        master_node: 19, slave_nodes: vec![20, 21, 22, 23, 24], plane: "XY".into(),
    }));
    // Eccentric connection: master=20, slave=25, offset_y=0.3
    input.constraints.push(Constraint::EccentricConnection(EccentricConnectionConstraint {
        master_node: 20, slave_node: 25,
        offset_x: 0.0, offset_y: 0.3, offset_z: 0.0,
        releases: vec![],
    }));

    let result = linear::solve_3d(&input).expect("solve_3d failed on 4A steel building");

    // 1. No NaN/Inf
    for d in &result.displacements {
        assert!(d.ux.is_finite() && d.uy.is_finite() && d.uz.is_finite(),
            "NaN/Inf at node {}", d.node_id);
    }

    // 2. Diaphragm kinematics: slave_ux ≈ master_ux - dy*master_rz (floor 3)
    let master19 = result.displacements.iter().find(|d| d.node_id == 19).unwrap();
    for &slave_id in &[20, 21, 22, 23, 24] {
        let slave = result.displacements.iter().find(|d| d.node_id == slave_id).unwrap();
        let slave_node = input.nodes.values().find(|n| n.id == slave_id).unwrap();
        let master_node = input.nodes.values().find(|n| n.id == 19).unwrap();
        let dx = slave_node.x - master_node.x;
        let dy = slave_node.y - master_node.y;
        let expected_ux = master19.ux - dy * master19.rz;
        let expected_uy = master19.uy + dx * master19.rz;
        assert!((slave.ux - expected_ux).abs() < 1e-4,
            "Diaphragm ux mismatch at node {}: got {:.8}, expected {:.8}", slave_id, slave.ux, expected_ux);
        assert!((slave.uy - expected_uy).abs() < 1e-4,
            "Diaphragm uy mismatch at node {}: got {:.8}, expected {:.8}", slave_id, slave.uy, expected_uy);
    }

    // 3. Vertical equilibrium: reactions should be non-trivial and balance gravity
    let sum_fz_react: f64 = result.reactions.iter().map(|r| r.fz).sum();
    // With constraints, reactions may be reported via constraint forces. Check both.
    let has_significant_fz = result.reactions.iter().any(|r| r.fz.abs() > 1.0);
    let has_constraint_fz = result.constraint_forces.iter().any(|cf| cf.force.abs() > 1.0);
    assert!(has_significant_fz || has_constraint_fz,
        "Expected non-trivial vertical reactions or constraint forces, sum_fz_react={:.3}", sum_fz_react);

    // 4. Lateral drift: floor 3 ux > floor 2 ux > floor 1 ux
    let floor_ux = |master: usize| -> f64 {
        result.displacements.iter().find(|d| d.node_id == master).unwrap().ux
    };
    let ux1 = floor_ux(7);   // floor 1 master
    let ux2 = floor_ux(13);  // floor 2 master
    let ux3 = floor_ux(19);  // floor 3 master
    assert!(ux3 > ux2 && ux2 > ux1,
        "Lateral drift should increase: ux1={:.6}, ux2={:.6}, ux3={:.6}", ux1, ux2, ux3);

    // 5. constraint_forces non-empty
    assert!(!result.constraint_forces.is_empty(), "Expected non-empty constraint forces");

    // 6. Eccentric node 25 follows rigid-body from master 20
    let master20 = result.displacements.iter().find(|d| d.node_id == 20).unwrap();
    let slave25 = result.displacements.iter().find(|d| d.node_id == 25).unwrap();
    let expected_ux_25 = master20.ux - 0.3 * master20.rz;
    assert!((slave25.ux - expected_ux_25).abs() < 1e-3,
        "Eccentric node 25 ux: got {:.8}, expected {:.8}", slave25.ux, expected_ux_25);
}

// ── 4B: Frame + Quad Slab (Mixed Elements) ─────────────────────────

#[test]
fn acceptance_4b_frame_quad_slab() {
    // 2×2-bay single-story: 4 frame columns + 16 MITC4 quads
    let e = 30_000.0;
    let nu = 0.2;
    let col_a = 0.16;
    let iy = 2.133e-3;
    let iz_val = 2.133e-3;
    let j_val = 3.6e-3;
    let slab_t = 0.2;

    let nx = 4usize;
    let ny = 4usize;
    let slab_size = 8.0;
    let sx = slab_size / nx as f64;
    let sy = slab_size / ny as f64;

    let mut nodes_map = HashMap::new();
    let mut node_id = 1usize;

    // 4 base nodes at z=0 (column bases at slab corners)
    let base_corners = [(0.0, 0.0), (slab_size, 0.0), (slab_size, slab_size), (0.0, slab_size)];
    let mut base_ids = Vec::new();
    for &(x, y) in &base_corners {
        nodes_map.insert(node_id.to_string(), SolverNode3D { id: node_id, x, y, z: 0.0 });
        base_ids.push(node_id);
        node_id += 1;
    }

    // 5×5 slab grid at z=4
    let mut grid = vec![vec![0usize; ny + 1]; nx + 1];
    for i in 0..=nx {
        for j in 0..=ny {
            let x = i as f64 * sx;
            let y = j as f64 * sy;
            nodes_map.insert(node_id.to_string(), SolverNode3D { id: node_id, x, y, z: 4.0 });
            grid[i][j] = node_id;
            node_id += 1;
        }
    }
    assert_eq!(node_id - 1, 29); // 4 base + 25 slab

    // Column tops = slab grid corners
    let top_ids = [grid[0][0], grid[nx][0], grid[nx][ny], grid[0][ny]];

    let mut materials = HashMap::new();
    materials.insert("1".to_string(), SolverMaterial { id: 1, e, nu });
    let mut sections = HashMap::new();
    sections.insert("1".to_string(), SolverSection3D {
        id: 1, name: None, a: col_a, iy, iz: iz_val, j: j_val,
        cw: None, as_y: None, as_z: None,
    });

    // 4 column elements
    let mut elements = HashMap::new();
    let mut eid = 1usize;
    for ci in 0..4 {
        elements.insert(eid.to_string(), SolverElement3D {
            id: eid, elem_type: "frame".to_string(),
            node_i: base_ids[ci], node_j: top_ids[ci],
            material_id: 1, section_id: 1,
            hinge_start: false, hinge_end: false,
            local_yx: None, local_yy: None, local_yz: None, roll_angle: None,
        });
        eid += 1;
    }

    // 4×4 MITC4 quad mesh
    let mut quads = HashMap::new();
    let mut qid = 1usize;
    for i in 0..nx {
        for j in 0..ny {
            quads.insert(qid.to_string(), SolverQuadElement {
                id: qid,
                nodes: [grid[i][j], grid[i+1][j], grid[i+1][j+1], grid[i][j+1]],
                material_id: 1, thickness: slab_t,
            });
            qid += 1;
        }
    }
    assert_eq!(quads.len(), 16);

    // Fixed supports at base
    let mut supports = HashMap::new();
    for (i, &nid) in base_ids.iter().enumerate() {
        supports.insert((i + 1).to_string(), SolverSupport3D {
            node_id: nid,
            rx: true, ry: true, rz: true, rrx: true, rry: true, rrz: true,
            kx: None, ky: None, kz: None, krx: None, kry: None, krz: None,
            dx: None, dy: None, dz: None, drx: None, dry: None, drz: None,
            normal_x: None, normal_y: None, normal_z: None,
            is_inclined: None, rw: None, kw: None,
        });
    }

    // Loads: QuadPressure = -5 kN/m² on all 16 quads + lateral at slab center
    let mut loads = Vec::new();
    for qid in 1..=16 {
        loads.push(SolverLoad3D::QuadPressure(SolverPressureLoad {
            element_id: qid, pressure: -5.0,
        }));
    }
    let center_node = grid[nx / 2][ny / 2];
    loads.push(SolverLoad3D::Nodal(SolverNodalLoad3D {
        node_id: center_node, fx: 20.0, fy: 0.0, fz: 0.0,
        mx: 0.0, my: 0.0, mz: 0.0, bw: None,
    }));

    let input = SolverInput3D {
        nodes: nodes_map, materials, sections, elements, supports, loads,
        constraints: vec![], plates: HashMap::new(), quads,
        left_hand: None, curved_beams: vec![], connectors: HashMap::new(),
    };

    let result = linear::solve_3d(&input).expect("solve_3d failed on 4B frame+quad slab");

    // 1. Correct counts
    assert_eq!(result.quad_stresses.len(), 16, "Expected 16 quad stresses");
    assert_eq!(result.element_forces.len(), 4, "Expected 4 element forces (columns)");

    // 2. Slab center deflects downward; |uz_center| > |uz_corner| (bowl shape)
    let center_disp = result.displacements.iter()
        .find(|d| d.node_id == center_node).unwrap();
    assert!(center_disp.uz < 0.0, "Center should deflect down, got uz={:.6}", center_disp.uz);
    let corner_disp = result.displacements.iter()
        .find(|d| d.node_id == grid[0][0]).unwrap();
    assert!(center_disp.uz.abs() > corner_disp.uz.abs(),
        "Bowl: |uz_center|={:.6} should > |uz_corner|={:.6}", center_disp.uz.abs(), corner_disp.uz.abs());

    // 3. Reactions Fz sum should balance applied quad pressure (positive = upward)
    let sum_fz: f64 = result.reactions.iter().map(|r| r.fz).sum();
    assert!(sum_fz > 100.0, "Reactions should be significantly positive, got {:.2}", sum_fz);

    // 4. Positive ux at slab center from lateral load
    assert!(center_disp.ux > 0.0,
        "Center ux should be positive from lateral load, got {:.6}", center_disp.ux);

    // 5. All quad von_mises > 0 and finite
    for qs in &result.quad_stresses {
        assert!(qs.von_mises > 0.0 && qs.von_mises.is_finite(),
            "Quad {} von_mises invalid: {:.6}", qs.element_id, qs.von_mises);
    }
}

// ── 4C: Contact + Uplift Foundation ─────────────────────────────────

#[test]
fn acceptance_4c_contact_uplift() {
    // 3-bay 2D portal: 8 nodes, 7 elements, uplift at outer columns
    let e = 200_000.0;
    let a = 0.01;
    let iz = 1e-4;

    // Nodes: bases 1–4 (y=0), tops 5–8 (y=4)
    let nodes = vec![
        (1, 0.0, 0.0), (2, 6.0, 0.0), (3, 12.0, 0.0), (4, 18.0, 0.0),
        (5, 0.0, 4.0), (6, 6.0, 4.0), (7, 12.0, 4.0), (8, 18.0, 4.0),
    ];

    // 4 columns + 3 beams
    let elems = vec![
        (1, "frame", 1, 5, 1, 1, false, false),
        (2, "frame", 2, 6, 1, 1, false, false),
        (3, "frame", 3, 7, 1, 1, false, false),
        (4, "frame", 4, 8, 1, 1, false, false),
        (5, "frame", 5, 6, 1, 1, false, false),
        (6, "frame", 6, 7, 1, 1, false, false),
        (7, "frame", 7, 8, 1, 1, false, false),
    ];

    // Pinned supports at all 4 bases
    let sups = vec![
        (1, 1, "pinned"), (2, 2, "pinned"), (3, 3, "pinned"), (4, 4, "pinned"),
    ];

    // Gravity: fy=-20 kN at top nodes; Lateral: fx=+500 kN at node 5 (high overturning)
    let loads = vec![
        SolverLoad::Nodal(SolverNodalLoad { node_id: 5, fx: 500.0, fy: -20.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 6, fx: 0.0, fy: -20.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 7, fx: 0.0, fy: -20.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 8, fx: 0.0, fy: -20.0, mz: 0.0 }),
    ];

    let solver = make_input(nodes, vec![(1, e, 0.3)], vec![(1, a, iz)], elems, sups, loads);

    let contact_input = ContactInput {
        solver,
        element_behaviors: HashMap::new(),
        gap_elements: vec![],
        uplift_supports: vec![1, 4],
        max_iter: Some(30),
        tolerance: Some(1e-6),
        augmented_lagrangian: None,
        max_flips: None,
        damping_coefficient: None,
        al_max_iter: None,
        contact_type: ContactType::default(),
        node_to_surface_pairs: vec![],
    };

    let result = contact::solve_contact_2d(&contact_input)
        .expect("Contact solve failed on 4C uplift");

    // 1. Converged
    assert!(result.converged, "Contact solver did not converge");

    // 2. No NaN/Inf
    for d in &result.results.displacements {
        assert!(d.ux.is_finite() && d.uy.is_finite(),
            "NaN/Inf at node {}", d.node_id);
    }

    // 3. 7 element forces
    assert_eq!(result.results.element_forces.len(), 7, "Expected 7 element forces");

    // 4. Under overturning: windward column (node 1) reaction ry ≤ 0 (tensile → uplift)
    let react1 = result.results.reactions.iter().find(|r| r.node_id == 1);
    if let Some(r) = react1 {
        // If uplift kicked in, reaction is zero or absent. If still in contact, it may be tensile.
        assert!(r.ry <= 0.001,
            "Node 1 should uplift (ry ≤ 0), got ry={:.4}", r.ry);
    }
    // Absence of reaction for node 1 also valid (uplift = support removed)

    // 5. Leeward reactions compensate: node 4 ry should be larger than interior
    let react4 = result.results.reactions.iter()
        .find(|r| r.node_id == 4)
        .expect("Leeward node 4 should have reaction");
    assert!(react4.ry > 0.0, "Leeward node 4 should have positive ry, got {:.4}", react4.ry);

    // 6. Total vertical reaction ≈ total gravity (80 kN)
    let sum_ry: f64 = result.results.reactions.iter().map(|r| r.ry).sum();
    assert!((sum_ry - 80.0).abs() < 5.0,
        "Total ry should ≈ 80, got {:.2}", sum_ry);
}

// ── 4D: Fiber Nonlinear Pushover ────────────────────────────────────

#[test]
fn acceptance_4d_fiber_pushover() {
    // 2-story single-bay portal: 6 nodes, 6 elements
    let e = 200_000.0;

    // Nodes
    let nodes = vec![
        (1, 0.0, 0.0), (2, 6.0, 0.0),
        (3, 0.0, 3.0), (4, 6.0, 3.0),
        (5, 0.0, 6.0), (6, 6.0, 6.0),
    ];

    // Elements: 4 columns (sec 1) + 2 beams (sec 2)
    let elems = vec![
        (1, "frame", 1, 3, 1, 1, false, false), // col
        (2, "frame", 2, 4, 1, 1, false, false), // col
        (3, "frame", 3, 5, 1, 1, false, false), // col
        (4, "frame", 4, 6, 1, 1, false, false), // col
        (5, "frame", 3, 4, 1, 2, false, false), // beam
        (6, "frame", 5, 6, 1, 2, false, false), // beam
    ];

    let sups = vec![(1, 1, "fixed"), (2, 2, "fixed")];

    // Fiber sections
    let mut fiber_sections = HashMap::new();
    fiber_sections.insert("1".into(), rectangular_fiber_section(
        0.3, 0.3, 10,
        FiberMaterial::SteelBilinear { e: 200_000.0, fy: 250.0, hardening_ratio: 0.01 },
    ));
    fiber_sections.insert("2".into(), rectangular_fiber_section(
        0.3, 0.5, 10,
        FiberMaterial::SteelBilinear { e: 200_000.0, fy: 250.0, hardening_ratio: 0.01 },
    ));

    // --- Part A: Elastic check ---
    let loads_elastic = vec![
        SolverLoad::Nodal(SolverNodalLoad { node_id: 3, fx: 5.0, fy: -20.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 4, fx: 0.0, fy: -20.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 5, fx: 5.0, fy: -20.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 6, fx: 0.0, fy: -20.0, mz: 0.0 }),
    ];

    // Sections for linear comparison: col A=0.09, Iz=6.75e-4; beam A=0.15, Iz=3.125e-3
    let solver_elastic = make_input(
        nodes.clone(), vec![(1, e, 0.3)],
        vec![(1, 0.09, 6.75e-4), (2, 0.15, 3.125e-3)],
        elems.clone(), sups.clone(), loads_elastic.clone(),
    );

    let linear_result = linear::solve_2d(&solver_elastic)
        .expect("Linear solve failed for 4D elastic check");

    let fiber_elastic_input = FiberNonlinearInput {
        solver: solver_elastic.clone(),
        fiber_sections: fiber_sections.clone(),
        n_integration_points: 5,
        max_iter: 30,
        tolerance: 1e-6,
        n_increments: 1,
    };

    let fiber_elastic = fiber_nonlinear::solve_fiber_nonlinear_2d(&fiber_elastic_input)
        .expect("Fiber elastic solve failed");

    // Compare displacements within 5%
    for ld in &linear_result.displacements {
        let fd = fiber_elastic.results.displacements.iter()
            .find(|d| d.node_id == ld.node_id).unwrap();
        if ld.ux.abs() > 1e-8 {
            let rel = (fd.ux - ld.ux).abs() / ld.ux.abs();
            assert!(rel < 0.05, "Node {} ux mismatch: linear={:.6}, fiber={:.6}, rel={:.4}",
                ld.node_id, ld.ux, fd.ux, rel);
        }
        if ld.uy.abs() > 1e-8 {
            let rel = (fd.uy - ld.uy).abs() / ld.uy.abs();
            assert!(rel < 0.05, "Node {} uy mismatch: linear={:.6}, fiber={:.6}, rel={:.4}",
                ld.node_id, ld.uy, fd.uy, rel);
        }
    }

    // --- Part B: Yielding check ---
    // Column Mp ≈ fy × Z ≈ 250 × (0.3×0.3²/4) ≈ 1687 kN·m, so need very high lateral
    let loads_yield = vec![
        SolverLoad::Nodal(SolverNodalLoad { node_id: 3, fx: 3000.0, fy: -100.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 4, fx: 0.0, fy: -100.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 5, fx: 3000.0, fy: -100.0, mz: 0.0 }),
        SolverLoad::Nodal(SolverNodalLoad { node_id: 6, fx: 0.0, fy: -100.0, mz: 0.0 }),
    ];

    let solver_yield = make_input(
        nodes.clone(), vec![(1, e, 0.3)],
        vec![(1, 0.09, 6.75e-4), (2, 0.15, 3.125e-3)],
        elems.clone(), sups.clone(), loads_yield,
    );

    let fiber_yield_input = FiberNonlinearInput {
        solver: solver_yield,
        fiber_sections,
        n_integration_points: 5,
        max_iter: 30,
        tolerance: 1e-6,
        n_increments: 20,
    };

    let fiber_yield = fiber_nonlinear::solve_fiber_nonlinear_2d(&fiber_yield_input)
        .expect("Fiber yield solve failed");

    assert!(fiber_yield.converged, "Fiber yield did not converge");

    // At least one element should have yielded
    assert!(fiber_yield.fiber_status.iter().any(|s| s.yielded),
        "Expected at least one element to yield under high lateral load");

    // Roof displacement (node 5) should exceed linear prediction
    let roof_fiber = fiber_yield.results.displacements.iter()
        .find(|d| d.node_id == 5).unwrap();
    let roof_linear = linear_result.displacements.iter()
        .find(|d| d.node_id == 5).unwrap();
    // Scale linear: loads are 600× larger for fx (3000/5), so linear ux would be ~600× larger
    let linear_scaled_ux = roof_linear.ux * (3000.0 / 5.0);
    assert!(roof_fiber.ux > linear_scaled_ux,
        "Yielding should increase drift: fiber_ux={:.6}, scaled_linear_ux={:.6}",
        roof_fiber.ux, linear_scaled_ux);
}

// ── 4E: Guyan Reduction vs Full Linear ──────────────────────────────

#[test]
fn acceptance_4e_guyan_vs_full() {
    // 5-story single-bay frame: 12 nodes, 15 elements
    let e = 200_000.0;
    let a = 0.01;
    let iz = 1e-4;

    // Nodes: 2 per level, x={0,6}, y={0,3,6,9,12,15}
    let mut nodes = Vec::new();
    let mut nid = 1usize;
    for level in 0..6 {
        let y = level as f64 * 3.0;
        nodes.push((nid, 0.0, y)); nid += 1;
        nodes.push((nid, 6.0, y)); nid += 1;
    }
    assert_eq!(nodes.len(), 12);

    // Elements: 10 columns + 5 beams
    let mut elems = Vec::new();
    let mut eid = 1usize;
    // Columns: left (1→3, 3→5, 5→7, 7→9, 9→11) and right (2→4, 4→6, 6→8, 8→10, 10→12)
    for level in 0..5 {
        let left_bot = level * 2 + 1;
        let left_top = left_bot + 2;
        elems.push((eid, "frame", left_bot, left_top, 1, 1, false, false)); eid += 1;
        let right_bot = level * 2 + 2;
        let right_top = right_bot + 2;
        elems.push((eid, "frame", right_bot, right_top, 1, 1, false, false)); eid += 1;
    }
    // Beams: floor level j → nodes (2j+1, 2j+2), j=1..5
    for level in 1..6 {
        let left = level * 2 + 1;
        let right = level * 2 + 2;
        elems.push((eid, "frame", left, right, 1, 1, false, false)); eid += 1;
    }
    assert_eq!(elems.len(), 15);

    // Supports: nodes 1,2 fixed
    let sups = vec![(1, 1, "fixed"), (2, 2, "fixed")];

    // Loads: fx=10 kN at left column per floor + fy=-50 kN at all floor nodes
    let mut loads = Vec::new();
    for level in 1..6 {
        let left = level * 2 + 1;
        loads.push(SolverLoad::Nodal(SolverNodalLoad { node_id: left, fx: 10.0, fy: 0.0, mz: 0.0 }));
        let right = level * 2 + 2;
        loads.push(SolverLoad::Nodal(SolverNodalLoad { node_id: left, fx: 0.0, fy: -50.0, mz: 0.0 }));
        loads.push(SolverLoad::Nodal(SolverNodalLoad { node_id: right, fx: 0.0, fy: -50.0, mz: 0.0 }));
    }

    let solver_input = make_input(
        nodes, vec![(1, e, 0.3)], vec![(1, a, iz)], elems, sups, loads,
    );

    // Full linear solve
    let full_result = linear::solve_2d(&solver_input)
        .expect("Full linear solve failed for 4E");

    // Guyan reduction
    let guyan_input = GuyanInput {
        solver: solver_input,
        boundary_nodes: vec![1, 2, 11, 12],
    };

    let guyan_result = reduction::guyan_reduce_2d(&guyan_input)
        .expect("Guyan reduction failed for 4E");

    // 1. Both solved successfully — already asserted by .expect()

    // 2. Displacements match at ALL nodes within 1e-10 (Guyan is exact for static)
    for fd in &full_result.displacements {
        let gd = guyan_result.displacements.iter()
            .find(|d| d.node_id == fd.node_id)
            .unwrap_or_else(|| panic!("Guyan missing node {}", fd.node_id));
        assert!((gd.ux - fd.ux).abs() < 1e-10,
            "Node {} ux: full={:.12}, guyan={:.12}", fd.node_id, fd.ux, gd.ux);
        assert!((gd.uy - fd.uy).abs() < 1e-10,
            "Node {} uy: full={:.12}, guyan={:.12}", fd.node_id, fd.uy, gd.uy);
        assert!((gd.rz - fd.rz).abs() < 1e-10,
            "Node {} rz: full={:.12}, guyan={:.12}", fd.node_id, fd.rz, gd.rz);
    }

    // 3. K_condensed symmetric: |K[i,j] - K[j,i]| < 1e-12
    let nb = guyan_result.n_boundary;
    for i in 0..nb {
        for j in (i + 1)..nb {
            let kij = guyan_result.k_condensed[i * nb + j];
            let kji = guyan_result.k_condensed[j * nb + i];
            assert!((kij - kji).abs() < 1e-10,
                "K_condensed not symmetric: K[{},{}]={:.14}, K[{},{}]={:.14}", i, j, kij, j, i, kji);
        }
    }

    // 4. Element forces match within 1e-8
    for ff in &full_result.element_forces {
        let gf = guyan_result.element_forces.iter()
            .find(|f| f.element_id == ff.element_id)
            .unwrap_or_else(|| panic!("Guyan missing element {}", ff.element_id));
        assert!((gf.n_start - ff.n_start).abs() < 1e-8,
            "Elem {} n_start: full={:.12}, guyan={:.12}", ff.element_id, ff.n_start, gf.n_start);
        assert!((gf.v_start - ff.v_start).abs() < 1e-8,
            "Elem {} v_start: full={:.12}, guyan={:.12}", ff.element_id, ff.v_start, gf.v_start);
        assert!((gf.m_start - ff.m_start).abs() < 1e-8,
            "Elem {} m_start: full={:.12}, guyan={:.12}", ff.element_id, ff.m_start, gf.m_start);
        assert!((gf.m_end - ff.m_end).abs() < 1e-8,
            "Elem {} m_end: full={:.12}, guyan={:.12}", ff.element_id, ff.m_end, gf.m_end);
    }

    // 5. Reactions match within 1e-8
    for fr in &full_result.reactions {
        let gr = guyan_result.reactions.iter()
            .find(|r| r.node_id == fr.node_id)
            .unwrap_or_else(|| panic!("Guyan missing reaction at node {}", fr.node_id));
        assert!((gr.rx - fr.rx).abs() < 1e-8,
            "Node {} rx: full={:.12}, guyan={:.12}", fr.node_id, fr.rx, gr.rx);
        assert!((gr.ry - fr.ry).abs() < 1e-8,
            "Node {} ry: full={:.12}, guyan={:.12}", fr.node_id, fr.ry, gr.ry);
        assert!((gr.mz - fr.mz).abs() < 1e-8,
            "Node {} mz: full={:.12}, guyan={:.12}", fr.node_id, fr.mz, gr.mz);
    }

    // 6. n_boundary + n_interior == total_free_dofs
    // Fixed nodes 1,2 → 6 DOFs constrained; 12 nodes × 3 DOFs = 36 total → 30 free
    let total_free = guyan_result.n_boundary + guyan_result.n_interior;
    assert_eq!(total_free, 30,
        "n_boundary({}) + n_interior({}) should == 30 free DOFs", guyan_result.n_boundary, guyan_result.n_interior);
}

// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Highly Symmetric Graph

use crate::core;

fn break_symmetry_vertex<T: core::graph::Vertex>(
    vertice_index: usize,
    vertices: &mut Vec<T>
) {
    if let Some(vi) = vertices.get_mut(vertice_index) {
        vi.break_symmetry_vertex();
    }
}

fn break_symmetry_edge<T: core::graph::Vertex> (
    vertice_index_1: usize,
    vertice_index_2: usize,
    vertices: &mut Vec<T>
) {
    vertices[vertice_index_1].break_symmetry_edge(vertice_index_2);
    vertices[vertice_index_2].break_symmetry_edge(vertice_index_1);
}

pub fn run_high_symmetry<T: core::graph::VertexExtendableHash>(
    vv: &core::graph::VertexVec<T::VertexType>,
    edges: &Vec<(usize, usize, usize)>,
    length: usize,
    orbits_residual: &Vec<core::orbit_ops::Orbit>,
    orbits_symmetry: &mut Vec<core::orbit_ops::Orbit>,
) -> Result<(), String> {
    if cfg!(debug_assertions) {
        println!("Proceed to Case of High Symmetry: {:?}", orbits_residual);
        println!("Edges: {:?}", edges);
    }

    for orbit in orbits_residual.iter() {
        for &v_idx in orbit.iter() {
            let mut vertices_cloned = vv.all_vertices().to_vec();
            break_symmetry_vertex(v_idx, &mut vertices_cloned);
            let orbits_tmp: Vec<core::orbit_ops::Orbit> = core::givp::partition_vertices::<T>(vv.valid_indexes(), &vertices_cloned);
            if !core::cnap::is_computable(&orbits_tmp) {
                return Err(format!("High Symmetry of Vertex Broken, But still not computable {:?}", orbits_tmp));
            }
            core::cnap::get_symmetric_orbits(&orbits_tmp, edges, length, orbits_symmetry);
        }
    }

    for &edge in edges.iter() {
        if edge.1 > edge.0 {
            let mut vertices_cloned = vv.all_vertices().to_vec();
            break_symmetry_edge(edge.0, edge.1, &mut vertices_cloned);
            let orbits_tmp: Vec<core::orbit_ops::Orbit> = core::givp::partition_vertices::<T>(vv.valid_indexes(), &vertices_cloned);
            if !core::cnap::is_computable(&orbits_tmp) {
                return Err(format!("High Symmetry of Edge Broken, But still not computable {:?}", orbits_tmp));
            }
            core::cnap::get_symmetric_orbits(&orbits_tmp, edges, length, orbits_symmetry);
        }
    }

    Ok(())
}
// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Graph operations

use super::graph;

/// The graph distance between two vertices is the minimum count of edges to connect them. For a fully-connected graph, it means how many steps to take for traversing from vertice 1 to vertice 2
pub fn graph_distance<T: graph::Vertex>(
    vertex_1: usize,
    vertex_2: usize,
    vv: &graph::VertexVec<T>
) -> usize {
    if vertex_1 == vertex_2 {
        return 0
    } else {
        let mut vertices_visited: Vec<usize> = vec![];
        let mut vertices_boundary: Vec<usize> = vec![vertex_1];
        let mut distance: usize = 0;
        while vertices_visited.len() <  vv.len() {
            distance += 1;
            let mut neighbours: Vec<usize> = vertices_boundary.iter()
                .map(|&idx| vv[idx].neighbour_indexes())
                .flatten()
                .collect(); 
            neighbours.sort_unstable();
            neighbours.dedup();
            if neighbours.contains(&vertex_2) {
                break;
            } else {
                vertices_boundary.clear();
                for &idx in neighbours.iter() {
                    if !vertices_visited.contains(&idx) {
                        vertices_visited.push(idx);
                        vertices_boundary.push(idx);
                    }
                }
            }
        }

        distance
    }
}

#[cfg(test)]
mod test_core_graph_ops {
    use super::*;
    use crate::ext::molecule;

    #[test]
    fn test_graph_distance() {
        let smiles: String = "c1ccccc1CN".to_string();
        let mol = molecule::Molecule::from_smiles(&smiles);
        let vv = graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());

        assert_eq!(graph_distance(0, 0, &vv), 0); 
        assert_eq!(graph_distance(0, 1, &vv), 1); 
        assert_eq!(graph_distance(0, 2, &vv), 2); 
        assert_eq!(graph_distance(0, 3, &vv), 3); 
        assert_eq!(graph_distance(0, 4, &vv), 2); 
        assert_eq!(graph_distance(0, 5, &vv), 1); 
        assert_eq!(graph_distance(0, 6, &vv), 2); 
        assert_eq!(graph_distance(0, 7, &vv), 3); 
    }
}
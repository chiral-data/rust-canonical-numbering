// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::core;
use super::mapping_ops;
use super::case_breakable;
use super::case_cyclic;

const MAX_REDUCIBLE_CYCLE_SIZE: usize = 8;

/// for (k, v) in returned HashMap, vertex k is removed and reduced into a custom vertex with index v
fn get_reduced_groups(
    mapping: &mapping_ops::VertexMapping,
    boundary_edges: &Vec<Vec<mapping_ops::SimpleEdge>>
) -> std::collections::HashMap<usize, usize> {
    let keys: Vec<usize> = boundary_edges.iter().map(|bes| bes[0].0).collect();
    let mut reduced_groups: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for idx in 0..keys.len() {
        for jdx in 0..mapping[idx].len() {
            reduced_groups.insert(mapping[idx][jdx], keys[idx]);
        }
    }

    reduced_groups
}

fn new_vertex_index(
    vertex_index: usize,
    reduced_groups: &std::collections::HashMap<usize, usize> 
) -> usize {
    match reduced_groups.get(&vertex_index) {
        Some(&target_index) => target_index,
        None => vertex_index
    }
}

fn new_boundary_edge(
    edge: &mapping_ops::SimpleEdge,
    reduced_groups: &std::collections::HashMap<usize, usize> 
) -> mapping_ops::SimpleEdge {
    (new_vertex_index(edge.0, reduced_groups), new_vertex_index(edge.1, reduced_groups))
}

fn get_new_boundary_edges(
    edges: &Vec<mapping_ops::SimpleEdge>,
    reduced_groups: &std::collections::HashMap<usize, usize> 
) -> Vec<mapping_ops::SimpleEdge> {
    edges.iter()
        .map(|edge| new_boundary_edge(edge, reduced_groups))
        .collect()
}

fn get_symmetric_index(
    vertex_index: usize,
    orbits_symmetry: &Vec<core::orbit_ops::Orbit>,
    reduced_groups: &std::collections::HashMap<usize, usize>, 
) -> usize {
    let vertex_new = match reduced_groups.get(&vertex_index) {
        Some(&reduced_vertex) => reduced_vertex,
        None => vertex_index
    };
    for orbit in orbits_symmetry.iter() {
        if orbit.contains(&vertex_new) {
            return orbit[0];
        }
    }

    vertex_new
}

/// ReducibleGraph: struct & implementation
#[derive(Debug, Clone)]
pub struct ReducibleGraph<T: core::graph::Vertex> {
    pub vv: core::graph::VertexVec<T>,
    pub mapping: mapping_ops::VertexMapping,
    pub boundary_edges: mapping_ops::BoundaryEdges, 
    pub orbits_after_partition: Vec<core::orbit_ops::Orbit>,
    pub numbering: Vec<usize>,
}


impl<T: core::graph::Vertex> ReducibleGraph<T> {
    pub fn init(
        vv: core::graph::VertexVec<T>,
        mapping: mapping_ops::VertexMapping,
        boundary_edges: Vec<Vec<(usize, usize)>>,
        orbits_after_partition: Vec<core::orbit_ops::Orbit>,
        numbering: Vec<usize>,
    ) -> Self {
        Self { vv, mapping, boundary_edges, orbits_after_partition, numbering }
    }

    pub fn create_mapping(&mut self) -> Result<(), String> {
        if cfg!(debug_assertions) {
            println!("\nOrbits after partition: {:?}, total: {}", self.orbits_after_partition, self.orbits_after_partition.iter().map(|orbit| orbit.len()).collect::<Vec<usize>>().iter().sum::<usize>());
            println!("All {} vertices: {:?}", self.vv.len(), self.vv.valid_indexes());
        }

        match case_breakable::create_mapping(&self.orbits_after_partition, &self.numbering, &self.vv) {
            Ok((mapping, boundary_edges)) => {
                self.mapping = mapping;
                self.boundary_edges = boundary_edges;
                return Ok(());
            },
            Err(_) => ()
        }

        match case_cyclic::create_mapping(&self.orbits_after_partition, &self.numbering, &self.vv, MAX_REDUCIBLE_CYCLE_SIZE) {
            Ok((mapping, boundary_edges)) => {
                self.mapping = mapping;
                self.boundary_edges = boundary_edges;
                return Ok(())
            },
            Err(_) => ()
        }

        Err("Reduce Mapping Error".to_string())
    }

    pub fn get_simplified_graph(&self, custom_marker: &mut usize) -> Self {
        if self.mapping.len() == 0 {
            panic!("Reduce mapping not done")
        }

        if self.mapping.len() == 1 {
            panic!("Cannot reduce only 1 vertex")
        }

        let removed_indexes: Vec<usize> = self.mapping.to_vec().into_iter().flatten().collect();
        let mut indexes_new: Vec<usize> = self.vv.valid_indexes().to_vec().into_iter()
            .filter(|idx| !removed_indexes.contains(idx))
            .collect();
        let mut vertices_new = self.vv.all_vertices().to_vec();

        let reduced_groups: std::collections::HashMap<usize, usize> = get_reduced_groups(&self.mapping, &self.boundary_edges);
        for &vi in indexes_new.iter() {
            vertices_new[vi].update_edges_in_reduced_graph(vi, &reduced_groups, &self.numbering);
        }

        *custom_marker += 1;
        for idx in 0..self.boundary_edges.len() {
            let vertex_index: usize = self.boundary_edges[idx][0].0;
            let custom_vertice = T::custom_new_in_reduced_graph(
                vertex_index, *custom_marker, &self.boundary_edges[idx], &get_new_boundary_edges(&self.boundary_edges[idx], &reduced_groups), self.vv.all_vertices(), &self.numbering
            ); 
            indexes_new.push(vertex_index);
            vertices_new[vertex_index] = custom_vertice;
        }

        let new_vv = core::graph::VertexVec::init(indexes_new, vertices_new);
        ReducibleGraph {
            vv: new_vv,
            mapping: vec![],
            boundary_edges: vec![],
            orbits_after_partition: vec![],
            numbering: vec![] 
        }
    }

    pub fn get_folded_subgraph(&self, 
        custom_marker: &mut usize,
        orbits_symmetry: &Vec<core::orbit_ops::Orbit>
    ) -> Self {
        let mut indexes_new: Vec<usize> = self.mapping[0].clone();
        let mut vertices_new = self.vv.all_vertices().to_vec();
        
        let mut mapping_sliced = self.mapping.clone(); mapping_sliced.remove(0);
        let mut boundary_edges_sliced = self.boundary_edges.clone(); boundary_edges_sliced.remove(0);
        let reduced_groups: std::collections::HashMap<usize, usize> = get_reduced_groups(&mapping_sliced, &boundary_edges_sliced);
        let mut edges_with_value: Vec<(usize, usize, usize)> = self.boundary_edges[0].iter()
            .map(|edge| (edge.0, edge.1, get_symmetric_index(edge.1, orbits_symmetry, &reduced_groups)))
            .collect();
        edges_with_value.sort_by_key(|edge| edge.2); 

        *custom_marker += 1;
        let mut current_symmetric_index = edges_with_value[0].2;
        for edge in edges_with_value.iter() {
            if edge.2 > current_symmetric_index {
                current_symmetric_index = edge.2;
                *custom_marker += 1;
            }

            indexes_new.push(edge.1);
            let custom_vertex = T::custom_new_in_folded_graph(*custom_marker, &(edge.0, edge.1), self.vv.all_vertices());
            vertices_new[edge.1] = custom_vertex;
        }

        let new_vv = core::graph::VertexVec::init(indexes_new, vertices_new);
        ReducibleGraph {
            vv: new_vv,
            mapping: vec![],
            boundary_edges: vec![],
            orbits_after_partition: vec![],
            numbering: vec![] 
        }
    }

    pub fn reduced_subgraph_indexes(&self) -> Vec<usize> {
        self.boundary_edges.iter()
            .map(|bes| bes[0].0)
            .collect()
    }

    pub fn symmetric_orbits_from_folded_graph(&self) -> Vec<core::orbit_ops::Orbit> {
        (0..self.mapping[0].len())
            .map(|idx| self.mapping.iter().map(|m| m[idx]).collect())
            .collect()
    }

    pub fn is_reducible(&self) -> bool {
        self.mapping[0].len() > 1
    }
}

impl<T: core::graph::Vertex> std::fmt::Display for ReducibleGraph<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut output_content: String = "".to_string();
        output_content += &self.vv.valid_indexes().iter().map(|vi| vi.to_string()).collect::<Vec<String>>().join("\t");
        output_content += "\n";
        output_content += &self.vv.valid_indexes().iter().map(|&vi| self.numbering[vi].to_string()).collect::<Vec<String>>().join("\t");
        output_content += "\n";
        write!(f, "{}", format!("{}", output_content))
    }
}

#[cfg(test)]
mod test_reduce_graph {
    use crate::ext::molecule;
    use super::*;

    #[test]
    fn test_reducible_graph() {
        let smiles_vec: Vec<String> = vec![
            "O=C(O)CCCC(=O)NCCCC[C@@H](C(=O)NCCCCCCCCCCCCNC(=O)[C@@H](CCCCNC(=O)CCCC(=O)O)N(Cc1ccc(OCc2ccccc2)cc1)Cc1ccc(OCc2ccccc2)cc1)N(Cc1ccc(OCc2ccccc2)cc1)Cc1ccc(OCc2ccccc2)cc1",
            "CCOP(=O)(OCC)OCc1ccc(S(=O)(=O)CC(C[C@H]2O[C@@H]3C[C@]4(C)CC[C@H]5[C@H](C)CC[C@@H]([C@H]2C)[C@]53OO4)C[C@H]2O[C@@H]3C[C@]4(C)CC[C@H]5[C@H](C)CC[C@@H]([C@H]2C)[C@]53OO4)cc1",
            "CC(C)(C)c1cc2c(OCCCCNC(=N)N)c(c1)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)C2",
            "O=c1cc(-c2ccc(OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOc3ccc(-c4cc(=O)c5ccccc5o4)cc3)cc2)oc2ccccc12",
            // "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)n5)n4)n3",
        ].iter().map(|s| s.to_string()).collect();
        let results = vec![
            (vec![22, 23], vec![22, 21, 20, 19, 18, 17, 16, 14, 15, 13, 12, 77, 11, 78, 93, 10, 79, 94, 9, 92, 80, 107, 95, 8, 91, 81, 106, 96, 6, 82, 97, 7, 5, 83, 98, 4, 84, 99, 3, 85, 100, 1, 90, 86, 105, 101, 2, 0, 89, 87, 104, 102, 88, 103, 23], 202),
            (vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 59, 60, 19, 39], vec![19, 20, 21, 34, 22, 35, 33, 23, 36, 32, 24, 37, 28, 31, 25, 26, 38, 27, 29, 30, 18], 202),
            (vec![0, 1, 2, 3, 4, 5, 6, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 79, 7, 29, 49, 69], vec![7, 8, 9, 10, 11, 12, 13, 14, 16, 15, 6, 17], 202),
            (vec![27, 28], vec![27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 65, 5, 66, 4, 3, 2, 67, 1, 68, 0, 73, 69, 72, 70, 71, 28], 202),
            // (vec![3, 8, 13, 18, 19, 14, 12, 4], vec![10, 11, 9, 12, 22, 8, 13], 202),
        ];

        for (idx, smiles) in smiles_vec.iter().enumerate() {
            let mol = molecule::Molecule::from_smiles(smiles);
            if cfg!(debug_assertions) {
                println!("{}", mol.smiles_with_index(smiles, &vec![]));
            }

            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut rg = ReducibleGraph {
                vv: vv,
                mapping: vec![],
                boundary_edges: vec![],
                orbits_after_partition: vec![],
                numbering: vec![]
            };
            
            let mut numbering_rg: Vec<usize> = vec![];
            let mut orbits_after_partition_rg: Vec<core::orbit_ops::Orbit> = vec![];
            core::givp::run::<molecule::AtomExtendable>(&rg.vv, &mut numbering_rg, &mut orbits_after_partition_rg);
            rg.numbering = numbering_rg;
            rg.orbits_after_partition = orbits_after_partition_rg;
            match rg.create_mapping() {
                Ok(_) => {
                    let mut custom_marker: usize = 200;
                    let sg = rg.get_simplified_graph(&mut custom_marker);
                    let fg = rg.get_folded_subgraph(&mut custom_marker, &rg.orbits_after_partition);
                    assert_eq!(
                        (sg.vv.valid_indexes().to_vec(), fg.vv.valid_indexes().to_vec(), custom_marker),
                        results[idx]
                    );
                },
                Err(_) => ()
            }
        }
    }
}
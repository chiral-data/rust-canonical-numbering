// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Case breakable

use crate::core;
use super::mapping_ops;

fn is_neighbour_breakable<T: core::graph::Vertex>(
    vertex_index: usize,
    neighbour_index: usize,
    vertices_unexpected: &Vec<usize>,
    vertex_vec: &core::graph::VertexVec<T>
) -> bool {
    let mut vertices_visited: Vec<usize> = vec![vertex_index];
    let mut vertices_boundary: Vec<usize> = vertex_vec[neighbour_index].neighbour_indexes().into_iter()
        .filter(|&idx| idx != vertex_index)
        .collect();

    for vn in vertices_unexpected.iter() {
        if vertices_boundary.contains(vn) {
            return false;
        }
    }
    
    while vertices_boundary.len() > 0 {
        let neighbours: Vec<usize> = vertices_boundary.iter()
            .map(|&idx| vertex_vec[idx].neighbour_indexes())
            .flatten()
            .filter(|idx| !vertices_visited.contains(idx))
            .filter(|idx| !vertices_boundary.contains(idx))
            .collect();

        for vn in vertices_unexpected.iter() {
            if neighbours.contains(vn) {
                return false;
            }
        }

        vertices_visited.append(&mut vertices_boundary);
        vertices_boundary = neighbours;
    }

    true
}

fn get_breakable_neighbours<T: core::graph::Vertex>(
    vertex_index: usize,
    vertices_unexpected: &Vec<usize>,
    vv: &core::graph::VertexVec<T>
) -> Vec<usize> {
    vv[vertex_index].neighbour_indexes().into_iter()
        .filter(|&vi| is_neighbour_breakable(vertex_index, vi, vertices_unexpected, vv))
        .collect()
}

/// The maximum graph distance between two vertices in an orbit
/// Normally every two-vertex pair should be considered.
/// As the orbit is resulted from givp process, we only consider the distances between the first vertex and the other vertices
fn max_graph_distance_inside_orbit<T: core::graph::Vertex>(
    orbit: &Vec<usize>,
    vv: &core::graph::VertexVec<T>
) -> usize {
    let mut distance: usize = 0;
    for idx in 1..orbit.len() {
        let distance_tmp: usize = core::graph_ops::graph_distance(orbit[0], orbit[idx], vv);
        if distance_tmp > distance {
            distance = distance_tmp;
        }
    }

    distance
}

/// return start orbit and mapping starting vertices
fn find_starting_orbit<T: core::graph::Vertex>(
    orbits_after_partition: &Vec<core::orbit_ops::Orbit>,
    vv: &core::graph::VertexVec<T>
) -> (core::orbit_ops::Orbit, Vec<Vec<usize>>) {
    let mut orbits_cloned = orbits_after_partition.to_vec();
    orbits_cloned.sort_by_cached_key(
        |orbit| max_graph_distance_inside_orbit(orbit, vv)
    );

    for orbit in orbits_cloned.iter() {
        if vv[orbit[0]].degree() == 1 {
            continue;
        }

        if get_breakable_neighbours(orbit[0], orbit, vv).len() > 0 {
            return (orbit.clone(), orbit.iter().map(|&vi| get_breakable_neighbours(vi, orbit, vv)).collect());
        }
    }

    (vec![], vec![])
}


fn get_boundary_edges_for_breakable_subgraph<T: core::graph::Vertex>(
    vertex_index: usize,
    starting_neighbours: &mapping_ops::NeighbourIndexes,
    vv: &core::graph::VertexVec<T>
) -> Vec<mapping_ops::SimpleEdge> {
    vv[vertex_index].neighbour_indexes().iter()
        .filter(|nvi| !starting_neighbours.contains(nvi))
        .map(|&nvi| (vertex_index, nvi))
        .collect()
}

pub fn create_mapping<T: core::graph::Vertex>(
    orbits_after_partition: &Vec<core::orbit_ops::Orbit>,
    numbering: &Vec<usize>,
    vv: &core::graph::VertexVec<T>
) -> Result<(mapping_ops::VertexMapping, mapping_ops::BoundaryEdges), String> {
    let (starting_orbit, starting_neighbours) = find_starting_orbit(orbits_after_partition, vv);
    if starting_orbit.len() == 0 {
        Err(format!("Breakable Mapping: cannot find a starting orbit from orbits: {:?}", orbits_after_partition))
    }
    else {
        if cfg!(debug_assertions) {
            println!("Breakable Mapping: starting orbit {:?}", starting_orbit);
            println!("Breakable Mapping: starting neighbours: {:?}", starting_neighbours);
        }

        let boundary_edges = (0..starting_orbit.len())
            .map(|idx| get_boundary_edges_for_breakable_subgraph(starting_orbit[idx], &starting_neighbours[idx], vv))
            .collect();
        let mut mapping: mapping_ops::VertexMapping = starting_orbit.iter()
            .map(|&v_idx| vec![v_idx])
            .collect();
        let mut mapping_stack: Vec<mapping_ops::MappingStatus> = vec![];
        let mut mapped_cur: usize = 0;

        match mapping_ops::mapping_neighbours(&mut starting_neighbours.clone(), numbering) {
            Ok(neighbour_mappings) => {
                mapping_ops::update_mapping_stack(mapped_cur, &mapping, &boundary_edges, &neighbour_mappings, &mut mapping_stack);
                if let Some((_, new_mapping, _)) = mapping_stack.pop() {
                    mapping = new_mapping;
                } else {
                    panic!("successful neighbour mapping should have at least one mapping!")
                }
            },
            Err(_) => {
                panic!("Error on mapping starting neighbours")
            }
        }

        // Mapping Process
        mapped_cur += 1;
        while mapped_cur < mapping[0].len() {
            let mut neighbours: Vec<mapping_ops::NeighbourIndexes> = (0..mapping.len())
                .map(|idx| mapping_ops::find_new_neighbours(mapping[idx][mapped_cur], &mapping[idx], vv))
                .collect();

            if neighbours[0].len() > 0 {
                match mapping_ops::mapping_neighbours(&mut neighbours, numbering) {
                    Ok(neighbour_mappings) => {
                        mapping_ops::update_mapping_stack(mapped_cur, &mapping, &boundary_edges, &neighbour_mappings, &mut mapping_stack);
                    },
                    Err(_) => {
                    }
                }

                if let Some((new_mapped_cur, new_mapping, _)) = mapping_stack.pop() {
                    mapping = new_mapping;
                    mapped_cur = new_mapped_cur;
                } else {
                    panic!("mapping status runs out! mapping fails!")
                }
            } else {
                mapped_cur += 1;
            }
        }

        Ok((mapping, boundary_edges))
    }
}


#[cfg(test)]
mod test_reduce_case_breakable {
    use super::*;
    use crate::ext::molecule;

    #[test]
    fn test_is_neighbour_breakable() {
        let smiles_vec: Vec<String> = vec![
            r#"COc1cc(Cc2cnc(/N=C3\C(=O)N(CN(Cc4ccccc4)Cc4ccccc4)c4ccc(Cl)cc43)nc2N)cc(OC)c1OC"#, // chembl 2209
            "CCCC[C@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)CCNC(=O)CCNC(=S)Nc1ccc(O)c(NC(=O)CNC(=O)CSC(c2ccccc2)(c2ccccc2)c2ccccc2)c1)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(N)=O", // chembl 2067
        ].iter().map(|s| s.to_string()).collect();
        let tests_params: Vec<Vec<(usize, usize, Vec<usize>, bool)>> = vec![
            vec![
                (2, 1, vec![2, 45, 42], true),
                (2, 3, vec![2, 45, 42], false),
                (2, 45, vec![2, 45, 42], false)
            ],
            vec![
                (61, 62, vec![61, 49, 55], true),
                (61, 66, vec![61, 49, 55], true),
                (61, 48, vec![61, 49, 55], false),
            ],
        ];

        for (idx, smiles) in smiles_vec.iter().enumerate() {
            let mol = molecule::Molecule::from_smiles(smiles);
             println!("{}", mol.smiles_with_index(smiles, &vec![]));
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            for tp in tests_params[idx].iter() {
                assert_eq!(is_neighbour_breakable(tp.0, tp.1, &tp.2, &vv), tp.3);
            }
        }
    }

    #[test]
    fn test_find_starting_orbit() {
        type InputType1 = String;
        type ReturnType = (core::orbit_ops::Orbit, Vec<Vec<usize>>);
        let test_data: Vec<(InputType1, ReturnType)> = vec![
            (
                "O=C(O)CCCC(=O)NCCCC[C@@H](C(=O)NCCCCCCCCCCCCNC(=O)[C@@H](CCCCNC(=O)CCCC(=O)O)N(Cc1ccc(OCc2ccccc2)cc1)Cc1ccc(OCc2ccccc2)cc1)N(Cc1ccc(OCc2ccccc2)cc1)Cc1ccc(OCc2ccccc2)cc1",
                (vec![22, 23], vec![vec![21], vec![24]]),
            ),
            (
                "CCOP(=O)(OCC)OCc1ccc(S(=O)(=O)CC(C[C@H]2O[C@@H]3C[C@]4(C)CC[C@H]5[C@H](C)CC[C@@H]([C@H]2C)[C@]53OO4)C[C@H]2O[C@@H]3C[C@]4(C)CC[C@H]5[C@H](C)CC[C@@H]([C@H]2C)[C@]53OO4)cc1",
                (vec![19, 39], vec![vec![20], vec![40]]),
            ),
            (
                "CC(C)(C)c1cc2c(OCCCCNC(=N)N)c(c1)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)C2",
                (vec![7, 29, 49, 69], vec![vec![8], vec![30], vec![50], vec![70]]),
            ),
            (
                "O=c1cc(-c2ccc(OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOc3ccc(-c4cc(=O)c5ccccc5o4)cc3)cc2)oc2ccccc12",
                (vec![27, 28], vec![vec![26], vec![29]]),
            ),
            (
                "CC1=C(CC[C@H](C)CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)O[C@H]2C[C@H]3[C@@H]4CC[C@@H]5C[C@@H](O[C@@H]6O[C@H](CO)[C@H](O)[C@H](O)[C@H]6O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)CC[C@]5(C)[C@H]4CC[C@]3(C)[C@@H]12", 
                (vec![], vec![]),
            ),
            (
                "CCn1c2ccc3cc2c2cc(ccc21)C(=O)c1ccc(cc1)Cn1cc[n+](c1)Cc1ccc(cc1)C(=O)c1ccc2c(c1)c1cc(ccc1n2CC)C(=O)c1ccc(cc1)C[n+]1ccn(c1)Cc1ccc(cc1)C3=O",
                (vec![15, 74], vec![vec![16], vec![75]]),
            ),
            (
                "CCCCCCCCOc1c2cc(C(=O)O)cc1Cc1cc(C(=O)O)cc(c1OCCCCCCCC)Cc1cc(C(=O)O)cc(c1OCCCCCCCC)Cc1cc(cc(C(=O)O)c1)C2",
                (vec![9, 46], vec![vec![8], vec![47]]),
            ),
            (
                "O=C(Nc1ccc([N+](=O)[O-])cc1)OCCN1CCN(CCOC(=O)Nc2ccc([N+](=O)[O-])cc2)CCN(CCOC(=O)Nc2ccc([N+](=O)[O-])cc2)CCN(CCOC(=O)Nc2ccc([N+](=O)[O-])cc2)CC1",
                (vec![15, 18, 36, 54], vec![vec![14], vec![19], vec![37], vec![55]]),
            ),
            (
                "CC(C)C[C@@H]1NC(=O)[C@H](CCCN)NC(=O)[C@H](C(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](C2c3ccccc3-c3ccccc32)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCN)NC(=O)[C@H](C(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](C2c3ccccc3-c3ccccc32)NC1=O",
                (vec![4, 47], vec![vec![3], vec![48]]),
            ),
            (
                "CC(=O)O[C@H]1[C@H](C2=CC(=O)c3ccccc3C2=O)O[C@H](Cn2cc(COC(=O)c3cccc(C(=O)OCc4cn(C[C@H]5O[C@@H](C6=CC(=O)c7ccccc7C6=O)[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@H]5OC(C)=O)nn4)c3)nn2)[C@H](OC(C)=O)[C@@H]1OC(C)=O",
                (vec![28, 32], vec![vec![26], vec![33]]),
            ),
            (
                "c1ccc(Cc2ccc[n+](CCCCCc3cc(CCCCC[n+]4cccc(Cc5ccccc5)c4)c(CCCCC[n+]4cccc(Cc5ccccc5)c4)cc3CCCCC[n+]3cccc(Cc4ccccc4)c3)c2)cc1",
                (vec![15, 17, 36, 56], vec![vec![14], vec![18], vec![37], vec![57]]),
            ),
            (
                "NCCCNCCCCN(CCCN)C(=O)CCCCNC(=O)c1ccc(-c2c3nc(c(-c4ccc(C(=O)NCCCCC(=O)N(CCCN)CCCCNCCCN)cc4)c4ccc([nH]4)c(-c4ccc(C(=O)NCCCCC(=O)N(CCCN)CCCCNCCCN)cc4)c4nc(c(-c5ccc(C(=O)NCCCCC(=O)N(CCCN)CCCCNCCCN)cc5)c5ccc2[nH]5)C=C4)C=C3)cc1",
                (vec![27, 31, 66, 99], vec![vec![26], vec![32], vec![67], vec![100]]),
            ),
            (
                "CC(C)C[C@H]1C(=O)N(C)CC(=O)N(C)[C@@H](C(C)C)C(=O)NC[C@H](NC(=O)c2ccc3ccccc3n2)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)CC(=O)N(C)[C@@H](C(C)C)C(=O)NC[C@H](NC(=O)c2ccc3ccccc3n2)C(=O)N1C",
                (vec![4, 40], vec![vec![3], vec![41]]),
            ),
            (
                "CCCC[C@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)CCNC(=O)CCNC(=S)Nc1ccc(O)c(NC(=O)CNC(=O)CSC(c2ccccc2)(c2ccccc2)c2ccccc2)c1)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(N)=O", // chembl 2067
                (vec![49, 55, 61], vec![vec![54, 50], vec![60, 56], vec![66, 62]]),
            ),
        ].iter().map(|td| (td.0.to_string(), td.1.clone())).collect();

        for td in test_data.iter() {
            let (smiles, results) = td;
            let mol = molecule::Molecule::from_smiles(smiles);
            println!("{}", mol.smiles_with_index(smiles, &vec![]));
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut numbering: Vec<usize> = vec![];
            let mut orbits_after_partition: Vec<core::orbit_ops::Orbit> = vec![];
            core::givp::run::<molecule::AtomExtendable>(&vv, &mut numbering, &mut orbits_after_partition);
            assert_eq!(find_starting_orbit(&orbits_after_partition, &vv), *results);
        }
    }

    #[test]
    fn test_mapping() {
        type InputType1 = String;
        type ReturnType = (mapping_ops::VertexMapping, mapping_ops::BoundaryEdges);
        let test_data: Vec<(InputType1, ReturnType)> = vec![
            (
                "O=C(O)CCCC(=O)NCCCC[C@@H](C(=O)NCCCCCCCCCCCCNC(=O)[C@@H](CCCCNC(=O)CCCC(=O)O)N(Cc1ccc(OCc2ccccc2)cc1)Cc1ccc(OCc2ccccc2)cc1)N(Cc1ccc(OCc2ccccc2)cc1)Cc1ccc(OCc2ccccc2)cc1",
                (vec![vec![22, 21, 20, 19, 18, 17, 16, 14, 15, 13, 12, 77, 11, 78, 93, 10, 79, 94, 9, 92, 80, 107, 95, 8, 91, 81, 106, 96, 6, 82, 97, 7, 5, 83, 98, 4, 84, 99, 3, 85, 100, 1, 90, 86, 105, 101, 2, 0, 89, 87, 104, 102, 88, 103], vec![23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 46, 34, 62, 47, 35, 63, 48, 36, 64, 76, 49, 61, 37, 65, 75, 50, 60, 38, 66, 51, 39, 40, 67, 52, 41, 68, 53, 42, 69, 54, 43, 70, 74, 55, 59, 45, 44, 71, 73, 56, 58, 72, 57]], vec![vec![(22, 23)], vec![(23, 22)]]),
            ),
            (
                "CCOP(=O)(OCC)OCc1ccc(S(=O)(=O)CC(C[C@H]2O[C@@H]3C[C@]4(C)CC[C@H]5[C@H](C)CC[C@@H]([C@H]2C)[C@]53OO4)C[C@H]2O[C@@H]3C[C@]4(C)CC[C@H]5[C@H](C)CC[C@@H]([C@H]2C)[C@]53OO4)cc1",
                (vec![vec![19, 20, 21, 34, 22, 35, 33, 23, 36, 32, 24, 37, 28, 31, 25, 26, 38, 27, 29, 30], vec![39, 40, 41, 54, 42, 55, 53, 43, 56, 52, 44, 57, 48, 51, 45, 46, 58, 47, 49, 50]], vec![vec![(19, 18)], vec![(39, 18)]]),
            ),
            (
                "CC(C)(C)c1cc2c(OCCCCNC(=N)N)c(c1)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)C2",
                (vec![vec![7, 8, 9, 10, 11, 12, 13, 14, 16, 15], vec![29, 30, 31, 32, 33, 34, 35, 36, 38, 37], vec![49, 50, 51, 52, 53, 54, 55, 56, 58, 57], vec![69, 70, 71, 72, 73, 74, 75, 76, 78, 77]], vec![vec![(7, 6), (7, 17)], vec![(29, 28), (29, 20)], vec![(49, 48), (49, 40)], vec![(69, 68), (69, 60)]]),
            ),
            (
                "O=c1cc(-c2ccc(OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOc3ccc(-c4cc(=O)c5ccccc5o4)cc3)cc2)oc2ccccc12",
                (vec![vec![27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 65, 5, 66, 4, 3, 2, 67, 1, 68, 0, 73, 69, 72, 70, 71], vec![28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 64, 50, 63, 51, 52, 53, 62, 54, 61, 55, 56, 60, 57, 59, 58]], vec![vec![(27, 28)], vec![(28, 27)]]),
            ),
            (
                "CCn1c2ccc3cc2c2cc(ccc21)C(=O)c1ccc(cc1)Cn1cc[n+](c1)Cc1ccc(cc1)C(=O)c1ccc2c(c1)c1cc(ccc1n2CC)C(=O)c1ccc(cc1)C[n+]1ccn(c1)Cc1ccc(cc1)C3=O",
                (vec![vec![15, 16], vec![74, 75]], vec![vec![(15, 11), (15, 17)], vec![(74, 71), (74, 6)]]),
            ),
            (
                "CCCCCCCCOc1c2cc(C(=O)O)cc1Cc1cc(C(=O)O)cc(c1OCCCCCCCC)Cc1cc(C(=O)O)cc(c1OCCCCCCCC)Cc1cc(cc(C(=O)O)c1)C2",
                (vec![vec![9, 8, 7, 6, 5, 4, 3, 2, 1, 0], vec![46, 47, 48, 49, 50, 51, 52, 53, 54, 55]], vec![vec![(9, 17), (9, 10)], vec![(46, 45), (46, 38)]]),
            ),
            (
                "O=C(Nc1ccc([N+](=O)[O-])cc1)OCCN1CCN(CCOC(=O)Nc2ccc([N+](=O)[O-])cc2)CCN(CCOC(=O)Nc2ccc([N+](=O)[O-])cc2)CCN(CCOC(=O)Nc2ccc([N+](=O)[O-])cc2)CC1",
                (vec![vec![15, 14, 13, 12, 1, 0, 2, 3, 11, 4, 10, 5, 6, 7, 9, 8], vec![18, 19, 20, 21, 22, 23, 24, 25, 26, 33, 27, 32, 28, 29, 31, 30], vec![36, 37, 38, 39, 40, 41, 42, 43, 51, 44, 50, 45, 46, 47, 49, 48], vec![54, 55, 56, 57, 58, 59, 60, 61, 62, 69, 63, 68, 64, 65, 67, 66]], vec![vec![(15, 71), (15, 16)], vec![(18, 17), (18, 34)], vec![(36, 35), (36, 52)], vec![(54, 53), (54, 70)]]),
            ),
            (
                "CC(C)C[C@@H]1NC(=O)[C@H](CCCN)NC(=O)[C@H](C(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](C2c3ccccc3-c3ccccc32)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCN)NC(=O)[C@H](C(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](C2c3ccccc3-c3ccccc32)NC1=O",
                (vec![vec![4, 3, 1, 0, 2], vec![47, 48, 49, 51, 50]], vec![vec![(4, 92), (4, 5)], vec![(47, 45), (47, 52)]]),
            ),
            (
                "CC(=O)O[C@H]1[C@H](C2=CC(=O)c3ccccc3C2=O)O[C@H](Cn2cc(COC(=O)c3cccc(C(=O)OCc4cn(C[C@H]5O[C@@H](C6=CC(=O)c7ccccc7C6=O)[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@H]5OC(C)=O)nn4)c3)nn2)[C@H](OC(C)=O)[C@@H]1OC(C)=O",
                (vec![vec![28, 26, 27, 25, 24, 23, 22, 74, 21, 75, 20, 19, 18, 76, 5, 77, 81, 4, 6, 78, 82, 3, 7, 16, 79, 80, 83, 1, 8, 17, 15, 84, 85, 0, 2, 9, 10, 14, 11, 13, 12], vec![32, 33, 34, 35, 36, 37, 38, 72, 39, 71, 40, 41, 42, 66, 43, 67, 61, 56, 44, 68, 62, 57, 45, 54, 69, 70, 63, 58, 46, 55, 53, 64, 65, 59, 60, 47, 48, 52, 49, 51, 50]], vec![vec![(28, 73), (28, 29)], vec![(32, 31), (32, 73)]]),
            ),
            (
                "c1ccc(Cc2ccc[n+](CCCCCc3cc(CCCCC[n+]4cccc(Cc5ccccc5)c4)c(CCCCC[n+]4cccc(Cc5ccccc5)c4)cc3CCCCC[n+]3cccc(Cc4ccccc4)c3)c2)cc1",
                (vec![vec![15, 14, 13, 12, 11, 10, 9, 8, 75, 7, 5, 6, 4, 3, 2, 76, 1, 77, 0], vec![17, 18, 19, 20, 21, 22, 23, 24, 35, 25, 27, 26, 28, 29, 30, 34, 31, 33, 32], vec![36, 37, 38, 39, 40, 41, 42, 43, 54, 44, 46, 45, 47, 48, 53, 49, 52, 50, 51], vec![56, 57, 58, 59, 60, 61, 62, 63, 74, 64, 66, 65, 67, 68, 69, 73, 70, 72, 71]], vec![vec![(15, 56), (15, 16)], vec![(17, 16), (17, 36)], vec![(36, 17), (36, 55)], vec![(56, 55), (56, 15)]]),
            ),
            (
                "NCCCNCCCCN(CCCN)C(=O)CCCCNC(=O)c1ccc(-c2c3nc(c(-c4ccc(C(=O)NCCCCC(=O)N(CCCN)CCCCNCCCN)cc4)c4ccc([nH]4)c(-c4ccc(C(=O)NCCCCC(=O)N(CCCN)CCCCNCCCN)cc4)c4nc(c(-c5ccc(C(=O)NCCCCC(=O)N(CCCN)CCCCNCCCN)cc5)c5ccc2[nH]5)C=C4)C=C3)cc1",
                (vec![vec![27, 26, 25, 138, 24, 139, 23, 21, 22, 20, 19, 18, 17, 16, 14, 15, 9, 10, 8, 11, 7, 12, 6, 13, 5, 4, 3, 2, 1, 0], vec![31, 32, 33, 60, 34, 59, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 50, 47, 51, 48, 52, 49, 53, 54, 55, 56, 57, 58], vec![66, 67, 95, 68, 94, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 85, 82, 86, 83, 87, 84, 88, 89, 90, 91, 92, 93], vec![99, 100, 101, 128, 102, 127, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 118, 115, 119, 116, 120, 117, 121, 122, 123, 124, 125, 126]], vec![vec![(27, 132), (27, 28)], vec![(31, 30), (31, 61)], vec![(66, 64), (66, 96)], vec![(99, 98), (99, 129)]]),
            ),
            (
                "CC(C)C[C@H]1C(=O)N(C)CC(=O)N(C)[C@@H](C(C)C)C(=O)NC[C@H](NC(=O)c2ccc3ccccc3n2)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)CC(=O)N(C)[C@@H](C(C)C)C(=O)NC[C@H](NC(=O)c2ccc3ccccc3n2)C(=O)N1C",
                (vec![vec![4, 3, 1, 0, 2], vec![40, 41, 42, 44, 43]], vec![vec![(4, 78), (4, 5)], vec![(40, 38), (40, 45)]]),
            )
        ].iter().map(|td| (td.0.to_string(), td.1.clone())).collect();

        for td in test_data.iter() {
            let (smiles, results) = td;
            let mol = molecule::Molecule::from_smiles(smiles);
            println!("{}", mol.smiles_with_index(smiles, &vec![]));
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut numbering: Vec<usize> = vec![];
            let mut orbits_after_partition: Vec<core::orbit_ops::Orbit> = vec![];
            core::givp::run::<molecule::AtomExtendable>(&vv, &mut numbering, &mut orbits_after_partition);

            match create_mapping(&orbits_after_partition, &numbering, &vv) {
                Ok(mapping_results) => { assert_eq!(mapping_results, *results); }
                Err(s) => { panic!("mapping failed: {:?}", s); }
            }
        }
    }

    #[test]
    fn test_molecules() {
        type ParamType1 = String;
        let test_data: Vec<ParamType1> = vec![
            "C1CCC2C(CC1)C1CCCCCC21", // an example for the paper
        ].iter().map(|td| td.to_string()).collect();

        for td in test_data.iter() {
            let smiles = td;
            let mol = molecule::Molecule::from_smiles(smiles);
            if cfg!(debug_assertions) {
                println!("{}", mol.smiles_with_index(smiles, &vec![]));
            }

            let mut orbits_partitioned: Vec<core::orbit_ops::Orbit> = vec![];
            let mut orbits_symmetry: Vec<core::orbit_ops::Orbit> = vec![];
            let mut numbering: Vec<usize> = vec![];
            molecule::canonical_numbering_and_symmetry_perception(&mol.atoms, &mut orbits_partitioned, &mut orbits_symmetry, &mut numbering);
            if cfg!(debug_assertions) {
                core::orbit_ops::orbits_sort(&mut orbits_partitioned);
                core::orbit_ops::orbits_sort(&mut orbits_symmetry);
                println!("GIVP vs CNAP\n{:?}\n{:?}", orbits_partitioned, orbits_symmetry);
            }
            assert_eq!(core::orbit_ops::orbits_equal(&orbits_partitioned, &orbits_symmetry), true);
        }
    }
}
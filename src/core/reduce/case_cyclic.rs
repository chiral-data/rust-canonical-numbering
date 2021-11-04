// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Case cyclic

use crate::core;
use super::mapping_ops;

fn filter_cyclic_orbits<T: core::graph::Vertex>(
    orbits: &Vec<core::orbit_ops::Orbit>,
    vv: &core::graph::VertexVec<T>,
) -> Vec<core::orbit_ops::Orbit> {
    orbits.to_vec().into_iter()
        .filter(|orbit| core::cycle_ops::cycle_size(orbit[0], vv) > 0)
        .collect()
}

fn get_cycle_index_for_vertex(
    vertex: usize,
    cycles: &Vec<Vec<usize>>
) -> usize {
    let indexes: Vec<usize> = (0..cycles.len())
        .filter(|&cycle_idx| cycles[cycle_idx].contains(&vertex))
        .collect();

    if indexes.len() == 0 {
        // panic!("Case cyclic: cannot find vertex {} in cycles", vertex);
        cycles.len() + 1
    } else if indexes.len() == 1 {
        indexes[0]
    } else {
        panic!("Case cyclic: more than 2 cylces contains vertex {}", vertex);
    }
}

fn divide_orbit_by_cycles(
    orbit: &core::orbit_ops::Orbit,
    cycles: &Vec<Vec<usize>>
) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    if get_cycle_index_for_vertex(orbit[0], cycles) == cycles.len() + 1 { // acyclic orbit
        return (vec![], vec![]);
    } 

    let cycle_indexes: Vec<usize> = orbit.iter().map(|&vi| get_cycle_index_for_vertex(vi, cycles)).collect();
    let mut vertices_by_cycle: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for idx in 0..orbit.len() {
        vertices_by_cycle.entry(cycle_indexes[idx]).or_insert(vec![]).push(orbit[idx]);
    }

    let mut orbit_divided: Vec<Vec<usize>> = vec![];
    let mut cycles_matched: Vec<Vec<usize>> = vec![];
    for (&key, val) in vertices_by_cycle.iter() {
        orbit_divided.push(val.to_vec());
        cycles_matched.push(cycles[key].clone());
    }

    (orbit_divided, cycles_matched)
}

fn find_starting_orbit<T: core::graph::Vertex>(
    orbits_after_partition: &Vec<core::orbit_ops::Orbit>,
    vv: &core::graph::VertexVec<T>,
    max_cycle_size: usize,
) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let orbits_cyclic = filter_cyclic_orbits(orbits_after_partition, vv);
    for orbit in orbits_cyclic.iter() {
        if core::cycle_ops::find_all_cycles(orbit[0], vv, max_cycle_size).len() > 1 { // skip the orbit who are intersections of cycles
            continue;
        }

        let mut cycles: Vec<Vec<usize>> = orbit.iter()
            .map(|&vi| core::cycle_ops::find_all_cycles(vi, vv, max_cycle_size))
            .flatten()
            .collect();
        core::orbit_ops::orbits_sort(&mut cycles);
        cycles.dedup();
        if cycles.len() == 1 {
            continue;
        }

        let mut cycle_overlapped: bool = false;
        for i in 0..cycles.len() {
            for j in (i+1)..cycles.len() {
                if core::orbit_ops::orbit_overlap(&cycles[i], &cycles[j]) {
                    cycle_overlapped = true;
                    break;
                } 
            }

            if cycle_overlapped { break; }
        }
        if cycle_overlapped {
            continue;
        }

        for cycle in cycles.iter_mut() {
            *cycle = cycle.clone().into_iter()
                .filter(|&vi| core::cycle_ops::find_all_cycles(vi, vv, max_cycle_size).len() == 1)
                .collect();
        }

        let (orbit_divided, cycles_matched) = divide_orbit_by_cycles(orbit, &cycles);
        if orbit_divided.len() > 1 {
            return (orbit_divided, cycles_matched)
        }
    }

    (vec![], vec![])
}


fn update_boundary_edges(
    mapping: &mapping_ops::VertexMapping,
    boundary_edges: &mut mapping_ops::BoundaryEdges,
    neighbours: &mut Vec<mapping_ops::NeighbourIndexes>,
    mapped_cur: usize,
    cycles: &Vec<Vec<usize>>,
) {
    for i in 0..mapping.len() {
        let mut indexes_to_remove: Vec<usize> = vec![];

        for j in 0..neighbours[i].len() {
            if cycles[i].contains(&neighbours[i][j]) {
                continue;
            } else {
                indexes_to_remove.push(j);
            }
        }
    
        indexes_to_remove.sort_unstable();
        indexes_to_remove.reverse();
        for &vi in indexes_to_remove.iter() {
            boundary_edges[i].push((mapping[i][mapped_cur], neighbours[i][vi]));
            neighbours[i].remove(vi);
        }
    }
}

pub fn create_mapping<T: core::graph::Vertex>(
    orbits_after_partition: &Vec<core::orbit_ops::Orbit>,
    numbering: &Vec<usize>,
    vv: &core::graph::VertexVec<T>,
    max_cycle_size: usize,
) -> Result<(mapping_ops::VertexMapping, mapping_ops::BoundaryEdges), String> {
    let (starting_orbit, cycles_to_reduce) = find_starting_orbit(orbits_after_partition, vv, max_cycle_size);
    if starting_orbit.len() == 0 {
        Err(format!("Cyclic Mapping: cannot find a starting orbit from orbits: {:?}", orbits_after_partition))
    } else {
        if cfg!(debug_assertions) {
            println!("Cyclic Mapping: starting orbit: {:?}", starting_orbit);
            println!("Cyclic Mapping: cycles to reduce: {:?}", cycles_to_reduce);
        }

        let mut mapping: mapping_ops::VertexMapping = starting_orbit.iter()
            .map(|od| vec![od[0]])
            .collect();
        let mut mapping_stack: Vec<mapping_ops::MappingStatus> = vec![];
        let mut boundary_edges: mapping_ops::BoundaryEdges = (0..mapping.len())
            .map(|_| vec![])
            .collect();

        let mut mapped_cur: usize = 0;
        while mapped_cur < mapping[0].len() {
            let mut neighbours: Vec<mapping_ops::NeighbourIndexes> = (0..mapping.len())
                .map(|idx| mapping_ops::find_new_neighbours(mapping[idx][mapped_cur], &mapping[idx], vv))
                .collect();

            update_boundary_edges(&mut mapping, &mut boundary_edges, &mut neighbours, mapped_cur, &cycles_to_reduce);

            if neighbours[0].len() > 0 {
                match mapping_ops::mapping_neighbours(&mut neighbours, numbering) {
                    Ok(neighbour_mappings) => {
                        mapping_ops::update_mapping_stack(mapped_cur, &mapping, &boundary_edges, &neighbour_mappings, &mut mapping_stack);
                    },
                    Err(_) => () 
                }

                if let Some((new_mapping_cur, new_mapping, new_boundary_edges)) = mapping_stack.pop() {
                    mapped_cur = new_mapping_cur;
                    mapping = new_mapping;
                    boundary_edges = new_boundary_edges;
                    mapped_cur += 1;
                } else {
                    return Err("Cyclic Mapping status runs out!!!".to_string());
                }
            } else {
                mapped_cur += 1;
            }

            if mapped_cur == mapping[0].len() {
                if core::orbit_ops::orbit_equal(&mapping[0], &mapping[1]) // self mapping
                || (mapping[0].len() == 1) // mapping ony 1 vertex, meaningless 
                {
                    if let Some((new_mapping_cur, new_mapping, new_boundary_edges)) = mapping_stack.pop() {
                        mapped_cur = new_mapping_cur;
                        mapping = new_mapping;
                        boundary_edges = new_boundary_edges;
                        mapped_cur += 1;
                    } else {
                        return Err("Cyclic Mapping status runs out!!!".to_string());
                    }
                }
            }
        }

        Ok((mapping, boundary_edges))
    }
}

#[cfg(test)]
mod test_reduce_case_cyclic {
    use super::*;
    use crate::ext::molecule;

    #[test]
    fn test_get_cycle_index_for_vertex() {
        type ResultType = usize;
        let params: Vec<(usize, Vec<Vec<usize>>, ResultType)> = vec![
            (5, vec![vec![1, 2], vec![3, 4], vec![5, 6], vec![7, 8]], 2),
        ];

        for param in params.iter() {
            assert_eq!(get_cycle_index_for_vertex(param.0, &param.1), param.2);
        }
    }

    #[test]
    fn test_divide_orbit_by_cycles() {
        type ParamType = (core::orbit_ops::Orbit, Vec<Vec<usize>>);
        type ReturnType = (Vec<Vec<usize>>, Vec<Vec<usize>>);
        let test_data: Vec<(ParamType, ReturnType)> = vec![
            (
                (vec![5, 18, 21, 27, 41, 47, 61, 67], vec![vec![4, 5, 6, 7, 17, 18], vec![20, 21, 22, 27, 28, 29], vec![40, 41, 42, 47, 48, 49], vec![60, 61, 62, 67, 68, 69]]),
                (vec![vec![5, 18], vec![21, 27], vec![41, 47], vec![61, 67]], vec![vec![4, 5, 6, 7, 17, 18], vec![20, 21, 22, 27, 28, 29], vec![40, 41, 42, 47, 48, 49], vec![60, 61, 62, 67, 68, 69]]),
            ),
        ];

        for td in test_data.iter() {
            let (mut orbit_divided, mut cycles_matched) = divide_orbit_by_cycles(&td.0.0, &td.0.1);
            core::orbit_ops::orbits_sort(&mut orbit_divided);
            core::orbit_ops::orbits_sort(&mut cycles_matched);
            assert_eq!(orbit_divided, td.1.0);
            assert_eq!(cycles_matched, td.1.1);
        }
    }

    #[test]
    fn test_update_boundary_edges() {
        type ParamType = (mapping_ops::VertexMapping, Vec<mapping_ops::NeighbourIndexes>, usize, Vec<Vec<usize>>);
        type ReturnType = (mapping_ops::BoundaryEdges, Vec<mapping_ops::NeighbourIndexes>);
        let test_data: Vec<(ParamType, ReturnType)> = vec![
            (
                (vec![vec![4], vec![22], vec![42], vec![62]], vec![vec![1, 5, 18], vec![21, 23, 27], vec![41, 43, 47], vec![61, 63, 67]], 0, vec![vec![4, 5, 6, 7, 17, 18], vec![20, 21, 22, 27, 28, 29], vec![40, 41, 42, 47, 48, 49], vec![60, 61, 62, 67, 68, 69]]),
                (vec![vec![(4, 1)], vec![(22, 23)], vec![(42, 43)], vec![(62, 63)]], vec![vec![5, 18], vec![21, 27], vec![41, 47], vec![61, 67]]),
            )
        ];

        for td in test_data.iter() {
            let mut boundary_edges: mapping_ops::BoundaryEdges = vec![vec![], vec![], vec![], vec![]];
            let mut neighbours: Vec<mapping_ops::NeighbourIndexes> = td.0.1.clone();
            update_boundary_edges(&td.0.0, &mut boundary_edges, &mut neighbours, td.0.2, &td.0.3);
            assert_eq!(boundary_edges, td.1.0);
            assert_eq!(neighbours, td.1.1);
        }
    }

    #[test]
    fn test_find_starting_orbit() {
        type ParamType1 = String;
        type ReturnType = (Vec<core::orbit_ops::Orbit>, Vec<Vec<usize>>);
        let test_data: Vec<(ParamType1, ReturnType)> = vec![
            (
                "O=P1([O-])OC2C3OP(=O)([O-])OP(=O)([O-])OC3C3OP(=O)([O-])OP(=O)([O-])OC3C2OP(=O)([O-])O1", // 168272
                (
                    vec![vec![1, 29], vec![7, 11], vec![18, 22]],
                    vec![vec![1, 3, 28, 29, 32], vec![6, 7, 10, 11, 14], vec![17, 18, 21, 22, 25]],
                ),
            ),
            (
                "CC(C)(C)c1cc2c(OCCCCNC(=N)N)c(c1)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)C2",
                (
                    vec![vec![4], vec![22], vec![42], vec![62]],
                    vec![vec![4, 5, 6, 7, 17, 18], vec![20, 21, 22, 27, 28, 29], vec![40, 41, 42, 47, 48, 49], vec![60, 61, 62, 67, 68, 69]],
                )
            ),
            (
                "CCn1c2ccc3cc2c2cc(ccc21)C(=O)c1ccc(cc1)Cn1cc[n+](c1)Cc1ccc(cc1)-c1cccc(c1C(=O)O)-c1ccc(cc1)C[n+]1ccn(c1)Cc1ccc(cc1)C3=O",
                (
                    vec![vec![24], vec![55]],
                    vec![vec![24, 25, 26, 27, 28], vec![52, 53, 54, 55, 56]],
                )
            )
        ].iter().map(|td| (td.0.to_string(), td.1.clone())).collect();

        for td in test_data.iter() {
            let (smiles, results) = td;
            let mol = molecule::Molecule::from_smiles(&smiles);
            println!("{}", mol.smiles_with_index(&smiles, &vec![]));
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut numbering: Vec<usize> = vec![];
            let mut orbits_after_partition: Vec<core::orbit_ops::Orbit> = vec![];
            core::givp::run::<molecule::AtomExtendable>(&vv, &mut numbering, &mut orbits_after_partition);
            let (mut starting_orbit, mut cycles) = find_starting_orbit(&orbits_after_partition, &vv, 10);
            core::orbit_ops::orbits_sort(&mut starting_orbit);
            core::orbit_ops::orbits_sort(&mut cycles);
            assert_eq!(starting_orbit, results.0);
            assert_eq!(cycles, results.1);
        }
    }

    #[test]
    fn test_mapping() {
        type ParamType1 = String;
        type ReturnType = (mapping_ops::VertexMapping, mapping_ops::BoundaryEdges);
        let test_data: Vec<(ParamType1, ReturnType)> = vec![
            (
                "CC(C)(C)c1cc2c(OCCCCNC(=N)N)c(c1)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)C2",
                (
                    vec![vec![4, 5, 6, 7, 17, 18], vec![20, 21, 22, 27, 28, 29], vec![40, 41, 42, 47, 48, 49], vec![60, 61, 62, 67, 68, 69]], 
                    vec![vec![(4, 1), (6, 79), (7, 8), (17, 19)], vec![(20, 19), (22, 23), (28, 39), (29, 30)], vec![(40, 39), (42, 43), (48, 59), (49, 50)], vec![(60, 59), (62, 63), (68, 79), (69, 70)]]
                ),
            ),
            (
                "CCn1c2ccc3cc2c2cc(ccc21)C(=O)c1ccc(cc1)Cn1cc[n+](c1)Cc1ccc(cc1)-c1cccc(c1C(=O)O)-c1ccc(cc1)C[n+]1ccn(c1)Cc1ccc(cc1)C3=O",  // chembl 15
                (
                    vec![vec![24, 25, 26, 27, 28], vec![52, 53, 54, 55, 56]], 
                    vec![vec![(24, 23), (27, 29)], vec![(52, 51), (55, 57)]]
                ),
            ),
        ].iter().map(|td| (td.0.to_string(), td.1.clone())).collect();

        for td in test_data.iter() {
            let mol = molecule::Molecule::from_smiles(&td.0);
            println!("{}", mol.smiles_with_index(&td.0, &vec![]));
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut numbering: Vec<usize> = vec![];
            let mut orbits_after_partition: Vec<core::orbit_ops::Orbit> = vec![];
            core::givp::run::<molecule::AtomExtendable>(&vv, &mut numbering, &mut orbits_after_partition);

            match create_mapping(&orbits_after_partition, &numbering, &vv, 8) {
                Ok((mut mapping, mut boundary_edges)) => { 
                    core::orbit_ops::orbits_sort(&mut mapping);
                    for be in boundary_edges.iter_mut() {
                        be.sort_unstable();
                    }
                    boundary_edges.sort_unstable();
                    assert_eq!((mapping, boundary_edges), td.1);
                },
                Err(_) => ()
            }
        }
    }

    #[test]
    fn test_molecules() {
        type ParamType1 = String;
        let test_data: Vec<ParamType1> = vec![
            //
            // *** SOLVED ***
            "O=P1([O-])OC2C3OP(=O)([O-])OP(=O)([O-])OC3C3OP(=O)([O-])OP(=O)([O-])OC3C2OP(=O)([O-])O1", // 171007
            "CCc1nc2c([nH]1)c1nc(CC)[nH]c1c1nc(CC)[nH]c21", // 189782
            "c1ccc2c(c1)c1ccccc1c1ccccc21", // 476955
            "C1CCC2C(CC1)C1CCCCCC21", // an example for the paper
            // 
            // *** UNSOLVED *** 
            // "OCCN=C1c2ccccc2C2C(c3ccccc31)C1c3ccccc3C(=NCCO)c3ccccc3C21", // 1295466
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

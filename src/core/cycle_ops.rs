// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Cycle related operations
//! 
use super::orbit_ops;
use super::graph;

/// Extend on edge in the graph, store all the cyclic routes and the open routes
fn extend_one_edge<T: graph::Vertex>(
    open_route: &Vec<usize>,
    vv: &graph::VertexVec<T>,
    cyclic_routes: &mut Vec<Vec<usize>>,
    open_routes: &mut Vec<Vec<usize>>
) {
    let length = open_route.len();
    for ni in vv[open_route[length - 1]].neighbour_indexes() {
        if open_route.contains(&ni) {
            if (ni == open_route[0]) && (open_route.len() > 2) {
                cyclic_routes.push(open_route.clone());
            }
        } else {
            let mut indexes_tmp = open_route.clone();
            indexes_tmp.push(ni);
            open_routes.push(indexes_tmp);
        }
    }
}

/// Find the size of the minimun cycle which the vertex belongs to. If the vertex is acylic, return 0 
pub fn cycle_size<T: graph::Vertex>(
    vertex_index: usize,
    vv: &graph::VertexVec<T>
) -> usize {
    let mut open_routes: Vec<Vec<usize>> = vec![vec![vertex_index]];
    let mut cyclic_routes: Vec<Vec<usize>> = vec![];

    while open_routes.len() > 0 {
        let open_routes_tmp = open_routes.clone();
        open_routes.clear();
        for route in open_routes_tmp.iter() {
            extend_one_edge(route, vv, &mut cyclic_routes, &mut open_routes);
            if cyclic_routes.len() > 0 {
                return cyclic_routes[0].len()
            }
        }
    }

    0
}

/// Find cycles that the two specified vertices belong to.
pub fn find_cycle_for_vertices<T: graph::Vertex>(
    vertice_idx_1: usize,
    vertice_idx_2: usize,
    vv: &graph::VertexVec<T>,
    max_steps: usize
) -> Option<Vec<usize>> {
    let mut open_routes: Vec<Vec<usize>> = vec![vec![vertice_idx_1]];
    let mut cyclic_routes: Vec<Vec<usize>> = vec![];
    for _ in 0..max_steps {
        let open_routes_tmp = open_routes.clone();
        open_routes.clear();
        for route in open_routes_tmp.iter() {
            extend_one_edge(route, vv, &mut cyclic_routes, &mut open_routes);
        }
    } 

    cyclic_routes.sort_by_key(|r| r.len());
    for route in cyclic_routes.iter() {
        if route.contains(&vertice_idx_2) {
            return Some(route.clone())
        }
    }

    return None
}

/// Find vertex cycles that include the specified vertex, with a maximum size. If one cycle is found, return
fn find_cycles<T: graph::Vertex>(
    vertex_index: usize,
    vv: &graph::VertexVec<T>,
    max_cycle_size: usize
) -> Vec<Vec<usize>> {
    let mut open_routes: Vec<Vec<usize>> = vec![vec![vertex_index]];
    let mut cyclic_routes: Vec<Vec<usize>> = vec![];

    while open_routes.len() > 0 && open_routes[0].len() <= max_cycle_size {
        let open_routes_tmp = open_routes.clone();
        open_routes.clear();
        for route in open_routes_tmp.iter() {
            extend_one_edge(route, vv, &mut cyclic_routes, &mut open_routes);
            if cyclic_routes.len() > 0 {
                return cyclic_routes
            }
        }
    }

    vec![]
}

/// Find all the cycles that include the specified vertex.
pub fn find_all_cycles<T: graph::Vertex>(
    vertex_index: usize,
    vv: &graph::VertexVec<T>,
    max_cycle_size: usize
) -> Vec<Vec<usize>> {
    let mut open_routes: Vec<Vec<usize>> = vec![vec![vertex_index]];
    let mut cyclic_routes: Vec<Vec<usize>> = vec![];

    while open_routes.len() > 0 && open_routes[0].len() <= max_cycle_size {
        let open_routes_tmp = open_routes.clone();
        open_routes.clear();
        for route in open_routes_tmp.iter() {
            extend_one_edge(route, vv, &mut cyclic_routes, &mut open_routes);
        }
    }

    orbit_ops::orbits_sort(&mut cyclic_routes);
    cyclic_routes.dedup();
    cyclic_routes
}

/// Find related cycles that include any vertex in the orbits
pub fn get_cycles_from_orbits<T: graph::Vertex>(
    orbits: &Vec<orbit_ops::Orbit>,
    vv: &graph::VertexVec<T>,
    max_cycle_size: usize,
) -> Vec<Vec<usize>> {
    let mut cycles: Vec<Vec<usize>> = vec![];

    for orbit in orbits.iter() {
        for &idx in orbit.iter() {
            cycles.append(&mut find_cycles(idx, vv, max_cycle_size));
        }
    }

    cycles
}

/// Merge cycles that share any vertex in a given set of cycles.
pub fn merge_cycles(
    cycles: &Vec<Vec<usize>>
) -> Vec<Vec<usize>> {
    let mut isolated_cycles: Vec<Vec<usize>> = vec![];
    let mut cycles_cloned = cycles.to_vec();

    while cycles_cloned.len() > 0 {
        if let Some(mut cycle) = cycles_cloned.pop() {
            let mut is_isolated: bool = true;
            for idx in 0..cycles_cloned.len() {
                if orbit_ops::orbit_overlap(&cycle, &cycles_cloned[idx]) {
                    is_isolated = false;
                    cycles_cloned[idx].append(&mut cycle);
                    cycles_cloned[idx].sort_unstable();
                    cycles_cloned[idx].dedup();
                    break;
                }
            }

            if is_isolated {
                isolated_cycles.push(cycle);
            }
        }
    }

    isolated_cycles
}

#[cfg(test)]
mod test_core_cycle_ops {
    use super::*;
    use crate::ext::molecule;

    #[test]
    fn test_extend_one_edge() {
        let smiles: String = "c1ccccc1CN".to_string();
        let mol = molecule::Molecule::from_smiles(&smiles);
        let vv = graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());

        let mut cyclic_routes: Vec<Vec<usize>> = vec![];
        let mut open_routes: Vec<Vec<usize>> = vec![];
        extend_one_edge(&vec![6], &vv, &mut cyclic_routes, &mut open_routes);
        assert_eq!(open_routes, vec![vec![6, 5], vec![6, 7]]);
        open_routes.clear();
        cyclic_routes.clear();
        extend_one_edge(&vec![6, 5], &vv, &mut cyclic_routes, &mut open_routes);
        assert_eq!(open_routes, vec![vec![6, 5, 4], vec![6, 5, 0]]);
    }

    #[test]
    fn test_cyclic_size_find_cycles() {
        let smiles: String = "NCCCNCCCCN(CCCN)C(=O)CCCCNC(=O)c1ccc(-c2c3nc(c(-c4ccc(C(=O)NCCCCC(=O)N(CCCN)CCCCNCCCN)cc4)c4ccc([nH]4)c(-c4ccc(C(=O)NCCCCC(=O)N(CCCN)CCCCNCCCN)cc4)c4nc(c(-c5ccc(C(=O)NCCCCC(=O)N(CCCN)CCCCNCCCN)cc5)c5ccc2[nH]5)C=C4)C=C3)cc1".to_string();
        let mol = molecule::Molecule::from_smiles(&smiles);
        let vv = graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());

        assert_eq!(cycle_size(65, &vv), 5); 
        assert_eq!(find_cycles(65, &vv, 4).len(), 0); 
        assert_eq!(find_cycles(65, &vv, 5)[0].len(), 5); 
        assert_eq!(find_cycles(65, &vv, 6)[0].len(), 5); 
        assert_eq!(cycle_size(31, &vv), 16); 
        assert_eq!(cycle_size(1, &vv), 0); 
    }

    #[test]
    fn tet_find_cycle_for_vertices() {
        let smiles: String = "c1ccccc1CN".to_string();
        let mol = molecule::Molecule::from_smiles(&smiles);
        let vv = graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());

        assert_eq!(find_cycle_for_vertices(0, 6, &vv, 6), None);
        assert_eq!(find_cycle_for_vertices(0, 5, &vv, 6), Some(vec![0, 5, 4, 3, 2, 1]));
        assert_eq!(find_cycle_for_vertices(0, 5, &vv, 5), None);
    }

    #[test]
    fn test_find_all_cycles() {
        type ParamType1 = String;
        type ParamType2 = usize;
        type ReturnType = usize;
        let test_data: Vec<(ParamType1, ParamType2, ReturnType)> = vec![
            ("c1ccc2cc3ccccc3cc2c1", 1, 1),
            ("c1ccc2cc3ccccc3cc2c1", 3, 2),
            ("c1ccc2cc3ccccc3cc2c1", 11, 1),
            ("CCOC12c3c4cccc3Oc3cccc(c31)Oc1cccc(c12)O4", 3, 3),
            ("CCOC12c3c4cccc3Oc3cccc(c31)Oc1cccc(c12)O4", 4, 3),
            ("CCOC12c3c4cccc3Oc3cccc(c31)Oc1cccc(c12)O4", 16, 3),
        ].iter().map(|td| (td.0.to_string(), td.1.clone(), td.2.clone())).collect();

        for td in test_data.iter() {
            let (smiles, vertex_index, cycle_count) = td;
            let mol = molecule::Molecule::from_smiles(&smiles);
            println!("{}", mol.smiles_with_index(&smiles, &vec![]));
            let vv = graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            assert_eq!(find_all_cycles(*vertex_index, &vv, 6).len(), *cycle_count);
        }
    }

    #[test]
    fn test_get_cycles_from_orbits() {
        let smiles_vec: Vec<String> = vec![
            "CC(C)(C)c1cc2c(OCCCCNC(=N)N)c(c1)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)C2".to_string()
        ].iter().map(|s| s.to_string()).collect();
        let params: Vec<(Vec<orbit_ops::Orbit>, Vec<orbit_ops::Orbit>)> = vec![
            (
                vec![vec![4, 22], vec![7, 28], vec![79, 39]],
                vec![vec![4, 5, 6, 7, 17, 18], vec![4, 5, 6, 7, 17, 18], vec![20, 21, 22, 27, 28, 29], vec![20, 21, 22, 27, 28, 29]]
            )
        ];


        for idx in 0..smiles_vec.len() {
            let mol = molecule::Molecule::from_smiles(&smiles_vec[idx]);
            if cfg!(debug_assertions) {
                println!("{}", mol.smiles_with_index(&smiles_vec[idx], &vec![]));
            }
            let vv = graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut cycles = get_cycles_from_orbits(&params[idx].0, &vv, 6);
            orbit_ops::orbits_sort(&mut cycles);
            assert_eq!(cycles, params[idx].1);
        }
    }

    #[test]
    fn test_merge_cycles() {
        let params: Vec<(Vec<Vec<usize>>, Vec<Vec<usize>>)> = vec![
            (vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8]], vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8]]),
            (vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8], vec![1, 9, 10]], vec![vec![1, 2, 3, 4, 9, 10], vec![5, 6, 7, 8]]),
            (vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8], vec![1, 9, 10], vec![5, 11, 12], vec![1, 5]], vec![vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]), 
        ];

        for param in params.iter() {
            let mut cycles_merged = merge_cycles(&param.0);
            orbit_ops::orbits_sort(&mut cycles_merged);
            assert_eq!(cycles_merged, param.1);
        }
    
    }

}
 
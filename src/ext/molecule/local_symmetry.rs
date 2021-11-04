// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::core;
use crate::core::graph::*;
use super::atom;

/// Seperate vertices in an orbit according to their neighbours. Only applicable for endpoint atoms
fn group_by_neighbour(
    orbit: &core::orbit_ops::Orbit,
    vv: &core::graph::VertexVec<atom::Atom>,
) -> Vec<core::orbit_ops::Orbit> {
    let mut atom_groups: std::collections::HashMap<usize, core::orbit_ops::Orbit> = std::collections::HashMap::new();
    for &ai in orbit.iter() {
        let nb_idx: usize = vv[ai].neighbour_indexes()[0]; // endpoint atom has only one neighbour
        atom_groups.entry(nb_idx).or_insert(vec![]).push(ai);
    }

    atom_groups.into_values().collect()
}

/// Acyclic Local Symmetry
fn is_acyclic_local_symmetric(
    orbit: &core::orbit_ops::Orbit,
    vv: &core::graph::VertexVec<atom::Atom>,
) -> bool {
    if vv[orbit[0]].degree() != 1 {
        return false;
    }

    let orbits_local_symmetry: Vec<core::orbit_ops::Orbit> = group_by_neighbour(&orbit, vv);
    orbits_local_symmetry.len() == 1
}

/// Find the cycle including the vertex with minimum size
fn find_cyclic_route(
    vertex: usize,
    vv: &core::graph::VertexVec<atom::Atom>,
    max_size: usize
) -> Vec<usize> {
    let mut cyclic_route: Vec<usize> = vec![0; max_size + 1];
    for &neighbour in vv[vertex].neighbour_indexes().iter() {
        match core::cycle_ops::find_cycle_for_vertices(vertex, neighbour, vv, max_size) {
            Some(cr) => {
                if cr.len() < cyclic_route.len() {
                    cyclic_route = cr
                }
            },
            None => ()
        }
    }

    cyclic_route
}

/// Cyclic Local Symmetry
fn is_cyclic_local_symmetric(
    orbit: &core::orbit_ops::Orbit,
    vv: &core::graph::VertexVec<atom::Atom>,
    max_size: usize
) -> bool {
    if orbit.len() != 2 {
        return false;
    }

    if vv[orbit[0]].degree() != 2 {
        return false;
    }

    // find the minimum cycle for orbit
    let cyclic_route: Vec<usize> = find_cyclic_route(orbit[0], vv, max_size); 
    let cyclic_route_1: Vec<usize> = find_cyclic_route(orbit[1], vv, max_size); 
    if !core::orbit_ops::orbit_equal(&cyclic_route, &cyclic_route_1) {
        return false;
    }

    if cyclic_route.len() > max_size {
        return false;
    }

    let boundary_vertices: Vec<usize> = cyclic_route.clone().into_iter()
        .filter(|&vi| vv[vi].degree() > 2)
        .collect();

    if cyclic_route.len() % 2 == 1 {
        boundary_vertices.len() == 1
    } else {
        (boundary_vertices.len() == 1 && core::graph_ops::graph_distance(boundary_vertices[0], orbit[0], vv) != cyclic_route.len() / 2)
        || (boundary_vertices.len() == 2 && core::graph_ops::graph_distance(boundary_vertices[0], boundary_vertices[1], vv) == cyclic_route.len() / 2 && core::graph_ops::graph_distance(boundary_vertices[0], orbit[0], vv) == core::graph_ops::graph_distance(boundary_vertices[0], orbit[1], vv))
    }
}

pub fn get_local_symmetric_orbits(
    vv: &core::graph::VertexVec<atom::Atom>,
    orbits_residual: &mut Vec<core::orbit_ops::Orbit>,
    orbits_symmetry: &mut Vec<core::orbit_ops::Orbit>
) {
    let mut orbits_tmp = orbits_residual.clone();
    orbits_residual.clear();

    while let Some(orbit) = orbits_tmp.pop() {
        if is_acyclic_local_symmetric(&orbit, vv) {
            orbits_symmetry.push(orbit);
        } else if is_cyclic_local_symmetric(&orbit, vv, 6) {
            orbits_symmetry.push(orbit);
        } else {
            orbits_residual.push(orbit);
        }
    }
}


#[cfg(test)]
mod test_ext_mol_local_symmetry {
    use super::*;
    use super::super::molecule;


    #[test]
    fn test_group_by_neighbour() {
        let test_data: Vec<(String, core::orbit_ops::Orbit, Vec<core::orbit_ops::Orbit>)> = vec![
            ("C(C)(C)CCN",
            vec![1, 2], vec![vec![1, 2]]),
            ("C(C)(C)CCNCCC(C)(C)",
            vec![1, 2, 9, 10], vec![vec![1, 2], vec![9, 10]])
        ].into_iter().map(|s| (s.0.to_string(), s.1, s.2)).collect();

        for td in test_data.iter() {
            let mol = molecule::Molecule::from_smiles(&td.0);
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            assert_eq!(core::orbit_ops::orbits_equal(&group_by_neighbour(&td.1, &vv), &td.2), true);
        }
    }

    #[test]
    fn test_is_acyclic_local_symmetric() {
        type InputType1 = String;
        type InputType2 = core::orbit_ops::Orbit;
        type ReturnType = bool;
        let test_data: Vec<(InputType1, InputType2, ReturnType)> = vec![
            (
                "C(C)(C)CCN",
                vec![1, 2], true 
            ),
            (
                "C(C)(C)CCNCCC(C)(C)",
                vec![1, 2, 9, 10], false 
            )
        ].into_iter().map(|s| (s.0.to_string(), s.1, s.2)).collect();

        for td in test_data.iter() {
            let mol = molecule::Molecule::from_smiles(&td.0);
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            assert_eq!(is_acyclic_local_symmetric(&td.1, &vv), td.2);
        }
    }

    #[test]
    fn test_is_cyclic_local_symmetric() {
        type InputType1 = String;
        type InputType2 = core::orbit_ops::Orbit;
        type ResultType = bool; 
        let test_data: Vec<(InputType1, Vec<(InputType2, ResultType)>)> = vec![
            (
                "C1C(N)C1",
                vec![
                    (vec![0, 3], true) 
                ]
            ),
            (
                "OC1CCC1",
                vec![
                    (vec![2, 4], true)
                ]
            ),
            (
                "OC1C(N)CC1",
                vec![
                    (vec![4, 5], false)
                ]
            ),
            (
                "OC1CC(N)C1",
                vec![
                    (vec![2, 5], true)
                ]
            ),
            (
                "OC1CCCC1",
                vec![
                    (vec![2, 5], true),
                    (vec![3, 4], true),
                ]
            ),
            (
                "OC1C(C)CCC1",
                vec![
                    (vec![4, 5], false)
                ]
            ),
            (
                "OC1CC(C)CC1",
                vec![
                    (vec![5, 6], false)
                ]
            ),
            (
                "Oc1ccccc1",
                vec![
                    (vec![2, 6], true),
                    (vec![3, 5], true),
                ]
            ),
            (
                "Oc1c(N)cccc1",
                vec![
                    (vec![5, 6], false),
                ]
            ),
            (
                "Oc1cc(N)ccc1",
                vec![
                    (vec![5, 6], false),
                ]
            ),
            (
                "Oc1ccc(N)cc1",
                vec![
                    (vec![2, 7], true),
                    (vec![3, 6], true),
                ]
            ),
        ].into_iter().map(|s| (s.0.to_string(), s.1)).collect();

        for td in test_data.iter() {
            let mol = molecule::Molecule::from_smiles(&td.0);
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            for test_param in td.1.iter() {
                assert_eq!(is_cyclic_local_symmetric(&test_param.0, &vv, 6), test_param.1);
            }
        }
    }

    #[test]
    fn test_get_local_symmetric_orbits() {
        type InputType1 = String;
        type InputType2 = Vec<core::orbit_ops::Orbit>;
        type ResultType = (Vec<core::orbit_ops::Orbit>, Vec<core::orbit_ops::Orbit>);
        let test_data: Vec<(InputType1, InputType2, ResultType)> = vec![
            (
                "FC(F)(F)c1ccc(C[N+]23CC[C@]45c6ccccc6N6[C@H]4[C@H]4[C@@H](C[C@@H]52)C(=CCO[C@H]4N2c4ccccc4[C@@]45CC[N+]7(Cc8ccc(C(F)(F)F)cc8)CC8=CCO[C@@H]6[C@@H]([C@H]24)[C@H]8C[C@@H]57)C3)cc1",
                vec![vec![6, 43, 51, 64]],
                (vec![vec![6, 43, 51, 64]], vec![]),
            ),
            (
                "FC(F)(F)c1ccc(C[N+]23CC[C@]45c6ccccc6N6[C@H]4[C@H]4[C@@H](C[C@@H]52)C(=CCO[C@H]4N2c4ccccc4[C@@]45CC[N+]7(Cc8ccc(C(F)(F)F)cc8)CC8=CCO[C@@H]6[C@@H]([C@H]24)[C@H]8C[C@@H]57)C3)cc1",
                vec![vec![0, 2, 3, 47, 48, 49]],
                (vec![vec![0, 2, 3, 47, 48, 49]], vec![]),
            ),
            (
                r#"COc1cc(Cc2cnc(/N=C3\C(=O)N(CN(Cc4ccccc4)Cc4ccccc4)c4ccc(Cl)cc43)nc2N)cc(OC)c1OC"#,
                vec![vec![0, 44], vec![1, 43], vec![2, 42], vec![3, 41], vec![17, 24], vec![18, 25], vec![19, 23, 26, 30], vec![20, 22, 27, 29], vec![21, 28]],
                (vec![vec![0, 44], vec![1, 43], vec![2, 42], vec![3, 41], vec![17, 24], vec![18, 25], vec![19, 23, 26, 30], vec![20, 22, 27, 29], vec![21, 28]], vec![])
            ),
        ].into_iter().map(|s| (s.0.to_string(), s.1, s.2)).collect();
        for td in test_data.iter() {
            let (smiles, mut orbits_residual, result) = td.clone();
            let mol = molecule::Molecule::from_smiles(&smiles);
            // println!("{}", mol.smiles_with_index(&smiles));
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut orbits_symmetry: Vec<core::orbit_ops::Orbit> = vec![];
            get_local_symmetric_orbits(&vv, &mut orbits_residual, &mut orbits_symmetry);
            assert_eq!(core::orbit_ops::orbits_equal(&orbits_residual, &result.0), true); 
            assert_eq!(core::orbit_ops::orbits_equal(&orbits_symmetry, &result.1), true); 
        }
    }

}
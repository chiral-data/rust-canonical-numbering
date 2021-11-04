// Copyr1ight 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Separable Graph


use crate::core;
use super::reducible_graph;

fn are_cycles_overlapping(
    cycles: &Vec<Vec<usize>>
) -> bool {
    for i in 0..cycles.len() {
        for j in (i+1)..cycles.len() {
            if core::orbit_ops::orbit_overlap(&cycles[i], &cycles[j]) {
                return true
            }
        }
    }

    false
}

fn is_orbit_in_cycle(
    orbit: &core::orbit_ops::Orbit,
    cycles: &Vec<Vec<usize>>
) -> bool {
    for cycle in cycles.iter() {
        if core::orbit_ops::orbit_cover(cycle, orbit) {
            return true;
        }
    }

    false
}

fn are_orbits_seprated_in_cycles(
    orbits: &Vec<core::orbit_ops::Orbit>,
    cycles: &Vec<Vec<usize>>
) -> bool {
    orbits.iter()
        .map(|orbit| is_orbit_in_cycle(orbit, cycles))
        .fold(true, |acc, x| acc & x)
}

pub fn find_cycles<T: core::graph::Vertex>(
    rg: &reducible_graph::ReducibleGraph<T>,
    max_cycle_size: usize,
) -> Vec<Vec<usize>> {
    core::cycle_ops::merge_cycles(&core::cycle_ops::get_cycles_from_orbits(&rg.orbits_after_partition, &rg.vv, max_cycle_size))
}

pub fn is_separable<T: core::graph::Vertex>(
    rg: &reducible_graph::ReducibleGraph<T>,
    cycles: &Vec<Vec<usize>>
) -> bool {
    if cfg!(debug_assertions) {
        let mut cycles_to_print = cycles.to_vec();
        core::orbit_ops::orbits_sort(&mut cycles_to_print);
        println!("Seperable Graphs: Cycles: \n{:?}", cycles_to_print);
    }
    cycles.len() > 1 && !are_cycles_overlapping(cycles) && are_orbits_seprated_in_cycles(&rg.orbits_after_partition, cycles)
}

fn find_boundary_vertices<T: core::graph::Vertex>(
    cycle: &Vec<usize>,
    vv: &core::graph::VertexVec<T>,
) -> Vec<usize> {
    let mut boundary_vertices: Vec<usize> = cycle.iter()
        .map(|&vi| vv[vi].neighbour_indexes())
        .flatten()
        .filter(|vi| !cycle.contains(vi))
        .collect();
    
    boundary_vertices.sort_unstable();
    boundary_vertices.dedup();
    boundary_vertices
}

pub fn construct_reducible_graph<T: core::graph::Vertex>(
    rg: &reducible_graph::ReducibleGraph<T>,
    cycle: &Vec<usize>,
    custom_maker: &mut usize,
) -> reducible_graph::ReducibleGraph<T> {
    let mut indexes_new = cycle.to_vec();
    let mut vertices_new = rg.vv.all_vertices().to_vec();

    let boundary_vertices = find_boundary_vertices(cycle, &rg.vv);
    let mut bv_numberings: Vec<usize> = boundary_vertices.iter().map(|&bv| rg.numbering[bv]).collect();
    bv_numberings.sort_unstable();
    bv_numberings.dedup();

    for &bv in boundary_vertices.iter() {
        let atomic_number_adding: usize = bv_numberings.iter().position(|&r| r == rg.numbering[bv]).unwrap() + 1;
        let custom_vertice = T::custom_new_in_separated_graph(*custom_maker + atomic_number_adding, bv, cycle, rg.vv.all_vertices());
        vertices_new[bv] = custom_vertice;
        indexes_new.push(bv);
    }
    *custom_maker += bv_numberings.len();

    let new_vv = core::graph::VertexVec::init(indexes_new, vertices_new);
    reducible_graph::ReducibleGraph {
        vv: new_vv,
        mapping: vec![],
        boundary_edges: vec![],
        orbits_after_partition: vec![],
        numbering: vec![] 
    }
}

#[cfg(test)]
mod test_reduce_graph_separable {
    use crate::ext::molecule;
    use super::*;

    #[test]
    fn test_is_seperable() {
        type InputType1 = String;
        type ReturnType = bool; 
        let test_data: Vec<(InputType1, ReturnType)> = vec![
            ("C1C2CC3CC1CC(C2)c1nc2c(nc13)C1CC3CC(C1)CC2C3", false), // 2064956
        ].iter().map(|td| (td.0.to_string(), td.1.clone())).collect();

        for td in test_data.iter() {
            let (smiles, result) = td; 
            let mol = molecule::Molecule::from_smiles(&smiles);
            println!("{}", mol.smiles_with_index(&smiles, &vec![]));
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut rg = core::reduce::reducible_graph::ReducibleGraph {
                vv: vv,
                mapping: vec![],
                boundary_edges: vec![],
                orbits_after_partition: vec![], 
                numbering: vec![] 
            };
            core::givp::run::<molecule::AtomExtendable>(&rg.vv, &mut rg.numbering, &mut rg.orbits_after_partition);
            let cycles = find_cycles(&rg, 6);

            assert_eq!(is_separable(&rg, &cycles), *result);
        }
    }

    #[test]
    fn test_construct_reducible_graph() {
        type ParamType1 = String;
        type ParamType2 = (Vec<Vec<usize>>, Vec<Vec<usize>>);
        type ReturnType = Vec<Vec<usize>>;
        let test_data: Vec<(ParamType1, ParamType2, ReturnType)> = vec![
            (
                "C[N+](C)(CCCCCC[N+](C)(C)CCCN1C(=O)C2C3c4ccccc4C(c4ccccc43)C2C1=O)CCCN1C(=O)c2ccccc2C1=O",
                (vec![], vec![]),
                vec![vec![38, 41, 49], vec![14, 17, 35]]
            ),
        ].iter().map(|td| (td.0.to_string(), td.1.clone(), td.2.clone())).collect();

        for td in test_data.iter() {
            let (smiles, _, results) = td;
            let mol = molecule::Molecule::from_smiles(&smiles);
            println!("{}", mol.smiles_with_index(&smiles, &vec![]));
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut rg = core::reduce::reducible_graph::ReducibleGraph {
                vv: vv,
                mapping: vec![],
                boundary_edges: vec![],
                orbits_after_partition: vec![], 
                numbering: vec![] 
            };
            core::givp::run::<molecule::AtomExtendable>(&rg.vv, &mut rg.numbering, &mut rg.orbits_after_partition);

            let cycles = find_cycles(&rg, 6);
            for (idx, cycle) in cycles.iter().enumerate() {
                let new_rg = construct_reducible_graph(&rg, cycle, &mut 200);
                let mut new_indexes = cycle.to_vec();
                new_indexes.append(&mut results[idx].clone());
                assert_eq!(core::orbit_ops::orbit_equal(new_rg.vv.valid_indexes(), &new_indexes), true); 
            }
        }
    }

    #[test]
    fn test_molecules() {
        type ParamType1 = String;
        let test_data: Vec<ParamType1> = vec![
            "S=C(Nc1ccc(C23CC4CC(CC(C4)C2)C3)cc1)NC12CC3CC(CC(C3)C1)C2", // 2014622
        ].iter().map(|td| td.to_string())
        .collect();

        for td in test_data.iter() {
            let smiles = td.clone();
            let mol = molecule::Molecule::from_smiles(&smiles);
            if cfg!(debug_assertions) {
                println!("{}", mol.smiles_with_index(&smiles, &vec![]));
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
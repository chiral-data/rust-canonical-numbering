// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::core;
use crate::core::graph::*;
use super::partition;

/// Quicksort by hash_values, that could be any type with trait PartialOrd implemented
pub fn partition_by_hash_values<T: std::cmp::PartialOrd>(
    orbit: &core::orbit_ops::Orbit,
    hash_values: &Vec<T>,
    current_numbering: usize,
    numbering: &mut Vec<usize>,
    residual_orbits: &mut Vec<core::orbit_ops::Orbit>
){
    partition::partition_recursively(orbit, hash_values, current_numbering, numbering, residual_orbits);
}

/// Quicksort and Extension algorithm for extendable hash_values 
pub fn partition_by_extendable_hash_values<T: core::graph::VertexExtendableHash>(
    orbits_to_be_partitioned: &Vec<core::orbit_ops::Orbit>,
    vertices: &Vec<T::VertexType>,
    fixed_hash_values: &Vec<core::graph::VertexFixedHashValue>,
    extendable_hash_entities: &mut Vec<T>,
    extendable_hash_values: &mut Vec<core::graph::VertexExtendableHashValue>,
    numbering: &mut Vec<usize>,
    orbits_residual: &mut Vec<core::orbit_ops::Orbit>
){
    let mut orbits_to_be_partitioned_tmp = orbits_to_be_partitioned.to_vec();
    while orbits_to_be_partitioned_tmp.len() > 0 {
        if let Some(orbit_tbp) = orbits_to_be_partitioned_tmp.pop() {
            if extendable_hash_entities[orbit_tbp[0]].is_completed() {
                orbits_residual.push(orbit_tbp);
            } else {
                for &idx in orbit_tbp.iter() {
                    extendable_hash_entities[idx].extend(vertices);
                    extendable_hash_values[idx] = extendable_hash_entities[idx].value(fixed_hash_values);
                }
                partition_by_hash_values::<core::graph::VertexExtendableHashValue>(&orbit_tbp, extendable_hash_values, numbering[orbit_tbp[0]], numbering, &mut orbits_to_be_partitioned_tmp);
            }
        }
    }
}

fn partition_vertices_with_fixed_hash_value<T: core::graph::VertexExtendableHash>(
    vertice_indexes: &Vec<usize>,
    vertices: &Vec<T::VertexType>,
    fixed_hash_values: &Vec<core::graph::VertexFixedHashValue>,
    numbering: &mut Vec<usize>,
    orbits_output: &mut Vec<core::orbit_ops::Orbit>
) {
    orbits_output.clear();
    partition_by_hash_values::<core::graph::VertexFixedHashValue>(vertice_indexes, &fixed_hash_values, vertices.len(), numbering, orbits_output); 
}

fn partition_vertices_with_extendable_hash_value<T: core::graph::VertexExtendableHash>(
    vertices: &Vec<T::VertexType>,
    fixed_hash_values: &Vec<core::graph::VertexFixedHashValue>,
    orbits_to_be_partitioned: &Vec<core::orbit_ops::Orbit>,
    numbering: &mut Vec<usize>,
    orbits_output: &mut Vec<core::orbit_ops::Orbit>
) {
    let mut extendable_hash_entities: Vec<T> = (0..vertices.len())
        .map(|idx| T::init(idx))
        .collect();
    let mut extendable_hash_values: Vec<core::graph::VertexExtendableHashValue> = extendable_hash_entities.iter()
        .map(|ent| ent.value(&fixed_hash_values))
        .collect();
    orbits_output.clear();
    partition_by_extendable_hash_values(
        orbits_to_be_partitioned, vertices, &fixed_hash_values, &mut extendable_hash_entities, &mut extendable_hash_values, numbering, orbits_output);
}

/// Partition vertices with fixed and extendable hash_values 
fn partition_vertices_once<T: core::graph::VertexExtendableHash>(
    vertice_indexes: &Vec<usize>,
    vertices: &Vec<T::VertexType>,
    fixed_hash_values: &Vec<core::graph::VertexFixedHashValue>,
    numbering: &mut Vec<usize>,
    orbits_output: &mut Vec<core::orbit_ops::Orbit>
) {
    partition_vertices_with_fixed_hash_value::<T>(vertice_indexes, vertices, fixed_hash_values, numbering, orbits_output);
    let orbits_to_be_partitioned: Vec<core::orbit_ops::Orbit> = orbits_output.clone();
    orbits_output.clear();
    partition_vertices_with_extendable_hash_value::<T>(vertices, fixed_hash_values, &orbits_to_be_partitioned, numbering, orbits_output);
}

pub fn partition_vertices<T: core::graph::VertexExtendableHash>(
    vertice_indexes: &Vec<usize>,
    vertices: &Vec<T::VertexType>
) -> Vec<core::orbit_ops::Orbit> {
    let fixed_hash_values: Vec<core::graph::VertexFixedHashValue> = vertices.iter()
        .map(|v| v.fixed_hash_value())
        .collect();
    let mut numbering: Vec<usize> = vec![vertices.len(); vertices.len()];
    let mut orbits_output: Vec<core::orbit_ops::Orbit> = vec![];
    partition_vertices_once::<T>(&vertice_indexes, vertices, &fixed_hash_values, &mut numbering, &mut orbits_output);
    
    orbits_output
}

/// GIVP procese
///     if vertex numbering is given, it will partition according to the extendable hash_values directly
pub fn run<T: core::graph::VertexExtendableHash>(
    vv: &core::graph::VertexVec<T::VertexType>,
    numbering: &mut Vec<usize>,
    orbits_output: &mut Vec<core::orbit_ops::Orbit>
) {
    if numbering.len() == 0 {
        let fixed_hash_values: Vec<core::graph::VertexFixedHashValue> = vv.all_vertices().iter()
            .map(|v| v.fixed_hash_value())
            .collect();
        *numbering = vec![vv.len(); vv.all_len()];
        orbits_output.clear();
        partition_vertices_once::<T>(vv.valid_indexes(), vv.all_vertices(), &fixed_hash_values, numbering, orbits_output);
    }

    // do GIAP again
    // case 469: COc1ccc(C(=O)C[n+]2c(C)n(Cc3c4c(cc5c3OCC5)OCC4)c3ccccc32)cc1
    let orbits_to_be_partitioned = orbits_output.clone();
    orbits_output.clear();
    partition_vertices_with_extendable_hash_value::<T>(vv.all_vertices(), &numbering.clone(), &orbits_to_be_partitioned, numbering, orbits_output);
}

#[cfg(test)]
mod test_givp_workflow {
    use super::*;
    use crate::core;
    use crate::ext::molecule;

    #[test]
    fn test_partition_molecule() {
        type InputType1 = String;
        type ReturnType = Vec<core::orbit_ops::Orbit>;
        let test_data: Vec<(InputType1, ReturnType)> = vec![
            ("c1ccccc1CN", vec![vec![0, 4], vec![1, 3]]),
            (
                "COc1ccc(C(=O)C[n+]2c(C)n(Cc3c4c(cc5c3OCC5)OCC4)c3ccccc32)cc1",
                vec![vec![3, 33], vec![4, 32]],  // chembl 469, if no 2nd round givp, there will be another orbit [21, 24]
            ),
            (
                r#"Cc1[nH]c2ccccc2c1/C=C1\SC(=N)N(c2nccs2)C1=O"#, // 926178 
                vec![],
            ),
            ("C(C)(C)CCNCCC(C)(C)", vec![vec![0, 8], vec![1, 2, 9, 10], vec![3, 7], vec![4, 6]]) // bug Issue #1
        ].into_iter().map(|td| (td.0.to_string(), td.1)).collect();


        for td in test_data.iter() {
            let (smiles, results) = td;
            let mol = molecule::Molecule::from_smiles(smiles);
            println!("{}", mol.smiles_with_index(smiles, &vec![]));
            let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
            let mut orbits: Vec<core::orbit_ops::Orbit> = vec![];
            let mut numbering: Vec<usize> = vec![];
            run::<molecule::AtomExtendable>(&vv, &mut numbering, &mut orbits);
            core::orbit_ops::orbits_sort(&mut orbits);
            assert_eq!(orbits, results.to_vec()); 
        }
    }
}
// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::core;
use super::combinatorial;

/// Types 

/// Permutation is a rearrangement of vertex order, then can be converted to a mapping
///     eg. permutation [1, 3, 2, 0] means vertex mapping (0 -> 1), (1 -> 3), (2 -> 2), (3 -> 0)
pub type Permuation = Vec<usize>;
/// Mapping stores how two elements from two sets pair
///     refer to the definition of Bijection: https://en.wikipedia.org/wiki/Bijection
pub type Mapping = std::collections::HashMap<usize, usize>;

/// Insert new elements into the mapping
fn update_mapping(
    vertex_from_indexes: &Vec<usize>,
    vertex_to_indexes: &Vec<usize>,
    mapping: &mut Mapping
) {
    for i in 0..vertex_from_indexes.len() {
        mapping.insert(vertex_from_indexes[i], vertex_to_indexes[i]);
    }
}

/// Conversion from Mapping to Permutation
fn permutation_from_mapping(
    mapping: &Mapping,
    length: usize
) -> Permuation {
    let mut permutation_vec: Permuation = vec![0; length];

    for idx in 0..length {
        match mapping.get(&idx) {
            Some(target_idx) => permutation_vec[idx] = *target_idx,
            None => permutation_vec[idx] = idx
        }
    }

    permutation_vec
}

/// Conversion from Mapping to a vector of vertex orbits, where each orbit contains two vertices 
fn orbits_from_mapping(
    mapping: &Mapping,
    length: usize
) -> Vec<core::orbit_ops::Orbit> {
    let mut orbits: Vec<core::orbit_ops::Orbit> = vec![];
    let mut indexes_visited: Vec<usize> = vec![];

    for idx in 0..length {
        match mapping.get(&idx) {
            Some(target_idx) => {
                if !indexes_visited.contains(&idx) {
                    orbits.push(vec![idx, *target_idx]);
                    indexes_visited.push(idx);
                    indexes_visited.push(*target_idx);
                }
            },
            None => () 
        }
    }

    orbits
}

/// Conversion from Permutation to a vector of vertex orbits 
pub fn orbits_from_permutation(
    permutation_vec: &Permuation,
    length: usize
) -> Vec<core::orbit_ops::Orbit> {
    orbits_from_mapping(
        &mapping_from_permutation(permutation_vec), length
    )
}

/// Conversion from Permutation to Mapping
pub fn mapping_from_permutation(
    permutation_vec: &Permuation
) -> Mapping {
    let mut mapping = Mapping::new();
    for i in 0..permutation_vec.len() {
        if i != permutation_vec[i] {
            mapping.insert(i, permutation_vec[i]);
        }
    }

    mapping
}

/// Auxiliary function for create combination of permutations from multiple vectors 
fn permutation_groups_of_orbit_indexes(
    v_counts: &Vec<usize>
) -> Vec<Vec<usize>> {
    let mut results: Vec<Vec<usize>> = vec![];

    if v_counts.len() == 1 {
        results = (0..v_counts[0]).map(|idx| vec![idx]).collect();
    } else if v_counts.len() > 1{
        let mut v_counts_down = v_counts.clone();
        if let Some(count) = v_counts_down.pop() {
            let results_down: Vec<Vec<usize>> = permutation_groups_of_orbit_indexes(&v_counts_down);
            for i in 0..count {
                let mut results_down_tmp = results_down.clone();
                while let Some(mut result_tmp) = results_down_tmp.pop() {
                    result_tmp.push(i);
                    results.push(result_tmp);
                }
            }
        }
    }

    results
}


/// Generate all possible permutations from the given orbits
pub fn generate_all_permutations(
    orbits: &Vec<core::orbit_ops::Orbit>,
    length: usize
) -> Vec<Permuation> {
    let factorials: Vec<Vec<Vec<usize>>> = orbits.iter().map(|rp| combinatorial::factorial_vec(rp)).collect();
    let counts: Vec<usize> = factorials.iter().map(|r| r.len()).collect();
    let combos: Vec<Vec<usize>> = permutation_groups_of_orbit_indexes(&counts);

    let mut results: Vec<Permuation> = vec![];
    for combo in combos.iter() {
        let mut mapping = Mapping::new();
        for i in 0..combo.len() {
            update_mapping(
                &orbits[i],
                &factorials[i][combo[i]],
                &mut mapping
            );
        }

        results.push(permutation_from_mapping(&mapping, length));
    }

    results
}

#[cfg(test)]
mod test_cnap_permutation {
    use super::*;

    #[test]
    fn test_update_mapping() {
        let mut mapping = Mapping::new();
        let vertex_from_indexes: Vec<usize> = vec![2, 5, 8, 10];
        let vertex_to_indexes: Vec<usize> = vec![5, 2, 10, 8];
        update_mapping(&vertex_from_indexes, &vertex_to_indexes, &mut mapping);
        assert_eq!(mapping.len(), 4);
        assert_eq!(mapping.get(&2), Some(&5));
        assert_eq!(mapping.get(&5), Some(&2));
        assert_eq!(mapping.get(&8), Some(&10));
        assert_eq!(mapping.get(&10), Some(&8));
    }

    #[test]
    fn test_orbits_from_permutation() {
        let test_data: Vec<(Permuation, Vec<core::orbit_ops::Orbit>)> = vec![
            (vec![0, 3, 5, 1, 4, 2, 6, 7], vec![vec![1, 3], vec![2, 5]]),
            (vec![0, 1, 2, 3, 4, 5, 6, 7], vec![])
        ];

        for td in test_data.iter() {
            let (permutation_vec, orbits) = td;
            assert_eq!(orbits_from_permutation(permutation_vec, 8), *orbits); 
        }
    }

    #[test]
    fn test_permutation_mapping_conversion() {
        let test_data: Vec<(Mapping, Permuation)> = vec![
            ([(1, 3), (3, 1), (2, 5), (5, 2)].iter().cloned().collect(), vec![0, 3, 5, 1, 4, 2, 6, 7]),
            (Mapping::new(), vec![0, 1, 2, 3, 4, 5, 6, 7])
        ];

        for td in test_data.iter() {
            let (mapping, permutation_vec) = td;
            let permutation_vec_converted = permutation_from_mapping(&mapping, 8);
            assert_eq!(permutation_vec_converted, *permutation_vec); 
            let mapping_converted = mapping_from_permutation(&permutation_vec);
            assert_eq!(mapping_converted, *mapping);
        }
    }

    #[test]
    fn test_permutation_groups_of_orbit_indexes() {
        let v_counts: Vec<usize> = vec![4, 2, 3, 2];
        let results: Vec<Vec<usize>> = permutation_groups_of_orbit_indexes(&v_counts);
        assert_eq!(results.len(), 4 * 2 * 3 * 2);
    }

    #[test]
    fn test_generate_all_permutations() {
        let residual_paritition: Vec<Vec<usize>> = vec![vec![1, 3], vec![2, 4, 6, 9], vec![0, 5, 8]];
        let all_permutations: Vec<Permuation> = generate_all_permutations(&residual_paritition, 7);
        assert_eq!(all_permutations.len(), 288);
    }
}
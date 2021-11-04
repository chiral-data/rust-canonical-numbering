// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Orbit Operations

/// Types
/// An orbit is the set of all atoms that are transformed from one into another by the actions of all automorphisms of a molecular graph.
/// Reference: Ovidiu Ivanciuc, Handbook of Chemoinformatics Wiley–VCH 2003, Chapter 12
pub type Orbit = Vec<usize>;

pub fn orbits_sort(orbits: &mut Vec<Orbit>) {
    for i in 0..orbits.len() {
        orbits[i].sort_unstable();
    }
    orbits.sort_unstable();
}

pub fn orbit_equal(orbit_1: &Orbit, orbit_2: &Orbit) -> bool {
    if orbit_1.len() != orbit_2.len() { return false }

    let mut orbit_1_cloned = orbit_1.to_vec();
    let mut orbit_2_cloned = orbit_2.to_vec();
    orbit_1_cloned.sort_unstable();
    orbit_2_cloned.sort_unstable();
    orbit_1_cloned == orbit_2_cloned
}

pub fn orbits_equal(orbits_1: &Vec<Orbit>, orbits_2: &Vec<Orbit>) -> bool {
    if orbits_1.len() != orbits_2.len() { return false }

    let mut orbits_1_cloned = orbits_1.to_vec();
    let mut orbits_2_cloned = orbits_2.to_vec();
    orbits_sort(&mut orbits_1_cloned);
    orbits_sort(&mut orbits_2_cloned);
    orbits_1_cloned == orbits_2_cloned
}


pub fn orbit_contain(orbits: &Vec<Orbit>, orbit: &Orbit) -> bool {
    let mut orbits_cloned = orbits.to_vec();
    for i in 0..orbits.len() {
        orbits_cloned[i].sort_unstable();
    }

    let mut orbit_cloned = orbit.to_vec();
    orbit_cloned.sort_unstable();

    orbits_cloned.contains(&orbit_cloned)
}

pub fn orbit_overlap(orbit_a: &Orbit, orbit_b: &Orbit) -> bool {
    for idx in orbit_a.iter() {
        if orbit_b.contains(idx) {
            return true
        }
    }

    return false
}

/// Where the large orbit contains all the elements from the small orbit
pub fn orbit_cover(large_orbit: &Orbit, small_orbit: &Orbit) -> bool {
    for idx in small_orbit.iter() {
        if !large_orbit.contains(idx) {
            return false
        }
    }

    true
}

/// Merge the orbits with common element, until all the orbits are not overlapping
pub fn orbits_self_merge(orbits: &mut Vec<Orbit>) {
    let mut are_orbits_overlapping = true;

    while are_orbits_overlapping {
        are_orbits_overlapping = false;
        for i in 0..orbits.len() {
            // let mut found_intersection = false;
            for j in (i+1)..orbits.len() {
                if  orbit_overlap(&orbits[i], &orbits[j]) {
                    are_orbits_overlapping = true;
                    let mut orbit_j = orbits[j].clone();
                    orbits[i].append(&mut orbit_j);
                    orbits[i].sort_unstable();
                    orbits[i].dedup();
                    orbits.remove(j);
                    break;
                }
            }

            if are_orbits_overlapping {
                break;
            }
        }
    }
}

#[cfg(test)]
mod test_core_orbit_ops {
    use super::*;

    #[test]
    fn test_orbits_self_merge() {
        type InputType1 = Vec<Orbit>;
        type InputType2 = Vec<Orbit>;
        type ReturnType = Vec<Orbit>;
        let test_data: Vec<(InputType1, InputType2, ReturnType)> = vec![
            (
                vec![vec![5, 7], vec![4, 8], vec![3, 9], vec![2, 10], vec![1, 11]],
                vec![vec![7, 9], vec![10, 11], vec![0, 6]],
                vec![vec![0, 6], vec![1, 2, 10, 11], vec![3, 5, 7, 9], vec![4, 8]]
            )
        ];

        for td in test_data.iter() {
            let (mut orbits_i1, orbits_i2, orbits_result) = td.clone();
            orbits_i1.append(&mut orbits_i2.clone());
            orbits_self_merge(&mut orbits_i1);
            orbits_sort(&mut orbits_i1);
            assert_eq!(orbits_i1, *orbits_result); 
        }
    }

    #[test]
    fn test_orbit_cover() {
        let params: Vec<(Orbit, Orbit, bool)> = vec![
            (vec![2, 3, 4], vec![1, 2], false),
            (vec![2, 3, 4, 5, 6], vec![5, 6, 3], true),
        ];

        for param in params.iter() {
            assert_eq!(orbit_cover(&param.0, &param.1), param.2);
        }
    }

    #[test]
    fn test_orbit_overlap() {
        let params: Vec<(Orbit, Orbit, bool)> = vec![
            (vec![2, 3, 4], vec![1, 2], true),
            (vec![2, 3, 4], vec![5, 6, 3], true),
            (vec![2, 8, 4], vec![5, 6, 3], false),
        ];

        for param in params.iter() {
            assert_eq!(orbit_overlap(&param.0, &param.1), param.2);
        }
    }

    #[test]
    fn test_orbits_equal() {
        let params: Vec<(Vec<Orbit>, Vec<Orbit>, bool)> = vec![
            (vec![vec![2], vec![3, 4]], vec![vec![1, 2]], false),
            (vec![vec![2], vec![3, 4]], vec![vec![3, 4], vec![1, 2]], false),
            (vec![vec![2, 1], vec![3, 4]], vec![vec![4, 3], vec![1, 2]], true),
        ];

        for param in params.iter() {
            assert_eq!(orbits_equal(&param.0, &param.1), param.2);
        }
    }

    #[test]
    fn test_orbit_contain() {
        let params: Vec<(Vec<Orbit>, Orbit, bool)> = vec![
            (vec![vec![2], vec![3, 4]], vec![1, 2], false),
            (vec![vec![2], vec![3, 4]], vec![3, 4], true),
            (vec![vec![2, 1], vec![3, 4]], vec![4, 3], true),
        ];

        for param in params.iter() {
            assert_eq!(orbit_contain(&param.0, &param.1), param.2);
        }
    }
}
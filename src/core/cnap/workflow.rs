// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::core;
use super::combinatorial;
use super::permutation;
use super::isomorphism;

/// Check whether brutal force computation is feasible or not accroding to the parameter COMPUTATION_POWER
pub fn is_computable(residual_orbits: &Vec<core::orbit_ops::Orbit>) -> bool {
    let mut computation: usize = 1;
    for rp in residual_orbits.iter() {
        if rp.len() > 20 { // factorial() will overflow
            return false
        }

        computation *= combinatorial::factorial(rp.len());
        if computation > core::config::COMPUTATION_POWER {
            return false
        }
    }

    true
}

/// Check whether any orbit contains only two elements or not
fn is_tow_folded_symmetry(residual_orbits: &Vec<core::orbit_ops::Orbit>) -> bool {
    if residual_orbits.len() == 0 {
        return false
    }

    let mut lengths: Vec<usize> = residual_orbits.iter()
        .map(|orbit| orbit.len())
        .collect();
    lengths.sort_unstable();

    lengths[lengths.len() - 1] == 2
}


/// Brutal-force checking to find out the symmetric orbits
pub fn get_symmetric_orbits(
    orbits_residual: &Vec<core::orbit_ops::Orbit>,
    edges: &Vec<(usize, usize, usize)>,
    length: usize,
    orbits_symmetry: &mut Vec<core::orbit_ops::Orbit>,
) {
    let all_permutations: Vec<permutation::Permuation> = permutation::generate_all_permutations(orbits_residual, length);
    for p in all_permutations.iter() {
        if isomorphism::is_automorphic(p, edges) {
            orbits_symmetry.append(&mut permutation::orbits_from_permutation(p, length));
            core::orbit_ops::orbits_self_merge(orbits_symmetry);
        }
    }
}

pub enum ErrorCNAP {
    ErrorIncomputable,
    ErrorTwoFolded,
    // ErrorHighSymmetry,
}

// The CNAP process
pub fn run<T: core::graph::VertexExtendableHash>(
    edges: &Vec<(usize, usize, usize)>,
    length: usize,
    orbits_residual: &Vec<core::orbit_ops::Orbit>,
    orbits_symmetry: &mut Vec<core::orbit_ops::Orbit>,
) -> Result<(), ErrorCNAP> {
    if is_computable(orbits_residual) {
        get_symmetric_orbits(orbits_residual, edges, length, orbits_symmetry);
        Ok(())
    } else {
        if is_tow_folded_symmetry(orbits_residual) { // Case two-folded symmetry 
            // two-folded symmetry cannot be handled by graph reduction
            // it is well worth trying automorphic checking on switching every two vertices inside each orbit
            let mut p: permutation::Permuation = (0..length).collect();
            for orbit in orbits_residual.iter() {
                p[orbit[0]] = orbit[1];
                p[orbit[1]] = orbit[0];
            }
            
            if isomorphism::is_automorphic(&p, edges) {
                orbits_symmetry.append(&mut orbits_residual.clone());
                core::orbit_ops::orbits_self_merge(orbits_symmetry);
                return Ok(());
            } else {
                return Err(ErrorCNAP::ErrorTwoFolded)
            }
        }

        Err(ErrorCNAP::ErrorIncomputable)
    }
}

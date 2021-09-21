use crate::core;
use super::combinatorial;
use super::permutation;
use super::isomorphism;


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

//
// Predicate: two folded symmetry
// 
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


//
// Brutal-force Calculation of Symmetric Orbits
//
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

//
// The Process
//
pub fn run<T: core::graph::VertexExtendableHash>(
    edges: &Vec<(usize, usize, usize)>,
    length: usize,
    orbits_residual: &Vec<core::orbit_ops::Orbit>,
    orbits_symmetry: &mut Vec<core::orbit_ops::Orbit>,
) -> Result<(), ErrorCNAP> {
    if is_computable(orbits_residual) { // Brutal-force CNAP
        if cfg!(debug_assertions) {
            println!("CNAP Computable {:?}", orbits_residual);
        }

        get_symmetric_orbits(orbits_residual, edges, length, orbits_symmetry);
        Ok(())
    } else {
        // Case Two-Folded Symmetry
        if is_tow_folded_symmetry(orbits_residual) {
            if cfg!(debug_assertions) {
                println!("Proceed to Two-folded Symmetry: {:?}", orbits_residual);
            }

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

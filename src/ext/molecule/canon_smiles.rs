// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//!
//! Generation of canonicalized SMILES for molecules
//!     Use GIVP orbits directly, assume that GIVP orbits were symmetric orbits
//! 

use crate::core;
use super::atom;
use super::bond;
use super::extendable_hash;
use super::molecule;

/// Implement the Trait required by smiles_writer from crate chem
impl chem::smiles_writer::TraitMoleculeForSMILES for molecule::Molecule {
    fn get_neighbours_of_atom(&self, atom: &usize) -> Vec<usize> {
        self.atoms[*atom].bonds.iter()
            .map(|b| b.tid)
            .collect()
    }

    fn get_bond_symbol(&self, atom_1: &usize, atom_2: &usize) -> String {
        let b: Vec<bond::Bond> = self.atoms[*atom_1].bonds.clone().into_iter()
            .filter(|b| b.tid == *atom_2)
            .collect();

        if b.len() > 0 {
            b[0].bond_char()
        } else {
            String::from("No such bond")
        }
    }

    fn get_atom_symbol(&self, atom: &usize) -> String {
        self.atoms[*atom].kind.to_string()
    }

    fn get_atom_ranking(&self, atom: &usize, rankings: &Vec<usize>) -> usize {
        rankings[*atom]
    }

    fn count_of_atoms(&self) -> usize {
        self.atoms.len()
    }
}

/// Break symmetry to get canonical numbering
fn get_canon_numbering(
    vv: &core::graph::VertexVec<atom::Atom>,
    orbits_givp: &Vec<core::orbit_ops::Orbit>,
    numbering_givp: &Vec<usize>
) -> Vec<usize> {
    if orbits_givp.len() == 0 {
        numbering_givp.clone()
    } else {
        let mut orbits = orbits_givp.clone();
        let mut numbering = numbering_givp.clone();
        while orbits.len() > 0 {
            orbits.sort_by_key(|ob| numbering[ob[0]]);
            orbits.reverse();
            for i in 1..orbits[0].len() {
                numbering[orbits[0][i]] = numbering[orbits[0][0]] - 1;
            }
            core::givp::run::<extendable_hash::AtomExtendable>(&vv, &mut numbering, &mut orbits);
        }

        numbering
    }
}

/// Generate canonical SMILES from the input SMILES string
///     stereochemistry is not taken into consideration at the moment
pub fn get_canon_smiles(smiles: &String) -> String {
    let mol = molecule::Molecule::from_smiles(smiles);
    let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());
    let mut orbits: Vec<core::orbit_ops::Orbit> = vec![];
    let mut numbering: Vec<usize> = vec![];
    core::givp::run::<extendable_hash::AtomExtendable>(&vv, &mut numbering, &mut orbits);
    numbering = get_canon_numbering(&vv, &orbits, &numbering);
    chem::smiles_writer::write_smiles_for_mol(&mol, &numbering)
}

#[cfg(test)]
mod test_ext_mol_canon_smiles {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_get_canon_smiles() {
        type InputType1 = String;
        type ReturnType = String; 
        let test_data: Vec<(InputType1, ReturnType)> = vec![
            (
                "c1ccccc1CN", 
                "NCc1ccccc1"
            ), 
            (
                "COc1ccc(C(=O)C[n+]2c(C)n(Cc3c4c(cc5c3OCC5)OCC4)c3ccccc32)cc1",
                "COc1ccc(cc1)C(=O)C[n+]1c2ccccc2n(Cc2c3CCOc3cc3CCOc23)c1C",
            ),
            (
                r#"Cc1[nH]c2ccccc2c1/C=C1\SC(=N)N(c2nccs2)C1=O"#, // 926178 
                "Cc1[nH]c2ccccc2c1C=C1SC(=N)N(C1=O)c1nccs1"
            ),
            (
                "C(C)(C)CCNCCC(C)(C)", 
                "CC(C)CCNCCC(C)C"
            ) 
        ].into_iter().map(|td| (td.0.to_string(), td.1.to_string())).collect();

        let mut rng = rand::thread_rng();

        for td in test_data.iter() {
            let (smiles, results) = td;
            assert_eq!(get_canon_smiles(smiles), *results);

            // random test
            for _ in 0..5 {
                let mol = molecule::Molecule::from_smiles(smiles);
                let times = rng.gen_range(0..5);
                let mut random_rankings: Vec<usize> = (0..mol.atoms.len()).collect();
                for _ in 0..times{
                    let atom_1 = rng.gen_range(0..mol.atoms.len());
                    let atom_2 = rng.gen_range(0..mol.atoms.len());
                    let temp = random_rankings[atom_1];
                    random_rankings[atom_1] = random_rankings[atom_2];
                    random_rankings[atom_2] = temp;
                }

                let random_smiles = chem::smiles_writer::write_smiles_for_mol(&mol, &random_rankings); 
                assert_eq!(get_canon_smiles(&random_smiles), *results);
            }
        }
    }
}
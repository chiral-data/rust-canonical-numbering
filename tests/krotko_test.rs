// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Test on examples from Dr. Krotko's paper: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00453-4

use graph_canonicalization::core;
use graph_canonicalization::ext::molecule;

#[test]
fn krotko_test() {
    type InputType1 = String;
    type InputType2 = String;
    let test_data: Vec<(InputType1, InputType2)> = vec![
        ("C12C3C1C3C1C3C2C13", "figure 2",),
        // ("C1(C2C3C4C15)C6C7C2C8C3C9C%10C4C%11C5C6C%12C%11C%10C%13C%12C7C8C9%13", "figure 3"),
        ("NC1=CC=C1", "figure 6",),
        ("C12C3C4C5C1C6C7C2C8C3C6C5C8C74", "Petersen Graph",),
        ("C1OC23COC45COC11COC67COC8(COC9(CO2)COC(CO1)(CO6)OCC(CO9)(OC4)OCC(CO5)(OC7)OC8)OC3", "Shelley and Munk",)
    ].into_iter().map(|s| (s.0.to_string(), s.1.to_string())).collect();

    for td in test_data.iter() {
        let (smiles, source) = td.clone();
        let mol = molecule::molecule::Molecule::from_smiles(&smiles);
        if cfg!(debug_assertions) {
            println!("\n\nTest for {} from {}", smiles, source);
            println!("============================================");
            println!("SMILES with index:\n{}", mol.smiles_with_index(&smiles, &vec![]));
        }

        let mut orbits_partitioned: Vec<core::orbit_ops::Orbit> = vec![];
        let mut orbits_symmetry: Vec<core::orbit_ops::Orbit> = vec![];
        let mut numbering: Vec<usize> = vec![];
        molecule::workflow::canonical_numbering_and_symmetry_perception(&mol.atoms, &mut orbits_partitioned, &mut orbits_symmetry, &mut numbering);
        println!("SMILES with numbering:\n{}", mol.smiles_with_index(&smiles, &numbering));
        if cfg!(debug_assertions) {
            core::orbit_ops::orbits_sort(&mut orbits_partitioned);
            core::orbit_ops::orbits_sort(&mut orbits_symmetry);
            println!("GIVP: {:?}\nCNAP: {:?}", orbits_partitioned, orbits_symmetry);
        }
        assert_eq!(core::orbit_ops::orbits_equal(&orbits_partitioned, &orbits_symmetry), true);
    }
}
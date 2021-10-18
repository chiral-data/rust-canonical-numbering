/// Test on Database ChemBL
///

#[cfg(test)]
mod test_chembl {
    use crate::ext::molecule;
    use crate::core;
    use std::io::prelude::*;

    static CHEMBL_FILE: &str = "/Users/qw/Documents/data/ChEMBL/chembl_28_chemreps.txt";

    fn read_lines<P>(filename: P) -> std::io::Result<std::io::Lines<std::io::BufReader<std::fs::File>>>
    where P: AsRef<std::path::Path>, {
        let file = std::fs::File::open(filename)?;
        Ok(std::io::BufReader::new(file).lines())
    }

    fn clean_smiles(smiles_origin: &String) -> String {
        let mut parts: Vec<&str> = smiles_origin.split('.').collect();
        match parts.len() {
            1 => String::from(parts[0]),
            0 => panic!("Invalid smiles!!!"),
            _ => {
                parts.sort_by(|a, b| a.len().cmp(&b.len()));
                String::from(parts[parts.len()-1])
            }
        }
    }

    #[test]
    fn chembl() {
        if let Ok(lines) = read_lines(CHEMBL_FILE) {
            let mut count: usize = 0;
            let mut count_large_molecule: usize = 0;
            let skip: usize = 0;

            for line in lines {
                count += 1;
                if count == 1 { // skip headline
                    continue;
                }

                if count <= skip {
                    continue;
                }

                if let Ok(chembl_line) = line {
                    let parts: Vec<&str> = chembl_line.split('\t').collect();
                    if parts.len() == 4 {
                        // let now = std::time::Instant::now();
                        let smiles: String = clean_smiles(&String::from(parts[1]));
                        let mol = molecule::molecule::Molecule::from_smiles(&smiles);
                        if mol.atoms.len() == 0 {
                            continue;
                        }

                        if mol.atoms.len() > 250 {
                            println!("large molecule");
                            count_large_molecule += 1;
                            continue; // ignore large molecule
                        }

                        println!("\n{}\n{}\n{}", count, smiles, mol.smiles_with_index(&smiles, &vec![]));
                        let mut orbits_partitioned: Vec<core::orbit_ops::Orbit> = vec![];
                        let mut orbits_symmetry: Vec<core::orbit_ops::Orbit> = vec![];
                        let mut numbering: Vec<usize> = vec![];

                        molecule::workflow::canonical_numbering_and_symmetry_perception(&mol.atoms, &mut orbits_partitioned, &mut orbits_symmetry, &mut numbering);

                        if !core::orbit_ops::orbits_equal(&orbits_partitioned, &orbits_symmetry) {
                            core::orbit_ops::orbits_sort(&mut orbits_partitioned);
                            core::orbit_ops::orbits_sort(&mut orbits_symmetry);
                            println!("GIAP failed {}:\n{}\nGIAP orbits {:?}\nCNAP orbits {:?}\n", count, mol.smiles_with_index(&smiles, &vec![]), orbits_partitioned, orbits_symmetry);
                        }
                    } else {
                        println!("Parsing Error on {}", chembl_line)
                    }
                }
            }

            println!("large molecule count: {}", count_large_molecule);
        }
    }
}
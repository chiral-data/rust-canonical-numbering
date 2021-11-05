// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use super::atom;

pub struct Molecule {
    pub atoms: Vec<atom::Atom>,
}

impl Molecule {
    pub fn from_smiles(smiles: &str) -> Self {
        let mut builder = purr::graph::Builder::new();

        match purr::read::read(smiles, &mut builder, None) {
            Ok(_) => {
                let mut atoms = builder.build().expect("atoms");
                // lib 'purr' indicates an aromatic bond as bondkind::elided instead of bondkind::aromatic 
                //      if bond symbol ':' does not occur explicitly in molecule smiles.
                // for the applications which requires a clear differentiation between single bond and aromatic bond,
                //      bondkind::elided is not enough.
                let aromatic_flags: Vec<bool> = atoms.iter().map(|atom| atom.is_aromatic()).collect();
                for atom_idx in 0..(atoms.len()) {
                    for bond in atoms[atom_idx].bonds.iter_mut() {
                        if aromatic_flags[atom_idx] && aromatic_flags[bond.tid] {
                            *bond = purr::graph::Bond::new(purr::feature::BondKind::Aromatic, bond.tid);
                        }
                    }
                }

                let new_atoms: Vec<atom::Atom> = atoms.iter().map(|a| atom::Atom::from_atom_purr(&a)).collect();
                Self { atoms: new_atoms }
            },
            Err(e) => {
                println!("error smiles parsing: {:?}", e);
                Self { atoms: vec![] }
            }
        }
    }
    
    pub fn smiles_with_index(&self, smiles: &String, numbering: &Vec<usize>) -> String {
        if self.atoms.len() == 0 {
            return "smiles parsing errro".to_string()
        }

        let mut new_smiles: String = String::from("");
        let mut cur: usize = 0;

        for i in 0..(self.atoms.len()-1) {
            let atom_string: String = self.atoms[i].kind.to_string();
            let mapped_number = match numbering.len() > 0 {
                true => numbering[i],
                false => i
            };
            if atom_string.as_bytes()[0] == "[".as_bytes()[0] {
                new_smiles += &format!("[{atom_string}:{index}]", atom_string=String::from(&atom_string[1..(atom_string.len()-1)]), index=mapped_number);
                while smiles.as_bytes()[cur] != "]".as_bytes()[0] {
                    cur += 1;
                }
                cur += 1;
            } else {
                new_smiles += &format!("[{atom_string}:{index}]", atom_string=atom_string, index=mapped_number);
                cur += atom_string.len();
            }

            let cur_last: usize = cur;
            let next_atom_string: String = self.atoms[i+1].kind.to_string();
            while next_atom_string.as_bytes()[0] != smiles.as_bytes()[cur] {
                cur += 1;
                if cur >= smiles.len() {
                    break;
                }
            }

            new_smiles += &smiles[cur_last..cur];
        }


        let mapped_number = match numbering.len() > 0 {
            true => numbering[self.atoms.len() - 1],
            false => self.atoms.len() - 1
        };
        let atom_string: String = self.atoms[self.atoms.len() - 1].kind.to_string();
        new_smiles += &format!("[{atom_string}:{index}]", atom_string=atom_string, index=mapped_number);
        cur += atom_string.len();
        new_smiles += &smiles[cur..smiles.len()];

        new_smiles
    }
}

#[cfg(test)]
mod test_ext_mol_molecule {
    use super::*;

    #[test]
    fn test_from_smiles() {
        let smiles: String = "c1ccccc1CN".to_string();
        let mol = Molecule::from_smiles(&smiles);
        assert_eq!(mol.atoms[0].kind, "c".to_string());
        assert_eq!(mol.atoms[0].bonds.len(), 2);
        assert_eq!(mol.atoms[1].kind, "c".to_string());
        assert_eq!(mol.atoms[1].bonds.len(), 2);
        assert_eq!(mol.atoms[2].kind, "c".to_string());
        assert_eq!(mol.atoms[2].bonds.len(), 2);
        assert_eq!(mol.atoms[3].kind, "c".to_string());
        assert_eq!(mol.atoms[3].bonds.len(), 2);
        assert_eq!(mol.atoms[4].kind, "c".to_string());
        assert_eq!(mol.atoms[4].bonds.len(), 2);
        assert_eq!(mol.atoms[5].kind, "c".to_string());
        assert_eq!(mol.atoms[5].bonds.len(), 3);
        assert_eq!(mol.atoms[6].kind, "C".to_string());
        assert_eq!(mol.atoms[6].bonds.len(), 2);
        assert_eq!(mol.atoms[7].kind, "N".to_string());
        assert_eq!(mol.atoms[7].bonds.len(), 1);
    }
}
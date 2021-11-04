// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::core;
use super::element;
use super::bond;

/// Hash atom charge
/// Postive charge value > 200, Neutral = 100, Negative < 100
const DEFAULT_CHARGE_HASH_VALUE: usize = 100;
fn get_charge_hash_value_from_atom_purr(atom_purr: &purr::graph::Atom) -> usize {
    match &atom_purr.kind {
        purr::feature::AtomKind::Star => DEFAULT_CHARGE_HASH_VALUE,
        purr::feature::AtomKind::Aliphatic(_) => DEFAULT_CHARGE_HASH_VALUE, 
        purr::feature::AtomKind::Aromatic(_) => DEFAULT_CHARGE_HASH_VALUE, 
        purr::feature::AtomKind::Bracket {
            isotope: _, symbol: _, hcount: _, configuration: _, charge, map: _
        } => match charge {
                Some(c) => {
                    let charge_of_atom: i8 = c.into();
                    match charge_of_atom > 0 {
                        true => 200 + charge_of_atom as usize,
                        false => (0 - charge_of_atom) as usize
                    }
                },
                None => DEFAULT_CHARGE_HASH_VALUE
            }
        } 
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub kind: String, // purr::feature::AtomKind of the current version does not support trait Clone, so use String instead
    pub bonds: Vec<bond::Bond>,
    charge: usize,
    atomic_number: usize,
    is_aromatic: bool,
    is_breaking_symmetry: bool,
}


impl Atom {
    pub fn from_atom_purr(atom_purr: &purr::graph::Atom) -> Self {
        let atomic_number: usize = element::atomic_number(&atom_purr.kind) as usize;
        let is_aromatic: bool = atom_purr.is_aromatic();
        let kind: String = atom_purr.kind.to_string();
        let bonds: Vec<bond::Bond> = atom_purr.bonds
            .iter()
            .map(|b| bond::Bond::from_bond_purr(&b))
            .collect();
        let is_breaking_symmetry: bool = false;
        let charge: usize = get_charge_hash_value_from_atom_purr(atom_purr);

        Self { kind, bonds, charge, atomic_number, is_aromatic, is_breaking_symmetry }
    }
}

impl core::graph::Vertex for Atom {
    fn fixed_hash_value(&self) -> core::graph::VertexFixedHashValue {
        self.charge * 10 * 1000 * 10 * 10 
        + self.bonds.len() * 10 * 1000 * 10  
        + self.is_aromatic as usize * 10 * 1000
        + self.atomic_number as usize * 10 
        + self.is_breaking_symmetry as usize
    }

    fn break_symmetry_vertex(&mut self) {
        self.is_breaking_symmetry = true;
    }

    fn break_symmetry_edge(&mut self, tid: usize) {
        for bond in self.bonds.iter_mut() {
            if bond.tid == tid {
                bond.break_symmetry();
            }
        }
    }

    fn degree(&self) -> usize {
        self.bonds.len()
    }
    
    fn neighbour_indexes(&self) -> Vec<usize> {
        self.bonds.iter()
            .map(|b| b.tid)
            .collect()
    }

    fn custom_new_in_reduced_graph(
        self_index: usize,
        atomic_number: usize, 
        edges_from: &Vec<(usize, usize)>,
        edges_to: &Vec<(usize, usize)>,
        atoms: &Vec<Atom>,
        numbering: &Vec<usize>
    ) -> Self {
        let kind: String = "CustomAtom".to_string();
        let mut bonds: Vec<bond::Bond> = vec![];
        let is_aromatic: bool = false;
        let is_breaking_symmetry: bool = false;
        let charge: usize = DEFAULT_CHARGE_HASH_VALUE;

        // bonds construction
        for idx in 0..edges_from.len() {
            let edge_from = edges_from[idx];
            let bonds_found: Vec<bond::Bond> = atoms[edge_from.0].bonds.clone().into_iter()
                .filter(|b| b.tid == edge_from.1)
                .collect();

            if bonds_found.len() != 1 {
                panic!("Create New Atom: cannot find bond for edges: {:?}\n bonds found: {:?}\n from atom bonds: {:?}", edges_from, bonds_found, atoms[edge_from.0].bonds);
            }

            let mut new_bond: bond::Bond = bonds_found[0].clone();
            new_bond.set_kind_value(numbering, &(self_index, new_bond.tid));
            new_bond.tid = edges_to[idx].1;
            bonds.push(new_bond);
        }

        Self { kind, bonds, charge, atomic_number, is_aromatic, is_breaking_symmetry }
    }

    fn update_edges_in_reduced_graph(&mut self,
        self_index: usize,
        reduced_groups: &std::collections::HashMap<usize, usize>,
        numbering: &Vec<usize>
    ) {
        for bond in self.bonds.iter_mut() {
            match reduced_groups.get(&bond.tid) {
                Some(&target_index) => {
                    bond.set_kind_value(numbering, &(self_index, bond.tid));
                    bond.tid = target_index;
                },
                None => ()
            }
        }
    }

    fn custom_new_in_separated_graph(
        atomic_number: usize,
        vertex: usize,
        valid_neighbours: &Vec<usize>,
        atoms: &Vec<Atom>
    ) -> Self {
        let kind: String = "CustomAtom".to_string();
        let is_aromatic: bool = false;
        let is_breaking_symmetry: bool = false;
        let bonds: Vec<bond::Bond> = atoms[vertex].bonds.clone().into_iter()
            .filter(|b| valid_neighbours.contains(&b.tid))
            .collect();
        let charge: usize = DEFAULT_CHARGE_HASH_VALUE;
            
        Self { kind, bonds, charge, atomic_number, is_aromatic, is_breaking_symmetry }
    }

    fn custom_new_in_folded_graph(
        atomic_number: usize,
        boundary_edge: &(usize, usize),
        atoms: &Vec<Atom>
    ) -> Self {
        let kind: String = "CustomAtom".to_string();
        let bonds: Vec<bond::Bond> = atoms[boundary_edge.1].bonds.clone().into_iter()
            .filter(|b| b.tid == boundary_edge.0)
            .collect();
        let is_aromatic: bool = false;
        let is_breaking_symmetry: bool = false;
        let charge: usize = DEFAULT_CHARGE_HASH_VALUE;

        Self { kind, bonds, charge, atomic_number, is_aromatic, is_breaking_symmetry }
    }

    fn debug_print(&self) {
        println!("Atomic number: {}, Bonds: {:?}", self.atomic_number, self.bonds);
    }
}

#[cfg(test)]
mod test_ext_mol_atom {
    use super::*;
    use crate::core::graph::*;
    use crate::ext::molecule::molecule;

    #[test]
    fn test_atom_fixed_hash_value() {
        let mol = molecule::Molecule::from_smiles("c1ccccc1CN");
        assert_eq!(mol.atoms[0].fixed_hash_value(), 100210060);
        assert_eq!(mol.atoms[5].fixed_hash_value(), 100310060);
        assert_eq!(mol.atoms[6].fixed_hash_value(), 100200060);
    }

    #[test]
    fn test_break_symmetry() {
        let mut mol = molecule::Molecule::from_smiles("c1ccccc1CN");
        mol.atoms[6].break_symmetry_vertex();
        assert_eq!(mol.atoms[6].fixed_hash_value(), 100200061);
        mol.atoms[6].break_symmetry_edge(7);
        assert_eq!(mol.atoms[6].bonds[0].fixed_hash_value(), 10);
        assert_eq!(mol.atoms[6].bonds[1].fixed_hash_value(), 11);
 
    }

    #[test]
    fn test_custom_new_in_reduce_graph() {
        let mol = molecule::Molecule::from_smiles("c1ccccc1CN");
        let new_custom = Atom::custom_new_in_reduced_graph(0, 200, &vec![(6, 7)], &vec![(1, 2)], &mol.atoms, &vec![0;8]);
        assert_eq!(new_custom.kind, "CustomAtom".to_string());
        assert_eq!(new_custom.bonds.len(), 1);
        assert_eq!(new_custom.bonds[0].tid, 2);
    }

    #[test]
    fn test_custom_new_in_folded_graph() {
        let mol = molecule::Molecule::from_smiles("c1ccccc1CN");
        let new_custom = Atom::custom_new_in_folded_graph(200, &(5, 6), &mol.atoms);
        assert_eq!(new_custom.kind, "CustomAtom".to_string());
        assert_eq!(new_custom.bonds.len(), 1);
        assert_eq!(new_custom.bonds[0].tid, 5);
    }

    #[test]
    fn test_update_edges_in_reduced_graph() {
        // let mol = molecule::Molecule::from_smiles("c1cc2cccc3c4cccc5cccc(c(c1)c23)c54"); // 2065990
        // let new_custom = Atom::custom_new_in_folded_graph(200, &(5, 6), &mol.atoms);
        // unimplemented!();
        // TBI
    }

    #[test]
    fn test_charge() {
        type InputType1 = String;
        type InputType2 = usize;
        type ReturnType = usize;
        let test_data: Vec<(InputType1, InputType2, ReturnType)> = vec![
            ("c1ccc2cc3ccccc3cc2c1", 1, DEFAULT_CHARGE_HASH_VALUE),
            ("c1[c+]cc2cc3ccccc3cc2c1", 1, 201),
            ("c1[n-]cc2cc3ccccc3cc2c1", 1, 1),
        ].iter().map(|td| (td.0.to_string(), td.1.clone(), td.2.clone())).collect();

        for td in test_data.iter() {
            let (smiles, vertex_index, charge_hash_value) = td;
            let mol = molecule::Molecule::from_smiles(&smiles);
            assert_eq!(mol.atoms[*vertex_index].charge, *charge_hash_value);
        }
        
    }
}
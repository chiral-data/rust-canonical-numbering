// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Extendable Hash Value for Atom

use crate::core;
use crate::core::graph::*;
use super::config;
use super::bond;
use super::atom;


#[derive(Clone)]
pub struct AtomExtendable {
    indexes_atom_visited: Vec<usize>,
    indexes_atom_boundary: Vec<usize>,
    bonds: Vec<bond::Bond>,
}

impl core::graph::VertexExtendableHash for AtomExtendable {
    type VertexType = atom::Atom;

    fn init(atom_idx: usize) -> Self {
        let indexes_atom_visited: Vec<usize> = vec![atom_idx];
        let indexes_atom_boundary: Vec<usize> = vec![atom_idx];
        let bonds: Vec<bond::Bond> = vec![]; 

        Self { indexes_atom_visited, indexes_atom_boundary, bonds }
    }

    fn extend(&mut self, atoms: &Vec<atom::Atom>) {
        let mut indexes_atom_boundary: Vec<usize> = vec![];
        let mut indexes_atom_visited: Vec<usize> = self.indexes_atom_visited.clone();
        self.bonds = vec![];

        for &idx in self.indexes_atom_boundary.iter() {
            for bond in atoms[idx].bonds.iter() {
                if !self.indexes_atom_visited.contains(&bond.tid) {
                    indexes_atom_visited.push(bond.tid);
                    indexes_atom_boundary.push(bond.tid);
                    self.bonds.push(bond.clone());
                }
            }
        }

        self.indexes_atom_visited = indexes_atom_visited;
        self.indexes_atom_boundary = indexes_atom_boundary;
    }

    fn is_completed(&self) -> bool {
        self.indexes_atom_boundary.len() == 0
    }

    fn value(&self, atom_invariant_values: &Vec<core::graph::VertexFixedHashValue>) -> core::graph::VertexExtendableHashValue { 
        let mut values: Vec<usize> = self.bonds.iter()
            .map(|b| b.fixed_hash_value() * config::MAX_OF_ATOM_FIXED_HASH_VALUE + atom_invariant_values[b.tid])
            .collect();
        values.sort_unstable();
        values.reverse();
        values
    }
}

#[cfg(test)]
mod test_ext_mol_atom_ext {
    use super::*;
    use super::super::molecule;
    use crate::core;

    #[test]
    fn test_atom_repr() {
        let mol = molecule::Molecule::from_smiles("c1ccccc1CN");
        let aivs: Vec<core::graph::VertexFixedHashValue> = mol.atoms.iter()
            .map(|a| a.fixed_hash_value())
            .collect();
        let mut ar = AtomExtendable::init(5);
        assert_eq!(ar.value(&aivs), vec![]);
        assert_eq!(ar.is_completed(), false);
        ar.extend(&mol.atoms);
        assert_eq!(ar.value(&aivs), vec![5100210060, 5100210060, 1100200060]);
        assert_eq!(ar.is_completed(), false);
        ar.extend(&mol.atoms);
        assert_eq!(ar.value(&aivs), vec![5100210060, 5100210060, 1100100070]);
        assert_eq!(ar.is_completed(), false);
        ar.extend(&mol.atoms);
        assert_eq!(ar.value(&aivs), vec![5100210060, 5100210060]);
        assert_eq!(ar.is_completed(), false);
    }
}
// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::core;
use super::config;

#[derive(Debug, Clone)]
pub struct Bond {
    kind: purr::feature::BondKind,
    kind_value: usize,
    pub tid: usize,
    is_breaking_symmetry: bool
}

impl Bond {
    pub fn from_bond_purr(bond_purr: &purr::graph::Bond) -> Self {
        let mut new_bond = Self { kind: bond_purr.kind.clone(), kind_value: 0, tid: bond_purr.tid, is_breaking_symmetry: false };
        new_bond.kind_value = new_bond.bond_type();
        new_bond
    }

    fn bond_type(&self) -> usize {
        match self.kind {
            purr::feature::BondKind::Elided => 1,
            purr::feature::BondKind::Single => 1,
            purr::feature::BondKind::Double => 2,
            purr::feature::BondKind::Triple => 3,
            purr::feature::BondKind::Quadruple => 4,
            purr::feature::BondKind::Aromatic => 5,
            purr::feature::BondKind::Up => 1,
            purr::feature::BondKind::Down => 1
        }
    }

    pub fn bond_char(&self) -> String {
        match self.kind {
            purr::feature::BondKind::Elided => String::from(""),
            purr::feature::BondKind::Single => String::from(""),
            purr::feature::BondKind::Double => String::from("="),
            purr::feature::BondKind::Triple => String::from("#"),
            purr::feature::BondKind::Quadruple => String::from("$"),
            purr::feature::BondKind::Aromatic => String::from(""),
            purr::feature::BondKind::Up => String::from(""),
            purr::feature::BondKind::Down => String::from("")
        }
    }

    pub fn break_symmetry(& mut self) {
        self.is_breaking_symmetry = true;
    }

    pub fn set_kind_value(&mut self, numbering: &Vec<usize>, vertex_pair: &(usize, usize)) { // for reduced graph, connection to different part of the custom group should be differenciated. Otherwise, new symmetry will be created.
        let vertex_weight: usize = match numbering[vertex_pair.0] > numbering[vertex_pair.1] {
            true => numbering[vertex_pair.1] * config::MAX_COUNT_OF_ATOMS_IN_A_MOLECULE + numbering[vertex_pair.0],
            false => numbering[vertex_pair.0] * config::MAX_COUNT_OF_ATOMS_IN_A_MOLECULE + numbering[vertex_pair.1]
        };
        self.kind_value = vertex_weight * 10 + self.bond_type();
    }
}

impl core::graph::Edge for Bond {
    fn fixed_hash_value(&self) -> core::graph::EdgeFixedHashValue {
        self.kind_value * 10 + self.is_breaking_symmetry as usize
    }
}


#[cfg(test)]
mod test_ext_mol_bond {
    use super::*;
    use crate::core::graph::*;

    #[test]
    fn test_bond() {
        let b: Bond = Bond { kind: purr::feature::BondKind::Aromatic, kind_value: 5, tid: 3, is_breaking_symmetry: false };
        assert_eq!(b.fixed_hash_value(), 50);
    }

    #[test]
    fn test_break_symmetry() {
        let mut b: Bond = Bond { kind: purr::feature::BondKind::Aromatic, kind_value: 5, tid: 3, is_breaking_symmetry: false };
        b.break_symmetry();
        assert_eq!(b.fixed_hash_value(), 51);
    }

    #[test]
    fn test_kind_value() {
        let mut b: Bond = Bond { kind: purr::feature::BondKind::Aromatic, kind_value: 1, tid: 3, is_breaking_symmetry: false };
        b.set_kind_value(&vec![2, 3], &(0, 1));
        assert_eq!(b.fixed_hash_value(), 200350);
        b.set_kind_value(&vec![3, 2], &(0, 1));
        assert_eq!(b.fixed_hash_value(), 200350);
    }
}
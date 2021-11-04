// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Graph Trait for Vertex & Edge & Vector of Vertices

///  Types
pub type VertexFixedHashValue = usize;
pub type VertexExtendableHashValue = Vec<usize>;
pub type EdgeFixedHashValue = usize;

///  Trait Vertex
pub trait Vertex: Clone {
    fn fixed_hash_value(&self) -> VertexFixedHashValue;
    fn break_symmetry_vertex(&mut self);
    fn break_symmetry_edge(&mut self, tid: usize);
    fn degree(&self) -> usize;
    fn neighbour_indexes(&self) -> Vec<usize>;
    fn custom_new_in_reduced_graph(
        self_index: usize,
        custom_marker: usize, 
        edges_from: &Vec<(usize, usize)>,
        edges_to: &Vec<(usize, usize)>,
        vertices: &Vec<Self>,
        numbering: &Vec<usize>
    ) -> Self;
    fn update_edges_in_reduced_graph(
        &mut self,
        self_index: usize,
        reduced_groups: &std::collections::HashMap<usize, usize>,
        numbering: &Vec<usize>
    );
    fn custom_new_in_folded_graph(
        custom_marker: usize,
        boundary_edge: &(usize, usize),
        vertices: &Vec<Self>
    ) -> Self;
    fn custom_new_in_separated_graph(
        custom_maker: usize,
        vertex: usize,
        valid_neighbours: &Vec<usize>,
        vertices: &Vec<Self>
    ) -> Self;
    fn debug_print(&self);
}

/// Trait Edge
pub trait Edge: Clone {
    fn fixed_hash_value(&self) -> EdgeFixedHashValue;
}

/// Trait Extendable Hash for Vertex
pub trait VertexExtendableHash: Clone {
    type VertexType: Vertex;

    fn init(vertice_idx: usize) -> Self;
    fn extend(&mut self, vertices: &Vec<Self::VertexType>);
    fn is_completed(&self) -> bool;
    fn value(&self, fixed_hash_values: &Vec<VertexFixedHashValue>) -> VertexExtendableHashValue; 
}

/// Trait Vertices
#[derive(Debug, Clone)]
pub struct VertexVec<T: Vertex> {
    indexes: Vec<usize>,
    vertices: Vec<T>
}

impl<T: Vertex> VertexVec<T> {
    pub fn init(indexes: Vec<usize>, vertices: Vec<T>) -> Self {
        Self { indexes, vertices }
    }
    
    pub fn len(&self) -> usize {
        self.indexes.len()
    }

    pub fn all_len(&self) -> usize {
        self.vertices.len()
    }

    pub fn valid_indexes(&self) -> &Vec<usize> {
        &self.indexes
    }

    pub fn all_vertices(&self) -> &Vec<T> {
        &self.vertices
    }

    pub fn all_vertices_mut(&mut self) -> &mut Vec<T> {
        &mut self.vertices
    }
}

impl<T: Vertex> std::ops::Index<usize> for VertexVec<T> {
    type Output = T;

    fn index(&self, i: usize) -> &Self::Output {
        if self.indexes.contains(&i) {
            &self.vertices[i]
        } else {
            panic!("Vertex {} does not belong to this VertexVec", i);
        }
    }
}

#[cfg(test)]
mod test_core_graph_trait {
    use super::*;
    use crate::ext::molecule;

    #[test]
    fn test_index() {
        let mol = molecule::Molecule::from_smiles("c1ccccc1CN");
        let vv = VertexVec::init(vec![0, 1, 3, 5], mol.atoms.clone());
        assert_eq!(vv[1].bonds.len(), 2);
    }

    #[test]
    #[should_panic]
    fn test_index_panic() {
        let mol = molecule::Molecule::from_smiles("c1ccccc1CN");
        let vv = VertexVec::init(vec![0, 1, 3, 5], mol.atoms.clone());
        let _ = vv[2].clone();
    }

}
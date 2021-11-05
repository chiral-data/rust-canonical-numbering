// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

mod config;
pub mod orbit_ops;
pub mod graph;
pub mod graph_ops;
pub mod cycle_ops;
pub mod givp;
pub mod cnap;
pub mod reduce;

pub fn symmetry_perception_by_graph_reduction<T: graph::VertexExtendableHash>(
    rg: &mut reduce::ReducibleGraph<T::VertexType>,
    orbits_symmetry: &mut Vec<orbit_ops::Orbit>,
    get_reduced_edges: fn(&Vec<T::VertexType>, &Vec<usize>) -> Vec<(usize, usize, usize)>,
    get_local_symmetric_orbits: fn(&graph::VertexVec<T::VertexType>, &mut Vec<orbit_ops::Orbit>, &mut Vec<orbit_ops::Orbit>),
    custom_maker: usize
) {
    let mut mut_custom_maker = custom_maker;
    reduce::run::<T>(&mut mut_custom_maker, rg, orbits_symmetry, get_reduced_edges, get_local_symmetric_orbits);
}
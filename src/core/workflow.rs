use super::orbit_ops;
use super::graph;
use super::reduce;

pub fn symmetry_perception<T: graph::VertexExtendableHash>(
    rg: &mut reduce::reducible_graph::ReducibleGraph<T::VertexType>,
    orbits_symmetry: &mut Vec<orbit_ops::Orbit>,
    get_reduced_edges: fn(&Vec<T::VertexType>, &Vec<usize>) -> Vec<(usize, usize, usize)>,
    get_local_symmetric_orbits: fn(&graph::VertexVec<T::VertexType>, &mut Vec<orbit_ops::Orbit>, &mut Vec<orbit_ops::Orbit>),
    custom_maker: usize
) {
    let mut mut_custom_maker = custom_maker;
    reduce::run::<T>(&mut mut_custom_maker, rg, orbits_symmetry, get_reduced_edges, get_local_symmetric_orbits);
}
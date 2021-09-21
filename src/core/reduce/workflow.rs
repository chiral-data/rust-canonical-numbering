use crate::core;
use crate::core::cnap::ErrorCNAP;
use super::reducible_graph;
use super::graph_separable;
use super::graph_high_symmetry;


fn get_symmetric_orbits<T: core::graph::VertexExtendableHash>(
    custom_maker: &mut usize,
    rg: &mut reducible_graph::ReducibleGraph<T::VertexType>,
    orbits_symmetry: &mut Vec<core::orbit_ops::Orbit>,
    get_reduced_edges: fn(&Vec<T::VertexType>, &Vec<usize>) -> Vec<(usize, usize, usize)>,
    get_local_symmetric_orbits: fn(&core::graph::VertexVec<T::VertexType>, &mut Vec<core::orbit_ops::Orbit>, &mut Vec<core::orbit_ops::Orbit>),
) {
    // givp for reducible graph
    core::givp::run::<T>(&rg.vv, &mut rg.numbering, &mut rg.orbits_after_partition);
    get_local_symmetric_orbits(&rg.vv, &mut rg.orbits_after_partition, orbits_symmetry);

    // cnap for reducible graph
    match core::cnap::run::<T>(
        &get_reduced_edges(rg.vv.all_vertices(), rg.vv.valid_indexes()),
        rg.vv.all_len(), &rg.orbits_after_partition, orbits_symmetry
    ) {
        Ok(_) => (),
        Err(err) => {
            match err {
                ErrorCNAP::ErrorTwoFolded => { 
                    println!("Case Two-folded Error: {:?}", rg.orbits_after_partition);
                },
                ErrorCNAP::ErrorIncomputable => {
                    match rg.create_mapping() {
                        Ok(_) => {
                            if cfg!(debug_assertions) {
                                println!("Mapping: {:?}", rg.mapping);
                                println!("Boundary Edges: {:?}", rg.boundary_edges);
                            }

                            let mut sg = rg.get_simplified_graph(custom_maker);
                            get_symmetric_orbits::<T>(custom_maker, &mut sg, orbits_symmetry, get_reduced_edges, get_local_symmetric_orbits);
                            let mut fg = rg.get_folded_subgraph(custom_maker, orbits_symmetry);
                            get_symmetric_orbits::<T>(custom_maker, &mut fg, orbits_symmetry, get_reduced_edges, get_local_symmetric_orbits);
                            
                            let mut mapping_symmetric_orbits: Vec<core::orbit_ops::Orbit> = vec![];
                            for idx in 0..rg.mapping[0].len() {
                                mapping_symmetric_orbits.push(
                                    rg.mapping.iter().map(|m| m[idx]).collect()
                                );
                            }
                            orbits_symmetry.append(&mut mapping_symmetric_orbits);
                            core::orbit_ops::orbits_self_merge(orbits_symmetry);
                        },
                        Err(s) => {
                            if cfg!(debug_assertions) {
                                println!("Reducible Graph Mapping Not Succeed: {}", s);
                            }

                            let cycles = graph_separable::find_cycles(&rg, core::config::MAX_REDUCIBLE_CYCLE_SIZE); 
                            if graph_separable::is_separable(&rg, &cycles) {
                                for cycle in cycles.iter() {
                                    let mut new_rg = graph_separable::construct_reducible_graph(&rg, cycle, custom_maker);
                                    get_symmetric_orbits::<T>(custom_maker, &mut new_rg, orbits_symmetry, get_reduced_edges, get_local_symmetric_orbits);
                                }
                            } else {
                                match graph_high_symmetry::run_high_symmetry::<T>(&rg.vv, &get_reduced_edges(rg.vv.all_vertices(), rg.vv.valid_indexes()), rg.vv.all_len(), &rg.orbits_after_partition, orbits_symmetry) {
                                    Ok(_) => (),
                                    Err(s) => {
                                        println!("Reduce CNAP Finally Failed! {:?}", s);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

pub fn run<T: core::graph::VertexExtendableHash>(
    custom_maker: &mut usize,
    rg: &mut reducible_graph::ReducibleGraph<T::VertexType>,
    orbits_symmetry: &mut Vec<core::orbit_ops::Orbit>,
    get_reduced_edges: fn(&Vec<T::VertexType>, &Vec<usize>) -> Vec<(usize, usize, usize)>,
    get_local_symmetric_orbits: fn(&core::graph::VertexVec<T::VertexType>, &mut Vec<core::orbit_ops::Orbit>, &mut Vec<core::orbit_ops::Orbit>),
){
    get_symmetric_orbits::<T>(
        custom_maker, rg, orbits_symmetry, get_reduced_edges, get_local_symmetric_orbits
    );
}

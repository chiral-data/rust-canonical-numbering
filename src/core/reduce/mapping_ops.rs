// Copyr1ight 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::core;

/// Types
pub type VertexMapping = Vec<Vec<usize>>;
pub type NeighbourIndexes = Vec<usize>;
pub type SimpleEdge = (usize, usize);
pub type BoundaryEdges = Vec<Vec<SimpleEdge>>;
pub type MappingStatus = (usize, VertexMapping, BoundaryEdges);

/// Find the neighbours of the current vertex in the mapping sequence;
///     Ignore the neighbours if they have been in the mapping
pub fn find_new_neighbours<T: core::graph::Vertex>(
    vertex_index: usize,
    mapping_index: &Vec<usize>,
    vv: &core::graph::VertexVec<T>
) -> Vec<usize> {
    vv[vertex_index].neighbour_indexes().into_iter()
        .filter(|nb_idx| !mapping_index.contains(&nb_idx))
        .collect()
}

/// Enumerate all possible mappings of the neighbours according to a numbering
fn enumerate_neighbour_mappings(
    neighbours: &NeighbourIndexes, 
    numbering: &Vec<usize>,
) -> Vec<Vec<usize>> {
    let mut neighbour_indexes = neighbours.to_vec();
    neighbour_indexes.sort_by_key(|&ni| numbering[ni]);

    let mut neighbour_indexes_by_group: Vec<Vec<usize>> = vec![vec![neighbour_indexes[0]]];
    let mut current_numbering: usize = numbering[neighbour_indexes[0]];
    for i in 1..neighbour_indexes.len() {
        let size: usize = neighbour_indexes_by_group.len();
        if numbering[neighbour_indexes[i]] == current_numbering {
            neighbour_indexes_by_group[size - 1].push(neighbour_indexes[i]);
        } else {
            current_numbering = numbering[neighbour_indexes[i]];
            neighbour_indexes_by_group.push(vec![neighbour_indexes[i]]);
        }
    }

    let mut permutations: Vec<Vec<usize>> = vec![vec![]];
    for i in 0..neighbour_indexes_by_group.len() {
        if neighbour_indexes_by_group[i].len() == 1 {
            for j in 0..permutations.len() {
                permutations[j].push(neighbour_indexes_by_group[i][0]);
            }
        } else {
            let permutated_neighbours = core::cnap::factorial_vec(&neighbour_indexes_by_group[i]);
            let mut permutations_tmp = permutations.clone();
            permutations.clear();
            while let Some(p) = permutations_tmp.pop() {
                let mut pn_tmp = permutated_neighbours.clone();
                while let Some(mut pn) = pn_tmp.pop() {
                    let mut p_cloned: Vec<usize> = p.clone();
                    p_cloned.append(&mut pn);
                    permutations.push(p_cloned);
                }
            }
        }
    }

    permutations
}

pub enum MappingError {
    MismatchSize,
    MismatchNumbering
}

pub fn mapping_neighbours(
    neighbour_groups: &mut Vec<NeighbourIndexes>,
    numbering: &Vec<usize>,
) -> Result<Vec<VertexMapping>, MappingError> {
    for idx in 1..neighbour_groups.len() {
        if neighbour_groups[0].len() != neighbour_groups[idx].len() {
            // return Err("neighbour counts are not equal".to_string())
            return Err(MappingError::MismatchSize);
        }
    }

    for idx in 0..neighbour_groups.len() {
        neighbour_groups[idx].sort_by_key(|&nb_idx| numbering[nb_idx]);
    }

    for idx in 1..neighbour_groups.len() {
        if neighbour_groups[0].iter().map(|&nb_idx| numbering[nb_idx]).collect::<Vec<usize>>()
            !=  neighbour_groups[idx].iter().map(|&nb_idx| numbering[nb_idx]).collect::<Vec<usize>>() {
            // return Err("mapping numberings are not equal".to_string())
            return Err(MappingError::MismatchNumbering);
        }
    }

    let mut new_mapping_status: Vec<VertexMapping> = vec![vec![neighbour_groups[0].clone()]];
    for i in 1..neighbour_groups.len() {
        let permutations = enumerate_neighbour_mappings(&mut neighbour_groups[i], numbering);
        if permutations.len() == 1 {
            for im in 0..new_mapping_status.len() {
                new_mapping_status[im].push(permutations[0].clone());
            }
        } else {
            let mut mappings_tmp = new_mapping_status.clone();
            new_mapping_status.clear();
            while let Some(m) = mappings_tmp.pop() {
                let mut p_tmp = permutations.clone();
                while let Some(p) = p_tmp.pop() {
                    let mut m_cloned: VertexMapping = m.clone();
                    m_cloned.push(p);
                    new_mapping_status.push(m_cloned);
                }
            } 
        } 
    }

    Ok(new_mapping_status)

    // Is it necessary to double check the edges of the mappings?
}

/// All possible mappings are stored in a mapping stack
/// Update this stack with the mapping result of neighbours
pub fn update_mapping_stack(
    mapping_cur: usize,
    current_mapping: &VertexMapping,
    current_boundary_edges: &BoundaryEdges,
    neighbour_mapping: &Vec<VertexMapping>,
    mapping_stack: &mut Vec<MappingStatus>
) {
    for nm in neighbour_mapping.iter() {
        let mut current_mapping_cloned = current_mapping.to_vec();
        for i in 0..current_mapping_cloned.len() {
            current_mapping_cloned[i].append(&mut nm[i].to_vec())
        }
        mapping_stack.push((mapping_cur, current_mapping_cloned, current_boundary_edges.to_vec()));
    }
}

#[cfg(test)]
mod test_reduce_mapping_ops {
    use crate::ext::molecule;
    use super::*;

    #[test]
    fn test_find_new_neighbours() {
        let smiles_vec: Vec<String> = vec![
            "CC(C)(C)c1cc2c(OCCCCNC(=N)N)c(c1)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)Cc1cc(C(C)(C)C)cc(c1OCCCCNC(=N)N)C2".to_string()
        ].iter().map(|s| s.to_string()).collect();
        let mol = molecule::Molecule::from_smiles(&smiles_vec[0]);
        let vv = core::graph::VertexVec::init((0..mol.atoms.len()).collect(), mol.atoms.clone());

        assert_eq!(find_new_neighbours(7, &vec![], &vv), vec![6, 8, 17]);
        assert_eq!(find_new_neighbours(7, &vec![8], &vv), vec![6, 17]);
        assert_eq!(find_new_neighbours(7, &vec![8, 17], &vv), vec![6]);
    }

    #[test]
    fn test_enumerate_neighbour_mappings() {
        let neighbours: NeighbourIndexes = vec![1, 3, 4];
        let numbering: Vec<usize> = vec![0, 3, 0, 2, 2];
        let all_mappings = enumerate_neighbour_mappings(&neighbours, &numbering);
        assert_eq!(all_mappings.contains(&vec![3, 4, 1]), true);
        assert_eq!(all_mappings.contains(&vec![4, 3, 1]), true);
    }
}
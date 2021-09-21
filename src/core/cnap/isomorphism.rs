use super::permutation;

pub fn is_automorphic(
    permutation: &permutation::Permuation,
    edges: &Vec<(usize, usize, usize)>
) -> bool {
    for (n1, n2, edge_type) in edges.iter() {
        let n1_target: usize = permutation[*n1];
        let n2_target: usize = permutation[*n2];
        if !edges.contains(&(n1_target, n2_target, *edge_type)) {
            return false
        }
    }

    return true
}

#[cfg(test)]
mod test_isomorphism {
    use super::*;

    #[test]
    fn test_is_automorphic() {
        let edges: Vec<(usize, usize, usize)> = vec![
            (0, 1, 9),
            (0, 5, 9),
            (1, 0, 9),
            (1, 2, 9),
            (2, 1, 9),
            (2, 3, 9),
            (3, 2, 9),
            (3, 4, 9),
            (4, 5, 9),
            (4, 3, 9),
            (5, 0, 9),
            (5, 4, 9),
            (5, 6, 1),
            (6, 5, 1),
            (6, 7, 1),
            (7, 6, 1)
        ];

        let mapping_1: Vec<usize> = vec![4, 3, 2, 1, 0, 5, 6, 7];
        assert_eq!(is_automorphic(&mapping_1, &edges), true); 
        let mapping_2: Vec<usize> = vec![4, 2, 3, 1, 0, 5, 6, 7];
        assert_eq!(is_automorphic(&mapping_2, &edges), false); 
    }
}
// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

fn partition_once<T: std::cmp::PartialOrd>(
    indexes: &Vec<usize>,
    values: &Vec<T>
) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let mut part_left: Vec<usize> = vec![]; 
    let mut part_mid: Vec<usize> = vec![]; 
    let mut part_right: Vec<usize> = vec![]; 

    part_mid.push(indexes[0]);
    for i in 1..indexes.len() {
        if values[indexes[i]] == values[indexes[0]] {
            part_mid.push(indexes[i]);
        } else if values[indexes[i]] > values[indexes[0]] {
            part_right.push(indexes[i]);
        } else {
            part_left.push(indexes[i]);
        }
    }

    (part_left, part_mid, part_right)
}


/// Quicksort Algorithm: https://en.wikipedia.org/wiki/Quicksort
pub fn partition_recursively<T: std::cmp::PartialOrd>(
    indexes: &Vec<usize>,
    values: &Vec<T>,
    current_numbering: usize,
    numbering: &mut Vec<usize>,
    orbits: &mut Vec<Vec<usize>>
) {
    let (part_left, part_mid, part_right) = partition_once::<T>(indexes, values);

    if part_mid.len() == 1 {
        numbering[part_mid[0]] = current_numbering - part_right.len();
    } else {
        for idx in part_mid.clone() {
            numbering[idx] = current_numbering - part_right.len();
        }
        orbits.push(part_mid.clone());
    }

    if part_left.len() > 0 {
        let current_numbering_left = current_numbering - part_right.len() - part_mid.len();
        if part_left.len() == 1 {
            numbering[part_left[0]] = current_numbering_left;
        } else {
            partition_recursively::<T>(&part_left, values, current_numbering_left, numbering, orbits);
        }
    }

    if part_right.len() > 0 {
        if part_right.len() > 1 {
            partition_recursively::<T>(&part_right, values, current_numbering, numbering, orbits);
        } else {
            numbering[part_right[0]] = current_numbering;
        }
    }
}

#[cfg(test)]
mod test_givp_partition {
    use super::*;

    #[test]
    fn test_partition_once() {
        let indexes: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let values: Vec<i32>    = vec![9, 8, 9, 4, 12, 2, 4, 11];
        let (part_left, part_mid, part_right) = partition_once::<i32>(&indexes, &values);
        assert_eq!(part_left, vec![1, 3, 5, 6]);
        assert_eq!(part_mid, vec![0, 2]);
        assert_eq!(part_right, vec![4, 7]);
   }

   #[test]
   fn test_partition_recursively() {
       let indexes: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6, 7];
       let values: Vec<i32>    = vec![9, 8, 9, 4, 12, 2, 4, 11];
       let mut numbering: Vec<usize> = vec![values.len(); values.len()];
       let mut orbits: Vec<Vec<usize>> = vec![];
       partition_recursively::<i32>(&indexes, &values, numbering[0], &mut numbering, &mut orbits);

       assert_eq!(numbering, vec![6, 4, 6, 3, 8, 1, 3, 7]);
       assert_eq!(orbits.len(), 2);
       assert!(orbits.contains(&vec![0, 2]));
       assert!(orbits.contains(&vec![3, 6]));
   }
}

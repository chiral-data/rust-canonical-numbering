pub fn factorial(n: usize) -> usize {
    if n > 2 {
        n * factorial(n - 1)
    } else {
        2
    }
}

pub fn factorial_vec(v: &Vec<usize>) -> Vec<Vec<usize>> {
    let mut results: Vec<Vec<usize>> = vec![];

    if v.len() == 2 {
        results.push(v.clone());
        results.push(vec![v[1], v[0]]);
    } else if v.len() > 2 {
        for i in 0..v.len() {
            let mut v_tmp: Vec<usize> = v.clone();
            v_tmp.remove(i);
            let mut results_down: Vec<Vec<usize>> = factorial_vec(&v_tmp);
            while let Some(mut result_tmp) = results_down.pop() {
                result_tmp.push(v[i]);
                results.push(result_tmp);
            }
        }
    }

    results
}


#[cfg(test)]
mod test_combinatorial {
    use super::*;

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(2), 2);
        assert_eq!(factorial(3), 6);
        assert_eq!(factorial(4), 24);
        assert_eq!(factorial(5), 120);
        assert_eq!(factorial(6), 720);
        // assert_eq!(factorial(20), 720); factorial maximun 20, otherwise overflow
    }

    #[test]
    fn test_factorial_vec() {
        let factorial_counts: Vec<usize> = vec![1, 1, 2, 6, 24, 120, 720, 5040, 40320];

        for i in 2..9 {
            let v: Vec<usize> = (0..i).collect();
            let results: Vec<Vec<usize>> = factorial_vec(&v);
            assert_eq!(results.len(), factorial_counts[i]);
        }
    }
}
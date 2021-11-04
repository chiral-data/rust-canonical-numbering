# rust-canonical-numbering
An algorithm of canonical numbering and symmetry perception for graph implemented in Rust

## Examples 

### SMILES canonicalization 
```rust 
use graph_canonicalization;

assert_eq!(graph_canonicalization::ext::molecule::get_canon_smiles(&String::from("c1ccccc1CN")), "NCc1ccccc1".to_string());
```

### Symmetry perception
```rust
use graph_canonicalization;

let atom_vec = graph_canonicalization::ext::molecule::smiles_to_atom_vec("C(C)(C)CCN");
let (orbits_givp, numbering) = graph_canonicalization::ext::molecule::symmetry_perception_givp(&atom_vec);
assert_eq!(orbits_givp, vec![vec![1, 2]]);
assert_eq!(numbering, vec![6, 2, 2, 5, 4, 3]);
let (orbits_cnap) = graph_canonicalization::ext::molecule::symmetry_perception_cnap(&atom_vec, &orbits_givp, &numbering);
assert_eq!(orbits_cnap, vec![vec![1, 2]]);
``` 

## Tests

### Unit tests
```sh
cargo test --lib
```

### Integration tests 
```sh
cargo test --test chembl_test
cargo test --test krotko_test
```
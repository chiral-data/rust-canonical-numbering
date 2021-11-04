# rust-canonical-numbering
An algorithm of canonical numbering and symmetry perception for graph implemented in Rust

## Examples 

### Generate canonical smiles
```rust
ext::molecule::get_canon_smiles(&String::from("c1ccccc1CN"))
```

## Tests

### Unit tests
```sh
cargo test --lib
```

### Integration tests 
```sh
cargo test --test chembl_test
```
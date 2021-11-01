// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub mod core;
pub mod ext;


#[cfg(test)]
mod test_lib {
    use super::*;

    #[test]
    fn test_module_ext() {
        assert_eq!(ext::molecule::get_canon_smiles(&String::from("c1ccccc1CN")), String::from("NCc1ccccc1"))
    }
}
// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub const MAX_COUNT_OF_ATOMS_IN_A_MOLECULE: usize = 1000;
// 1 digit for symmetry
// 3 digits for atomic number
// 1 digit for aromacity
// 1 digit for bond count
// 3 digits for charge
pub const MAX_OF_ATOM_FIXED_HASH_VALUE: usize = 1000 * 10 * 10 * 1000 * 1;
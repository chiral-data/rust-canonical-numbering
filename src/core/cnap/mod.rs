// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! CNAP (Canonical Numbering by Automorphism Permutation)
//!      Symmetry perception with theoretical completeness

mod combinatorial;
mod permutation;
mod isomorphism;
mod workflow;

pub use workflow::ErrorCNAP;
pub use workflow::run;
pub use workflow::is_computable;
pub use workflow::get_symmetric_orbits;
pub use combinatorial::factorial_vec;
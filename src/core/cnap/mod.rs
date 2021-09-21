// CNAP: Canonical Numbering by Automorphism Permutation
// 

mod combinatorial;
mod permutation;
mod isomorphism;
mod workflow;

pub use workflow::ErrorCNAP;
pub use workflow::run;
pub use workflow::is_computable;
pub use workflow::get_symmetric_orbits;
pub use combinatorial::factorial_vec;
// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! GIVP: Graph Invariant Vertex Partitioning
//!     Fast symmetry perception

mod partition;
mod workflow;

pub use workflow::partition_vertices;
pub use workflow::run;
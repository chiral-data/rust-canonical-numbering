/// GIVP: Graph Invariant Vertex Partitioning
/// 

mod partition;
mod workflow;

pub use workflow::partition_vertices;
pub use workflow::partition_vertices_with;
pub use workflow::run;
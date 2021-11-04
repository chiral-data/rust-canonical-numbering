// Copyright 2021 Chiral Ltd.
// Licensed under the Apache-2.0 license (https://opensource.org/licenses/Apache-2.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

/// Computation power for isomorphism checking
/// It can be adjusted according to the hardware.
pub const COMPUTATION_POWER: usize = 4096 * 4;

/// The maximum size of cycles to be indentified
pub const MAX_REDUCIBLE_CYCLE_SIZE: usize = 8;
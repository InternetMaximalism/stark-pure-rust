pub mod r1csfile;
mod reader;

pub use r1csfile::{Coefficient, Constraint, Constraints, Factor, Header, R1csContents, Version};
pub use reader::read_bytes;

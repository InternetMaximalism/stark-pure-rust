use std::fmt::Debug;

pub trait Digest: AsRef<[u8]> + Clone + Default + Sync + Send + Eq + Debug {
    fn hash(message: &[u8]) -> Self;
}

use std::fmt;
use std::error::Error;
use nalgebra;

// Error struct for failed `nalgebra` operations
#[derive(Debug)]
pub struct Kalman2DError{typ: u16}
impl fmt::Display for Kalman2DError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.typ {
            1 => write!(f, "Can inverse matrix"),
            _ => write!(f, "Undefined error")
        }
    }
}
impl Error for Kalman2DError {}
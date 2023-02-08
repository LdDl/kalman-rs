pub mod kalman;
pub mod blob;


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        let result = blob::add(2, 2);
        assert_eq!(result, 4);
    }
}

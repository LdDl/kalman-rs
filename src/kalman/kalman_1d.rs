use std::fmt;
use std::error::Error;
use nalgebra;

// Error struct for failed `nalgebra` operations
#[derive(Debug)]
pub struct Kalman1DError{typ: u16}
impl fmt::Display for Kalman1DError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.typ {
            1 => write!(f, "Can inverse matrix"),
            _ => write!(f, "Undefined error")
        }
    }
}
impl Error for Kalman1DError {}


// Identity matrix. See the ref. https://en.wikipedia.org/wiki/Identity_matrix
const I: nalgebra::SMatrix::<f32, 2, 2> = nalgebra::SMatrix::<f32, 2, 2>::new(
    1.0, 0.0,
    0.0, 1.0,
);

/// Implementation of Discrete Kalman filter for case when there is only on variable X.
#[derive(Debug, Clone)]
pub struct Kalman1D {
    // Single cycle time
    dt: f32,
    // Control input
    u: f32,
    // Standart deviation of acceleration
    std_dev_a: f32,
    // Standart deviation of measurement
    std_dev_m: f32,
    // Transition matrix
    A: nalgebra::SMatrix<f32, 2, 2>,
    // Control matrix
    B: nalgebra::SMatrix<f32, 2, 1>,
    // Transformation (observation) matrix
    H: nalgebra::SMatrix<f32, 1, 2>,
    // Process noise covariance matrix
    Q: nalgebra::SMatrix<f32, 2, 2>,
    // Measurement noise covariance matrix
    R: nalgebra::SMatrix<f32, 1, 1>,
    // Error covariance matrix
    P: nalgebra::SMatrix<f32, 2, 2>,
    // State vector: x, vx
    x: nalgebra::SVector<f32, 2>,
}

impl Kalman1D {
    /// Creates new `Kalman1D`
    /// 
    /// Basic usage:
    /// 
    /// ```
    /// use kalman_rust::kalman::Kalman1D;
    /// let dt = 0.1; // Single cycle time
    /// let u = 2.0; // Control input
    /// let std_dev_a = 0.25; // Standart deviation of acceleration
    /// let std_dev_m = 1.2; // Standart deviation of measurement
    /// let mut kalman = Kalman1D::new(dt, u, std_dev_a, std_dev_m);
    /// ```
    pub fn new(dt: f32, u: f32, std_dev_a: f32, std_dev_m: f32) -> Self {
        Kalman1D {
            dt,
            u,
            std_dev_a,
            std_dev_m,
            // Ref.: Eq.(17)
            A: nalgebra::SMatrix::<f32, 2, 2>::new(
                1.0, dt,
                0.0, 1.0,
            ),
            // Ref.: Eq.(18)
            B: nalgebra::SMatrix::<f32, 2, 1>::new(
                0.5 * dt.powi(2),
                dt,
            ),
            // Ref.: Eq.(20)
            H: nalgebra::SMatrix::<f32, 1, 2>::new(
                1.0, 0.0,
            ),
            // Ref.: Eq.(25)
            Q: nalgebra::SMatrix::<f32, 2, 2>::new(
                0.25 * dt.powi(4), 0.5 * dt.powi(3),
                0.5 * dt.powi(3), dt.powi(2),
            )*std_dev_a.powi(2),
            // Ref.: Eq.(26)
            R: nalgebra::SMatrix::<f32, 1, 1>::new(
                std_dev_m.powi(2),
            ),
            P: nalgebra::SMatrix::<f32, 2, 2>::new(
                1.0, 0.0,
                0.0, 1.0,
            ),
            x: nalgebra::SVector::<f32, 2>::new(
                0.0,
                0.0,
            ),
        }
    }
    /// Projects the state and the error covariance ahead
    /// Mutates the state vector and the error covariance matrix
    /// 
    /// Basic usage:
    /// 
    /// ```
    /// use kalman_rust::kalman::Kalman1D;
    /// let dt = 0.1; // Single cycle time
    /// let u = 2.0; // Control input
    /// let std_dev_a = 0.25; // Standart deviation of acceleration
    /// let std_dev_m = 1.2; // Standart deviation of measurement
    /// let mut kalman = Kalman1D::new(dt, u, std_dev_a, std_dev_m);
    /// let measurements = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    /// for x in measurements.iter() {
    ///     // get measurement
    ///     kalman.predict();
    ///     // then do update 
    /// }
    /// ```
    pub fn predict(&mut self) {
        // Ref.: Eq.(5)
        self.x = (self.A*self.x) + (self.B*self.u);
        // Ref.: Eq.(6)
        self.P = self.A*self.P*self.A.transpose() + self.Q;
    }
    /// Computes the Kalman gain and then updates the state vector and the error covariance matrix
    /// Mutates the state vector and the error covariance matrix.
    /// 
    /// Basic usage:
    /// 
    /// ```
    /// use kalman_rust::kalman::Kalman1D;
    /// let dt = 0.1; // Single cycle time
    /// let u = 2.0; // Control input
    /// let std_dev_a = 0.25; // Standart deviation of acceleration
    /// let std_dev_m = 1.2; // Standart deviation of measurement
    /// let mut kalman = Kalman1D::new(dt, u, std_dev_a, std_dev_m);
    /// let measurements = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    /// for x in measurements {
    ///     kalman.predict();
    ///     kalman.update(x).unwrap(); // assuming that there is noise in measurement
    /// }
    /// ```
    pub fn update(&mut self, _z: f32) -> Result<(), Kalman1DError> {
        // Ref.: Eq.(7)
        let gain = match (self.H*self.P*self.H.transpose() + self.R).try_inverse() {
            Some(inv) => self.P*self.H.transpose()*inv,
            None => return Err(Kalman1DError{typ: 1}),
        };
        // Ref.: Eq.(8)
        let z = nalgebra::SMatrix::<f32, 1, 1>::new(_z);
        let r = z - self.H*self.x;
        // Ref.: Eq.(9)
        self.x = self.x + gain*r;
        // Ref.: Eq.(10)
        self.P = (I - gain*self.H)*self.P;
        Ok(())
    }
    /// Returns the current state (only X, not Vx)
    pub fn get_state(&self) -> f32 {
        self.x[0]
    }
    /// Returns the current state (both X and Vx)
    pub fn get_vector_state(&self) -> nalgebra::SVector::<f32, 2> {
        self.x
    }
}

fn float_loop(start: f32, threshold: f32, step_size: f32) -> impl Iterator<Item = f32> {
    std::iter::successors(Some(start), move |&prev| {
        let next = prev + step_size;
        (next < threshold).then_some(next)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;
    use rand_distr::StandardNormal;
    #[test]
    fn test_1d_kalman() {
        // Just and adoptation of https://machinelearningspace.com/object-tracking-python/
        let dt = 0.1;
        let u = 2.0;
        let std_dev_a = 0.25;
        let std_dev_m = 1.2;

        let t: nalgebra::SVector::<f32, 1000> = nalgebra::SVector::<f32, 1000>::from_iterator(float_loop(0.0, 100.0, dt));
        // let t:  = (0..100).map(|t| t as f32).collect();
        let track = t.map(|t| dt*(t*t - t));

        let mut kalman = Kalman1D::new(dt, u, std_dev_a, std_dev_m);
        let mut measurement: Vec<f32> = vec![];
        let mut predictions: Vec<f32>= vec![];
        for (t, x) in t.iter().zip(track.iter()) {
            // Add some noise to perfect track
            let v: f32 = StdRng::from_os_rng().sample::<f32, StandardNormal>(StandardNormal) * (50.0+50.0) - 50.0; // Generate noise in [-50, 50)
            let z = kalman.H.x * x + v;
            measurement.push(z);

            // Predict stage
            kalman.predict();
            let state = kalman.get_vector_state();
            predictions.push(state.x);

            // Update stage
            kalman.update(z).unwrap();
        }
        // println!("time;perfect;measurement;prediction");
        // for i in 0..track.len() {
        //     println!("{};{};{};{}", t[i], track[i], measurement[i], predictions[i]);
        // }
    }
}

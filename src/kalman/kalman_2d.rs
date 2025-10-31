use nalgebra;
use std::error::Error;
use std::fmt;

// Error struct for failed `nalgebra` operations
#[derive(Debug)]
pub struct Kalman2DError {
    typ: u16,
}
impl fmt::Display for Kalman2DError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.typ {
            1 => write!(f, "Can't inverse matrix"),
            _ => write!(f, "Undefined error"),
        }
    }
}
impl Error for Kalman2DError {}

// Identity matrix. See the ref. https://en.wikipedia.org/wiki/Identity_matrix
const I: nalgebra::SMatrix<f32, 4, 4> = nalgebra::SMatrix::<f32, 4, 4>::new(
    1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
);

/// Implementation of Discrete Kalman filter for case when there are two variables: X and Y.
#[derive(Debug, Clone)]
pub struct Kalman2D {
    // Single cycle time
    dt: f32,
    // Control input
    u: nalgebra::SMatrix<f32, 2, 1>,
    // Standart deviation of acceleration
    std_dev_a: f32,
    // Standart deviation of measurement for X
    std_dev_mx: f32,
    // Standart deviation of measurement for Y
    std_dev_my: f32,
    // Transition matrix
    A: nalgebra::SMatrix<f32, 4, 4>,
    // Control matrix
    B: nalgebra::SMatrix<f32, 4, 2>,
    // Transformation (observation) matrix
    H: nalgebra::SMatrix<f32, 2, 4>,
    // Process noise covariance matrix
    Q: nalgebra::SMatrix<f32, 4, 4>,
    // Measurement noise covariance matrix
    R: nalgebra::SMatrix<f32, 2, 2>,
    // Error covariance matrix
    P: nalgebra::SMatrix<f32, 4, 4>,
    // State vector: x, y, vx, vy
    x: nalgebra::SVector<f32, 4>,
}

impl Kalman2D {
    /// Creates new `Kalman2D`
    ///
    /// Basic usage:
    ///
    /// ```
    /// use kalman_rust::kalman::Kalman2D;
    /// let dt = 0.1; // Single cycle time
    /// let ux = 2.0; // Control input for X
    /// let uy = 2.0; // Control input for Y
    /// let std_dev_a = 0.25; // Standart deviation of acceleration
    /// let std_dev_mx = 1.2; // Standart deviation of measurement for X
    /// let std_dev_my = 1.2; // Standart deviation of measurement for Y
    /// let mut kalman = Kalman2D::new(dt, ux, uy, std_dev_a, std_dev_mx, std_dev_my);
    /// ```
    pub fn new(
        dt: f32,
        ux: f32,
        uy: f32,
        std_dev_a: f32,
        std_dev_mx: f32,
        std_dev_my: f32,
    ) -> Self {
        Kalman2D {
            dt,
            u: nalgebra::SMatrix::<f32, 2, 1>::new(ux, uy),
            std_dev_a,
            std_dev_mx,
            std_dev_my,
            // Ref.: Eq.(31)
            A: nalgebra::SMatrix::<f32, 4, 4>::new(
                1.0, 0.0, dt, 0.0, 0.0, 1.0, 0.0, dt, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ),
            // Ref.: Eq.(32)
            B: nalgebra::SMatrix::<f32, 4, 2>::new(
                0.5 * dt.powi(2),
                0.0,
                0.0,
                0.5 * dt.powi(2),
                dt,
                0.0,
                0.0,
                dt,
            ),
            // Ref.: Eq.(34)
            H: nalgebra::SMatrix::<f32, 2, 4>::new(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            // Ref.: Eq.(40)
            Q: nalgebra::SMatrix::<f32, 4, 4>::new(
                0.25 * dt.powi(4),
                0.0,
                0.5 * dt.powi(3),
                0.0,
                0.0,
                0.25 * dt.powi(4),
                0.0,
                0.5 * dt.powi(3),
                0.5 * dt.powi(3),
                0.0,
                dt.powi(2),
                0.0,
                0.0,
                0.5 * dt.powi(3),
                0.0,
                dt.powi(2),
            ) * std_dev_a.powi(2),
            // Ref.: Eq.(41)
            R: nalgebra::SMatrix::<f32, 2, 2>::new(
                std_dev_mx.powi(2),
                0.0,
                0.0,
                std_dev_my.powi(2),
            ),
            P: nalgebra::SMatrix::<f32, 4, 4>::new(
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ),
            x: nalgebra::SVector::<f32, 4>::new(0.0, 0.0, 0.0, 0.0),
        }
    }
    /// Creates new `Kalman2D` with initial state
    ///
    /// Why is it needed to set the initial state to the actual first observed coordinates of an object (sometimes)?
    /// When the first state vector is initialized with zeros, it assumes that the object is at the origin
    /// and the filter needs to estimate the position of the object from scratch, which can result in some initial inaccuracies.
    /// On the other hand, initializing the first state vector with the actual observed coordinates of the object can provide
    /// a more accurate estimate from the beginning, which can improve the overall tracking performance of the filter
    ///
    ///
    /// Basic usage:
    ///
    /// ```
    /// use kalman_rust::kalman::Kalman2D;
    /// let dt = 0.1; // Single cycle time
    /// let ux = 2.0; // Control input for X
    /// let uy = 2.0; // Control input for Y
    /// let std_dev_a = 0.25; // Standart deviation of acceleration
    /// let std_dev_mx = 1.2; // Standart deviation of measurement for X
    /// let std_dev_my = 1.2; // Standart deviation of measurement for Y
    /// let ix = 1.0; // Initial state for X
    /// let iy = 5.0; // Initial state for Y
    /// let mut kalman = Kalman2D::new_with_state(dt, ux, uy, std_dev_a, std_dev_mx, std_dev_my, ix, iy);
    /// ```
    pub fn new_with_state(
        dt: f32,
        ux: f32,
        uy: f32,
        std_dev_a: f32,
        std_dev_mx: f32,
        std_dev_my: f32,
        x: f32,
        y: f32,
    ) -> Self {
        Kalman2D {
            dt,
            u: nalgebra::SMatrix::<f32, 2, 1>::new(ux, uy),
            std_dev_a,
            std_dev_mx,
            std_dev_my,
            // Ref.: Eq.(31)
            A: nalgebra::SMatrix::<f32, 4, 4>::new(
                1.0, 0.0, dt, 0.0, 0.0, 1.0, 0.0, dt, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ),
            // Ref.: Eq.(32)
            B: nalgebra::SMatrix::<f32, 4, 2>::new(
                0.5 * dt.powi(2),
                0.0,
                0.0,
                0.5 * dt.powi(2),
                dt,
                0.0,
                0.0,
                dt,
            ),
            // Ref.: Eq.(34)
            H: nalgebra::SMatrix::<f32, 2, 4>::new(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            // Ref.: Eq.(40)
            Q: nalgebra::SMatrix::<f32, 4, 4>::new(
                0.25 * dt.powi(4),
                0.0,
                0.5 * dt.powi(3),
                0.0,
                0.0,
                0.25 * dt.powi(4),
                0.0,
                0.5 * dt.powi(3),
                0.5 * dt.powi(3),
                0.0,
                dt.powi(2),
                0.0,
                0.0,
                0.5 * dt.powi(3),
                0.0,
                dt.powi(2),
            ) * std_dev_a.powi(2),
            // Ref.: Eq.(41)
            R: nalgebra::SMatrix::<f32, 2, 2>::new(
                std_dev_mx.powi(2),
                0.0,
                0.0,
                std_dev_my.powi(2),
            ),
            P: nalgebra::SMatrix::<f32, 4, 4>::new(
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ),
            x: nalgebra::SVector::<f32, 4>::new(x, y, 0.0, 0.0),
        }
    }
    /// Projects the state and the error covariance ahead
    /// Mutates the state vector and the error covariance matrix
    ///
    /// Basic usage:
    ///
    /// ```
    /// use kalman_rust::kalman::Kalman2D;
    /// let dt = 0.1; // Single cycle time
    /// let ux = 2.0; // Control input for X
    /// let uy = 2.0; // Control input for Y
    /// let std_dev_a = 0.25; // Standart deviation of acceleration
    /// let std_dev_mx = 1.2; // Standart deviation of measurement for X
    /// let std_dev_my = 1.2; // Standart deviation of measurement for Y
    /// let mut kalman = Kalman2D::new(dt, ux, uy, std_dev_a, std_dev_mx, std_dev_my);
    /// let measurements = vec![(1.0, 2.0), (2.0, 3.0), (3.0, 4.0)];
    /// for x in measurements.iter() {
    ///     // get measurement
    ///     kalman.predict();
    ///     // then do update
    /// }
    /// ```
    pub fn predict(&mut self) {
        // Ref.: Eq.(5)
        self.x = (self.A * self.x) + (self.B * self.u);
        // Ref.: Eq.(6)
        self.P = self.A * self.P * self.A.transpose() + self.Q;
    }
    /// Computes the Kalman gain and then updates the state vector and the error covariance matrix
    /// Mutates the state vector and the error covariance matrix.
    ///
    /// Basic usage:
    ///
    /// ```
    /// use kalman_rust::kalman::Kalman2D;
    /// let dt = 0.1; // Single cycle time
    /// let ux = 2.0; // Control input for X
    /// let uy = 2.0; // Control input for Y
    /// let std_dev_a = 0.25; // Standart deviation of acceleration
    /// let std_dev_mx = 1.2; // Standart deviation of measurement for X
    /// let std_dev_my = 1.2; // Standart deviation of measurement for Y
    /// let mut kalman = Kalman2D::new(dt, ux, uy, std_dev_a, std_dev_mx, std_dev_my);
    /// let measurements = vec![(1.0, 2.0), (2.0, 3.0), (3.0, 4.0)];
    /// for (mx, my) in measurements.iter() {
    ///     kalman.predict();
    ///     kalman.update(*mx, *my).unwrap(); // assuming that there is noise in measurement
    /// }
    /// ```
    pub fn update(&mut self, _zx: f32, _zy: f32) -> Result<(), Kalman2DError> {
        // Ref.: Eq.(7)
        let gain = match (self.H * self.P * self.H.transpose() + self.R).try_inverse() {
            Some(inv) => self.P * self.H.transpose() * inv,
            None => return Err(Kalman2DError { typ: 1 }),
        };
        // Ref.: Eq.(8)
        let z = nalgebra::SMatrix::<f32, 2, 1>::new(_zx, _zy);
        let r = z - self.H * self.x;
        // Ref.: Eq.(9)
        self.x = self.x + gain * r;
        // Ref.: Eq.(10)
        self.P = (I - gain * self.H) * self.P;
        Ok(())
    }
    /// Returns the current state (only X and Y, not Vx and Vy)
    pub fn get_state(&self) -> (f32, f32) {
        (self.x[0], self.x[1])
    }
    /// Returns the current velocity (only Vx and Vy, not X and Y)
    pub fn get_velocity(&self) -> (f32, f32) {
        (self.x[2], self.x[3])
    }
    /// Returns prediction without mutating the state vector and the error covariance matrix
    pub fn get_predicted_position(&self) -> (f32, f32) {
        let x_pred = (self.A * self.x) + (self.B * self.u);
        (x_pred[0], x_pred[1])
    }
    /// Returns position uncertainty from P matrix
    pub fn get_position_uncertainty(&self) -> f32 {
        (self.P[(0, 0)].powi(2) + self.P[(1, 1)].powi(2)).sqrt()
    }
    /// Returns the current state (both (X, Y) and (Vx, Vy))
    pub fn get_vector_state(&self) -> nalgebra::SVector<f32, 4> {
        self.x
    }
}

mod tests {
    use super::*;
    #[test]
    fn test_2d_kalman() {
        let dt = 0.04; // 1/25 = 25 fps - just an example
        let ux = 1.0;
        let uy = 1.0;
        let std_dev_a = 2.0;
        let std_dev_mx = 0.1;
        let std_dev_my = 0.1;

        // Sample measurements
        // Note: in this example Y-axis going from up to down
        let xs = vec![
            311, 312, 313, 311, 311, 312, 312, 313, 312, 312, 312, 312, 312, 312, 312, 312, 312,
            312, 311, 311, 311, 311, 311, 310, 311, 311, 311, 310, 310, 308, 307, 308, 308, 308,
            307, 307, 307, 308, 307, 307, 307, 307, 307, 308, 307, 309, 306, 307, 306, 307, 308,
            306, 306, 306, 305, 307, 307, 307, 306, 306, 306, 307, 307, 308, 307, 307, 308, 307,
            306, 308, 309, 309, 309, 309, 308, 309, 309, 309, 308, 311, 311, 307, 311, 307, 313,
            311, 307, 311, 311, 306, 312, 312, 312, 312, 312, 312, 312, 312, 312, 312, 312, 312,
            312, 312, 312, 312, 312, 312, 312, 312, 312, 312,
        ];
        let ys = vec![
            5, 6, 8, 10, 11, 12, 12, 13, 16, 16, 18, 18, 19, 19, 20, 20, 22, 22, 23, 23, 24, 24,
            28, 30, 32, 35, 39, 42, 44, 46, 56, 58, 70, 60, 52, 64, 51, 70, 70, 70, 66, 83, 80, 85,
            80, 98, 79, 98, 61, 94, 101, 94, 104, 94, 107, 112, 108, 108, 109, 109, 121, 108, 108,
            120, 122, 122, 128, 130, 122, 140, 122, 122, 140, 122, 134, 141, 136, 136, 154, 155,
            155, 150, 161, 162, 169, 171, 181, 175, 175, 163, 178, 178, 178, 178, 178, 178, 178,
            178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178,
        ];

        // Assume that initial X,Y coordinates match the first measurement
        let ix = xs[0] as f32; // Initial state for X
        let iy = ys[0] as f32; // Initial state for Y
        let mut kalman = Kalman2D::new_with_state(dt, ux, uy, std_dev_a, std_dev_mx, std_dev_my, ix, iy);
        kalman.x.x = xs[0] as f32;
        kalman.x.y = ys[0] as f32;
        let mut predictions: Vec<Vec<f32>> = vec![];
        let mut updated_states: Vec<Vec<f32>> = vec![];
        for (x, y) in xs.iter().zip(ys.iter()) {
            // Considering that the measurements are noisy
            let mx = *x as f32;
            let my = *y as f32;

            // Predict stage
            kalman.predict();
            let state = kalman.get_vector_state();
            predictions.push(vec![state.x, state.y]);

            // Update stage
            kalman.update(mx, my).unwrap();
            let updated_state = kalman.get_vector_state();
            updated_states.push(vec![updated_state.x, updated_state.y]);
        }

        // println!("measurement X;measurement Y;prediction X;prediction Y;updated X;updated Y");
        // for i in 0..xs.len() {
        //     println!("{};{};{};{};{};{}", xs[i], ys[i], predictions[i][0], predictions[i][1], updated_states[i][0], updated_states[i][1]);
        // }
    }
}

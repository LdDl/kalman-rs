use nalgebra;
use std::error::Error;
use std::fmt;

// Error struct for failed `nalgebra` operations
#[derive(Debug)]
pub struct KalmanBBoxError {
    typ: u16,
}
impl fmt::Display for KalmanBBoxError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.typ {
            1 => write!(f, "Can't inverse matrix"),
            _ => write!(f, "Undefined error"),
        }
    }
}
impl Error for KalmanBBoxError {}

/// Implementation of Discrete Kalman filter for bounding box tracking.
/// State vector contains: center_x, center_y, width, height, velocity_cx, velocity_cy, velocity_w, velocity_h
#[derive(Debug, Clone)]
pub struct KalmanBBox {
    // Single cycle time
    dt: f32,
    // Control input (accelerations for cx, cy, w, h)
    u: nalgebra::SMatrix<f32, 4, 1>,
    // Standart deviation of acceleration
    std_dev_a: f32,
    // Standart deviation of measurement for center X
    std_dev_m_cx: f32,
    // Standart deviation of measurement for center Y
    std_dev_m_cy: f32,
    // Standart deviation of measurement for width
    std_dev_m_w: f32,
    // Standart deviation of measurement for height
    std_dev_m_h: f32,
    // Transition matrix
    A: nalgebra::SMatrix<f32, 8, 8>,
    // Control matrix
    B: nalgebra::SMatrix<f32, 8, 4>,
    // Transformation (observation) matrix
    H: nalgebra::SMatrix<f32, 4, 8>,
    // Process noise covariance matrix
    Q: nalgebra::SMatrix<f32, 8, 8>,
    // Measurement noise covariance matrix
    R: nalgebra::SMatrix<f32, 4, 4>,
    // Error covariance matrix
    P: nalgebra::SMatrix<f32, 8, 8>,
    // State vector: cx, cy, w, h, vx, vy, vw, vh
    x: nalgebra::SVector<f32, 8>,
}

impl KalmanBBox {
    /// Creates new `KalmanBBox`
    ///
    /// Basic usage:
    ///
    /// ```
    /// use kalman_rust::kalman::KalmanBBox;
    /// let dt = 0.04; // Single cycle time (1/25 fps)
    /// let u_cx = 1.0; // Control input for center X
    /// let u_cy = 1.0; // Control input for center Y
    /// let u_w = 0.0; // Control input for width
    /// let u_h = 0.0; // Control input for height
    /// let std_dev_a = 2.0; // Standart deviation of acceleration
    /// let std_dev_m_cx = 0.1; // Standart deviation of measurement for center X
    /// let std_dev_m_cy = 0.1; // Standart deviation of measurement for center Y
    /// let std_dev_m_w = 0.1; // Standart deviation of measurement for width
    /// let std_dev_m_h = 0.1; // Standart deviation of measurement for height
    /// let mut kalman = KalmanBBox::new(dt, u_cx, u_cy, u_w, u_h, std_dev_a, std_dev_m_cx, std_dev_m_cy, std_dev_m_w, std_dev_m_h);
    /// ```
    pub fn new(
        dt: f32,
        u_cx: f32,
        u_cy: f32,
        u_w: f32,
        u_h: f32,
        std_dev_a: f32,
        std_dev_m_cx: f32,
        std_dev_m_cy: f32,
        std_dev_m_w: f32,
        std_dev_m_h: f32,
    ) -> Self {
        let dt2 = dt * dt;
        let dt3 = dt2 * dt;
        let dt4 = dt3 * dt;
        let std_dev_a2 = std_dev_a * std_dev_a;

        KalmanBBox {
            dt,
            u: nalgebra::SMatrix::<f32, 4, 1>::new(u_cx, u_cy, u_w, u_h),
            std_dev_a,
            std_dev_m_cx,
            std_dev_m_cy,
            std_dev_m_w,
            std_dev_m_h,
            // Ref.: Eq.(53)
            #[rustfmt::skip]
            A: nalgebra::SMatrix::<f32, 8, 8>::from_row_slice(&[
                1.0, 0.0, 0.0, 0.0, dt, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0, 0.0, dt, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, dt, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, dt,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ]),
            // Ref.: Eq.(54)
            #[rustfmt::skip]
            B: nalgebra::SMatrix::<f32, 8, 4>::from_row_slice(&[
                0.5 * dt2, 0.0, 0.0, 0.0,
                0.0, 0.5 * dt2, 0.0, 0.0,
                0.0, 0.0, 0.5 * dt2, 0.0,
                0.0, 0.0, 0.0, 0.5 * dt2,
                dt, 0.0, 0.0, 0.0,
                0.0, dt, 0.0, 0.0,
                0.0, 0.0, dt, 0.0,
                0.0, 0.0, 0.0, dt,
            ]),
            // Ref.: Eq.(56)
            #[rustfmt::skip]
            H: nalgebra::SMatrix::<f32, 4, 8>::from_row_slice(&[
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            ]),
            // Ref.: Eq.(62)
            #[rustfmt::skip]
            Q: nalgebra::SMatrix::<f32, 8, 8>::from_row_slice(&[
                0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0,
                0.0, 0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0,
                0.0, 0.0, 0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0,
                0.0, 0.0, 0.0, 0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2,
                0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2, 0.0, 0.0, 0.0,
                0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2, 0.0, 0.0,
                0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2, 0.0,
                0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2,
            ]),
            // Ref.: Eq.(63)
            R: nalgebra::SMatrix::<f32, 4, 4>::new(
                std_dev_m_cx * std_dev_m_cx, 0.0, 0.0, 0.0,
                0.0, std_dev_m_cy * std_dev_m_cy, 0.0, 0.0,
                0.0, 0.0, std_dev_m_w * std_dev_m_w, 0.0,
                0.0, 0.0, 0.0, std_dev_m_h * std_dev_m_h,
            ),
            #[rustfmt::skip]
            P: nalgebra::SMatrix::<f32, 8, 8>::from_row_slice(&[
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ]),
            #[rustfmt::skip]
            x: nalgebra::SVector::<f32, 8>::from_row_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        }
    }

    /// Creates new `KalmanBBox` with initial state
    ///
    /// Why is it needed to set the initial state to the actual first observed bounding box (sometimes)?
    /// When the first state vector is initialized with zeros, it assumes that the object is at the origin
    /// and the filter needs to estimate the position of the object from scratch, which can result in some initial inaccuracies.
    /// On the other hand, initializing the first state vector with the actual observed bounding box can provide
    /// a more accurate estimate from the beginning, which can improve the overall tracking performance of the filter
    ///
    /// Basic usage:
    ///
    /// ```
    /// use kalman_rust::kalman::KalmanBBox;
    /// let dt = 0.04; // Single cycle time (1/25 fps)
    /// let u_cx = 1.0; // Control input for center X
    /// let u_cy = 1.0; // Control input for center Y
    /// let u_w = 0.0; // Control input for width
    /// let u_h = 0.0; // Control input for height
    /// let std_dev_a = 2.0; // Standart deviation of acceleration
    /// let std_dev_m_cx = 0.1; // Standart deviation of measurement for center X
    /// let std_dev_m_cy = 0.1; // Standart deviation of measurement for center Y
    /// let std_dev_m_w = 0.1; // Standart deviation of measurement for width
    /// let std_dev_m_h = 0.1; // Standart deviation of measurement for height
    /// let i_cx = 100.0; // Initial center X
    /// let i_cy = 50.0; // Initial center Y
    /// let i_w = 40.0; // Initial width
    /// let i_h = 80.0; // Initial height
    /// let mut kalman = KalmanBBox::new_with_state(dt, u_cx, u_cy, u_w, u_h, std_dev_a, std_dev_m_cx, std_dev_m_cy, std_dev_m_w, std_dev_m_h, i_cx, i_cy, i_w, i_h);
    /// ```
    pub fn new_with_state(
        dt: f32,
        u_cx: f32,
        u_cy: f32,
        u_w: f32,
        u_h: f32,
        std_dev_a: f32,
        std_dev_m_cx: f32,
        std_dev_m_cy: f32,
        std_dev_m_w: f32,
        std_dev_m_h: f32,
        cx: f32,
        cy: f32,
        w: f32,
        h: f32,
    ) -> Self {
        let dt2 = dt * dt;
        let dt3 = dt2 * dt;
        let dt4 = dt3 * dt;
        let std_dev_a2 = std_dev_a * std_dev_a;

        KalmanBBox {
            dt,
            u: nalgebra::SMatrix::<f32, 4, 1>::new(u_cx, u_cy, u_w, u_h),
            std_dev_a,
            std_dev_m_cx,
            std_dev_m_cy,
            std_dev_m_w,
            std_dev_m_h,
            // Ref.: Eq.(53)
            #[rustfmt::skip]
            A: nalgebra::SMatrix::<f32, 8, 8>::from_row_slice(&[
                1.0, 0.0, 0.0, 0.0, dt, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0, 0.0, dt, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, dt, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, dt,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ]),
            // Ref.: Eq.(54)
            #[rustfmt::skip]
            B: nalgebra::SMatrix::<f32, 8, 4>::from_row_slice(&[
                0.5 * dt2, 0.0, 0.0, 0.0,
                0.0, 0.5 * dt2, 0.0, 0.0,
                0.0, 0.0, 0.5 * dt2, 0.0,
                0.0, 0.0, 0.0, 0.5 * dt2,
                dt, 0.0, 0.0, 0.0,
                0.0, dt, 0.0, 0.0,
                0.0, 0.0, dt, 0.0,
                0.0, 0.0, 0.0, dt,
            ]),
            // Ref.: Eq.(56)
            #[rustfmt::skip]
            H: nalgebra::SMatrix::<f32, 4, 8>::from_row_slice(&[
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            ]),
            // Ref.: Eq.(62)
            #[rustfmt::skip]
            Q: nalgebra::SMatrix::<f32, 8, 8>::from_row_slice(&[
                0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0,
                0.0, 0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0,
                0.0, 0.0, 0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0,
                0.0, 0.0, 0.0, 0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2,
                0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2, 0.0, 0.0, 0.0,
                0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2, 0.0, 0.0,
                0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2, 0.0,
                0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2,
            ]),
            // Ref.: Eq.(63)
            R: nalgebra::SMatrix::<f32, 4, 4>::new(
                std_dev_m_cx * std_dev_m_cx, 0.0, 0.0, 0.0,
                0.0, std_dev_m_cy * std_dev_m_cy, 0.0, 0.0,
                0.0, 0.0, std_dev_m_w * std_dev_m_w, 0.0,
                0.0, 0.0, 0.0, std_dev_m_h * std_dev_m_h,
            ),
            #[rustfmt::skip]
            P: nalgebra::SMatrix::<f32, 8, 8>::from_row_slice(&[
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ]),
            #[rustfmt::skip]
            x: nalgebra::SVector::<f32, 8>::from_row_slice(&[cx, cy, w, h, 0.0, 0.0, 0.0, 0.0]),
        }
    }

    /// Projects the state and the error covariance ahead
    /// Mutates the state vector and the error covariance matrix
    ///
    /// Basic usage:
    ///
    /// ```
    /// use kalman_rust::kalman::KalmanBBox;
    /// let dt = 0.04;
    /// let mut kalman = KalmanBBox::new(dt, 1.0, 1.0, 0.0, 0.0, 2.0, 0.1, 0.1, 0.1, 0.1);
    /// let measurements = vec![(100.0, 50.0, 40.0, 80.0), (102.0, 52.0, 41.0, 81.0)];
    /// for _ in measurements.iter() {
    ///     // get measurement
    ///     kalman.predict();
    ///     // then do update
    /// }
    /// ```
    /// Rebuilds the time-dependent matrices for a new cycle time.
    ///
    /// See `Kalman2D::set_dt` for the reasoning: a discrete filter is only valid
    /// for the sampling interval its `A` and `Q` were built for, so an uneven
    /// frame arrival has to be fed in here before `predict`.
    pub fn set_dt(&mut self, dt: f32) {
        if dt <= 0.0 || (dt - self.dt).abs() < f32::EPSILON {
            return;
        }
        self.dt = dt;
        let dt2 = dt * dt;
        let dt3 = dt2 * dt;
        let dt4 = dt3 * dt;
        let std_dev_a2 = self.std_dev_a * self.std_dev_a;
        // Ref.: Eq.(53)
        self.A = nalgebra::SMatrix::<f32, 8, 8>::from_row_slice(&[
            1.0, 0.0, 0.0, 0.0, dt, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, dt, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, dt, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, dt,
            0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ]);
        // Ref.: Eq.(54)
        self.B = nalgebra::SMatrix::<f32, 8, 4>::from_row_slice(&[
            0.5 * dt2, 0.0, 0.0, 0.0,
            0.0, 0.5 * dt2, 0.0, 0.0,
            0.0, 0.0, 0.5 * dt2, 0.0,
            0.0, 0.0, 0.0, 0.5 * dt2,
            dt, 0.0, 0.0, 0.0,
            0.0, dt, 0.0, 0.0,
            0.0, 0.0, dt, 0.0,
            0.0, 0.0, 0.0, dt,
        ]);
        // Ref.: Eq.(62)
        self.Q = nalgebra::SMatrix::<f32, 8, 8>::from_row_slice(&[
            0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0,
            0.0, 0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0,
            0.0, 0.0, 0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0,
            0.0, 0.0, 0.0, 0.25 * dt4 * std_dev_a2, 0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2,
            0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2, 0.0, 0.0, 0.0,
            0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2, 0.0, 0.0,
            0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2, 0.0,
            0.0, 0.0, 0.0, 0.5 * dt3 * std_dev_a2, 0.0, 0.0, 0.0, dt2 * std_dev_a2,
        ]);
    }
    /// Returns the cycle time the filter is currently built for
    pub fn get_dt(&self) -> f32 {
        self.dt
    }
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
    /// use kalman_rust::kalman::KalmanBBox;
    /// let dt = 0.04;
    /// let mut kalman = KalmanBBox::new(dt, 1.0, 1.0, 0.0, 0.0, 2.0, 0.1, 0.1, 0.1, 0.1);
    /// let measurements = vec![(100.0, 50.0, 40.0, 80.0), (102.0, 52.0, 41.0, 81.0)];
    /// for (cx, cy, w, h) in measurements.iter() {
    ///     kalman.predict();
    ///     kalman.update(*cx, *cy, *w, *h).unwrap();
    /// }
    /// ```
    pub fn update(&mut self, _z_cx: f32, _z_cy: f32, _z_w: f32, _z_h: f32) -> Result<(), KalmanBBoxError> {
        // Ref.: Eq.(7)
        let S = self.H * self.P * self.H.transpose() + self.R;
        let gain = match S.try_inverse() {
            Some(inv) => self.P * self.H.transpose() * inv,
            None => return Err(KalmanBBoxError { typ: 1 }),
        };
        // Ref.: Eq.(8)
        let z = nalgebra::SMatrix::<f32, 4, 1>::new(_z_cx, _z_cy, _z_w, _z_h);
        let r = z - self.H * self.x;
        // Ref.: Eq.(9)
        self.x = self.x + gain * r;
        // Ref.: Eq.(10)
        // Identity matrix for 8x8
        #[rustfmt::skip]
        let I: nalgebra::SMatrix<f32, 8, 8> = nalgebra::SMatrix::<f32, 8, 8>::from_row_slice(&[
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ]);
        self.P = (I - gain * self.H) * self.P;
        Ok(())
    }

    /// Returns the current state (only cx, cy, w, h - not velocities)
    pub fn get_state(&self) -> (f32, f32, f32, f32) {
        (self.x[0], self.x[1], self.x[2], self.x[3])
    }

    /// Returns the current velocity (only vx, vy, vw, vh)
    pub fn get_velocity(&self) -> (f32, f32, f32, f32) {
        (self.x[4], self.x[5], self.x[6], self.x[7])
    }

    /// Returns prediction without mutating the state vector and the error covariance matrix
    pub fn get_predicted_state(&self) -> (f32, f32, f32, f32) {
        let x_pred = (self.A * self.x) + (self.B * self.u);
        (x_pred[0], x_pred[1], x_pred[2], x_pred[3])
    }

    /// Returns position uncertainty from P matrix (for center coordinates)
    pub fn get_position_uncertainty(&self) -> f32 {
        (self.P[(0, 0)].powi(2) + self.P[(1, 1)].powi(2)).sqrt()
    }

    /// Returns the current state (full 8-element vector)
    pub fn get_vector_state(&self) -> nalgebra::SVector<f32, 8> {
        self.x
    }

    /// Computes the squared Mahalanobis distance between measurement and predicted state
    /// Ref.: Eq.(66)
    ///
    /// This is useful for data association in multi-object tracking.
    /// A detection passes the gate if mahalanobis_distance_squared < threshold
    /// Common threshold for 4 DOF at 95% confidence: 9.488
    pub fn mahalanobis_distance_squared(&self, z_cx: f32, z_cy: f32, z_w: f32, z_h: f32) -> Result<f32, KalmanBBoxError> {
        // Ref.: Eq.(64) - Innovation (residual)
        let z = nalgebra::SMatrix::<f32, 4, 1>::new(z_cx, z_cy, z_w, z_h);
        let y = z - self.H * self.x;

        // Ref.: Eq.(65) - Innovation covariance
        let S = self.H * self.P * self.H.transpose() + self.R;

        // Ref.: Eq.(66) - Squared Mahalanobis distance
        let S_inv = match S.try_inverse() {
            Some(inv) => inv,
            None => return Err(KalmanBBoxError { typ: 1 }),
        };

        let d_squared = (y.transpose() * S_inv * y)[(0, 0)];
        Ok(d_squared)
    }

    /// Computes the Mahalanobis distance (square root of squared distance)
    /// Ref.: Eq.(67)
    pub fn mahalanobis_distance(&self, z_cx: f32, z_cy: f32, z_w: f32, z_h: f32) -> Result<f32, KalmanBBoxError> {
        let d_squared = self.mahalanobis_distance_squared(z_cx, z_cy, z_w, z_h)?;
        Ok(d_squared.sqrt())
    }

    /// Returns the innovation covariance matrix S
    /// Ref.: Eq.(65)
    ///
    /// S = H * P * H^T + R
    pub fn get_innovation_covariance(&self) -> nalgebra::SMatrix<f32, 4, 4> {
        self.H * self.P * self.H.transpose() + self.R
    }

    /// Checks if a measurement passes the gating threshold
    /// Ref.: Eq.(70)
    ///
    /// Common thresholds for alpha = 0.95:
    /// - 4 DOF (bounding box): 9.488
    /// - 2 DOF (position only): 5.991
    pub fn gating_check(&self, z_cx: f32, z_cy: f32, z_w: f32, z_h: f32, threshold: f32) -> bool {
        match self.mahalanobis_distance_squared(z_cx, z_cy, z_w, z_h) {
            Ok(d_sq) => d_sq < threshold,
            Err(_) => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bbox_kalman() {
        let dt = 0.04; // 1/25 = 25 fps - just an example
        let u_cx = 1.0;
        let u_cy = 1.0;
        let u_w = 0.0; // No control input for size
        let u_h = 0.0;
        let std_dev_a = 2.0;
        let std_dev_m_cx = 0.1;
        let std_dev_m_cy = 0.1;
        let std_dev_m_w = 0.1;
        let std_dev_m_h = 0.1;

        // Sample measurements for center coordinates (same as test_2d_kalman)
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

        // Generate bounding boxes around centroids
        // Base width/height with some variation
        let base_w = 40.0_f32;
        let base_h = 80.0_f32;

        // Generate width/height with slight variations (simulating object size changes)
        let ws: Vec<f32> = xs.iter().enumerate().map(|(i, _)| {
            // Slight oscillation
            base_w + (i as f32 * 0.1).sin() * 5.0
        }).collect();

        let hs: Vec<f32> = ys.iter().enumerate().map(|(i, _)| {
            // Slight oscillation
            base_h + (i as f32 * 0.15).sin() * 8.0
        }).collect();

        // Assume that initial state matches the first measurement
        let i_cx = xs[0] as f32;
        let i_cy = ys[0] as f32;
        let i_w = ws[0];
        let i_h = hs[0];

        let mut kalman = KalmanBBox::new_with_state(
            dt, u_cx, u_cy, u_w, u_h,
            std_dev_a,
            std_dev_m_cx, std_dev_m_cy, std_dev_m_w, std_dev_m_h,
            i_cx, i_cy, i_w, i_h
        );

        let mut predictions: Vec<(f32, f32, f32, f32)> = vec![];
        let mut updated_states: Vec<(f32, f32, f32, f32)> = vec![];

        for i in 0..xs.len() {
            let m_cx = xs[i] as f32;
            let m_cy = ys[i] as f32;
            let m_w = ws[i];
            let m_h = hs[i];

            // Predict stage
            kalman.predict();
            predictions.push(kalman.get_state());

            // Update stage
            kalman.update(m_cx, m_cy, m_w, m_h).unwrap();
            updated_states.push(kalman.get_state());
        }

        // Verify that predictions and updates are reasonable
        assert!(predictions.len() == xs.len());
        assert!(updated_states.len() == xs.len());

        // Check that final state is close to final measurement
        let (final_cx, final_cy, final_w, final_h) = updated_states.last().unwrap();
        let last_m_cx = *xs.last().unwrap() as f32;
        let last_m_cy = *ys.last().unwrap() as f32;
        let last_m_w = *ws.last().unwrap();
        let last_m_h = *hs.last().unwrap();

        assert!((final_cx - last_m_cx).abs() < 5.0, "Final cx too far from measurement");
        assert!((final_cy - last_m_cy).abs() < 5.0, "Final cy too far from measurement");
        assert!((final_w - last_m_w).abs() < 5.0, "Final w too far from measurement");
        assert!((final_h - last_m_h).abs() < 5.0, "Final h too far from measurement");

        // println!("measurement cx;measurement cy;measurement w;measurement h;prediction cx;prediction cy;prediction w;prediction h;updated cx;updated cy;updated w;updated h");
        // for i in 0..xs.len() {
        //     let (p_cx, p_cy, p_w, p_h) = predictions[i];
        //     let (u_cx, u_cy, u_w, u_h) = updated_states[i];
        //     println!("{};{};{};{};{};{};{};{};{};{};{};{}",
        //         xs[i], ys[i], ws[i], hs[i],
        //         p_cx, p_cy, p_w, p_h,
        //         u_cx, u_cy, u_w, u_h
        //     );
        // }
    }

    #[test]
    fn test_mahalanobis_distance() {
        let dt = 0.04;
        let mut kalman = KalmanBBox::new_with_state(
            dt, 1.0, 1.0, 0.0, 0.0,
            2.0, 0.1, 0.1, 0.1, 0.1,
            100.0, 50.0, 40.0, 80.0
        );

        // After prediction, covariance increases
        kalman.predict();

        // Measurement close to predicted state should have small Mahalanobis distance
        let d_close = kalman.mahalanobis_distance_squared(100.5, 50.5, 40.0, 80.0).unwrap();

        // Measurement far from predicted state should have large Mahalanobis distance
        let d_far = kalman.mahalanobis_distance_squared(200.0, 150.0, 60.0, 120.0).unwrap();

        assert!(d_close < d_far, "Close measurement should have smaller distance");

        // Check gating
        let threshold = 9.488; // 95% confidence for 4 DOF
        assert!(kalman.gating_check(100.5, 50.5, 40.0, 80.0, threshold), "Close measurement should pass gate");
    }

    #[test]
    fn test_velocity_tracking() {
        let dt = 0.04;
        let mut kalman = KalmanBBox::new_with_state(
            dt, 0.0, 0.0, 0.0, 0.0, // No control input
            2.0, 0.1, 0.1, 0.1, 0.1,
            100.0, 50.0, 40.0, 80.0
        );

        // Simulate object moving right and down with growing size
        for i in 1..20 {
            let cx = 100.0 + i as f32 * 2.0; // Moving right
            let cy = 50.0 + i as f32 * 1.5;  // Moving down
            let w = 40.0 + i as f32 * 0.5;   // Growing width
            let h = 80.0 + i as f32 * 0.3;   // Growing height

            kalman.predict();
            kalman.update(cx, cy, w, h).unwrap();
        }

        // After tracking, velocity should reflect the motion
        let (vx, vy, vw, vh) = kalman.get_velocity();

        // Velocity should be positive (moving right/down, growing)
        assert!(vx > 0.0, "Velocity X should be positive");
        assert!(vy > 0.0, "Velocity Y should be positive");
        assert!(vw > 0.0, "Velocity W should be positive");
        assert!(vh > 0.0, "Velocity H should be positive");
    }
}

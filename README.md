# Implementation of Discrete Kalman filter for object tracking purposes

The Kalman filter estimates the state of a system at time $k$ via the linear stochastic difference equation considering the state of a system at time $k$ is evolved from the previous state at time $k-1$. See the ref. https://en.wikipedia.org/wiki/Kalman_filter

In other words, the purpose of Kalman filter is to predict the next state via using prior knowledge of the current state. 

In this repository Hybrid Kalman filter is implemented considering continuous-time model while discrete-time measurements. See the ref. - https://en.wikipedia.org/wiki/Kalman_filter#Hybrid_Kalman_filter

## Main algorithm and equations

Define mentioned _linear stochastic difference equation_:

$$x_{k} = A⋅x_{k-1} + B⋅u_{k-1} + w_{k-1} \tag{1}$$

Define measurement model:
$$z_{k} = H⋅x_{k} + v_{k}\tag{2}$$

Let's denote variables:

* $A$ (sometimes it's written as $F$, but I prefer to stick with $A$) - [Transition matrix](https://en.wikipedia.org/wiki/State-transition_matrix) of size $n \times n$ relating state $k-1$ to state $k$
* $B$ - Control input matrix of size $n \times l$ which is applied to *optional* control input $u_{k-1}$
* $H$ - Transformation matrix of size $m \times n$
* $u_{k}$ - Control input
* $w_{k}$ - Process noise vector with covariance $Q$. Gaussian noise with the normal probability distribution:
$$ w(t) \sim N(0, Q) \tag{3} $$
* $v_{k}$ - Measurement noise vector (uncertainty) with covariance $R$. Gaussian noise with the normal probability distribution:
$$ v(t) \sim N(0, R) \tag{4} $$

### Prediction

Let's use the dash sign "$-$" as superscript to indicate the a priory state.

A priory state in matrix notation is defined as
$$\hat{x}^-_{k} = A⋅\hat{x}_{k-1} + B⋅u_{k-1} \tag{5}$$

$$\text{, where $\hat{x}^-_{k}$ - a priory state (a.k.a. predicted),  $\hat{x}_{k-1}$ - a posteriory state (a.k.a. previous)} $$

__Note: A posteriory state $\hat{x}_{k-1}$ on 0-th time step (initial) should be *guessed*__

Error covariance matrix is defined as
$$P^-_{k} =  A⋅P_{k-1}⋅A^{T} + Q \tag{6}$$
$$\text{, where $P_{k-1}$ - previously estimated error covariance matrix, Q - process noise covariance}$$

__Note: $P_{k-1}$ on 0-th time step (initial) should be *guessed*__

### Correction

The Kalman gain (which minimizes the estimate variance) in matrix notation is defined as:
$$K_{k} = P^-_{k}⋅H^{T}⋅(H⋅P^-_{k}⋅H^{T}+R)^{-1} \tag{7}$$

$$\text{, where H - observation matrix, R - measurement noise covariance}$$

After evaluating the Kalman gain we need to update a priory state $\hat{x}^-_{k}$. In order to do that we need to calculate measurement residual:
$$r_{k} = z_{k} - H⋅\hat{x}^-_{k} \tag{8}$$

$$\text{, where $z_{k}$ - true measurement, $H⋅\hat{x}^-_{k}$ - previously estimated measurement}$$

Then we can update predicted state $\hat{x}_{k}$:
$$\hat{x}_{k} = \hat{x}^-_{k} + K_{k}⋅r_{k}$$
$$\text{or} \tag{9}$$
$$\hat{x}_{k} = \hat{x}^-_{k} + K_{k}⋅(z_{k} - H⋅\hat{x}^-_{k})$$

After that we should update error covariance matrix P_{k} which will be used in next time stap (an so on):
$$P_{k} = (I - K_{k}⋅H)⋅P^-_{k}\tag{10}$$
$$\text{, where $I$ - identity matrix (square matrix with ones on the main diagonal and zeros elsewhere)}$$


### Overall
The whole algorithm can be described as high-level diagram:
<p align="center">
<img src="diagram.png" width="720" >
<p align="center">Fig 1. Operation of the Kalman filter. Welch & Bishop, 'An Introduction to the Kalman Filter'</p>
</p>

## 1-D Kalman filter
@todo: physical model / text / code / plots

## 2-D Kalman filter
@todo: physical model / text / code / plots

# Refrences
* [Greg Welch and Gary Bishop, ‘An Introduction to the Kalman Filter’, July 24, 2006](https://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf)
* [Introducion to the Kalman Filter by Alex Becker](https://www.kalmanfilter.net/default.aspx)
* [Kalman filter on wikipedia](https://en.wikipedia.org/wiki/Kalman_filter)
* [State-transition matrix](https://en.wikipedia.org/wiki/State-transition_matrix)

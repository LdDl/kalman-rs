# Implementation of Hybrid Kalman filter for object tracking purposes

The Kalman filter estimates the state of a system at time $k$ via the linear stochastic difference equation considering the state of a system at time $k$ is evolved from the previous state at time $k-1$. See the ref. https://en.wikipedia.org/wiki/Kalman_filter

In other words, the purpose of Kalman filter is to predict the next state via using prior knowledge of the current state. 

In this repository Hybrid Kalman filter is implemented considering continuous-time model while discrete-time measurements. See the ref. - https://en.wikipedia.org/wiki/Kalman_filter#Hybrid_Kalman_filter
## Main algorithm

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
* $v_{k}$ - Measurement noise vector with covariance $R$. Gaussian noise with the normal probability distribution:
$$ v(t) \sim N(0, R) \tag{4} $$

@todo: prediction and update stages
## 1-D Kalman filter
@todo: physical model / text / code / plots

## 2-D Kalman filter
@todo: physical model / text / code / plots
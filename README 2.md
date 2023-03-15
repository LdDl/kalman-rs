Considering acceleration motion again let's write down its equations:

Considering the same physical model as in $(13)$-$(14)$ let's write down state vector $\chi_{k}$:

$$\chi_{k} = \begin{bmatrix}
x_{k} \\
y_{k} \\
x'_{k} \\
y'_{k} \end{bmatrix} = \begin{bmatrix}
x_{k-1} + x'_{k-1}\Delta t + \frac{x''_{k-1}(\Delta t^2)}{2} \\
y_{k-1} + y'_{k-1}\Delta t + \frac{y''_{k-1}(\Delta t^2)}{2} \\
x'_{k-1} + x''_{k-1}\Delta t \\
y'_{k-1} + y''_{k-1}\Delta t
\end{bmatrix} \tag{26}$$

Matrix form of $\chi_{k}$:

$$\chi_{k} = \begin{bmatrix} x_{k} \\
y_{k} \\
x'_{k} \\
y'_{k}
\end{bmatrix} = \begin{bmatrix} 1 & 0 & \Delta t & 0 \\
0 & 1 & 0 & \Delta t \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \end{bmatrix} ⋅ \begin{bmatrix} x_{k-1} \\
y_{k-1} \\
x'_{k-1} \\
y'_{k-1} \end{bmatrix} + \begin{bmatrix} \frac{\Delta t^2}{2} & 0 \\
0 & \frac{\Delta t^2}{2} \\
\Delta t & 0 \\
0 & \Delta t \end{bmatrix} ⋅ \begin{bmatrix} x''_{k-1} \\
y''_{k-1} \end{bmatrix} = \begin{bmatrix} 1 & 0 & \Delta t & 0 \\
0 & 1 & 0 & \Delta t \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \end{bmatrix} ⋅ \chi_{k-1} + \begin{bmatrix} \frac{\Delta t^2}{2} & 0 \\
0 & \frac{\Delta t^2}{2} \\
\Delta t & 0 \\
0 & \Delta t \end{bmatrix} ⋅ \begin{bmatrix} x''_{k-1} \\
y''_{k-1} \end{bmatrix} \tag{27}$$

$$ \text{Assuming that $x''$ and $y''$ - is acceleration $a$, }$$

$$ a_{k-1} = \begin{bmatrix} x''_{k-1} \\
y''_{k-1} \end{bmatrix} \tag{28}$$

$$\chi_{k} = \begin{bmatrix} x_{k} \\
y_{k} \\
x'_{k} \\
y'_{k}
\end{bmatrix} = \begin{bmatrix} 1 & 0 & \Delta t & 0 \\
0 & 1 & 0 & \Delta t \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \end{bmatrix} ⋅ \chi_{k-1} + \begin{bmatrix} \frac{\Delta t^2}{2} & 0 \\
0 & \frac{\Delta t^2}{2} \\
\Delta t & 0 \\
0 & \Delta t \end{bmatrix} ⋅ a_{k-1} \tag{29}$$


Taking close look on $(16)$ and $(1)$ we can write transition matrix $A$ and control input matrix $B$ as follows:

$$A = \begin{bmatrix} 1 & 0 & \Delta t & 0 \\
0 & 1 & 0 & \Delta t \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \end{bmatrix} \tag{30}$$

$$B = \begin{bmatrix} \frac{\Delta t^2}{2} & 0 \\
0 & \frac{\Delta t^2}{2} \\
\Delta t & 0 \\
0 & \Delta t \end{bmatrix} \tag{31}$$

Let's find transformation matrix $H$. According to $(2)$:

$$z_{k} = H⋅\chi_{k} + v_{k} = \begin{bmatrix} 1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \end{bmatrix} ⋅\begin{bmatrix} x_{k} \\
y_{k} \\
{x'_{k}} \\
{y'_{k}} \end{bmatrix} + v_{k} \tag{32}$$

$$ H = \begin{bmatrix} 1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \end{bmatrix} \tag{33}$$

Process noise covariance matrix $Q$:

$$Q = \begin{matrix}
 & \begin{matrix}x && y && x' && y'\end{matrix} \\
\begin{matrix}x \\
y \\
x' \\
y'\end{matrix} & 
  \begin{bmatrix} \sigma^2_{x} & 0 & \sigma_{x} \sigma_{x'} & 0 \\
0 & \sigma^2_{y} & 0 & \sigma_{y} \sigma_{y'} \\
\sigma_{x'} \sigma_{x} & 0 & \sigma^2_{x'} & 0 \\
0 & \sigma_{y'} \sigma_{y} & 0 & \sigma^2_{y'}\end{bmatrix}
 \\\\
\end{matrix} \tag{34}$$

$$\text{, where} $$

$$ \text{$\sigma_{x}$ - standart deviation of position for $x$ component} $$

$$ \text{$\sigma_{y}$ - standart deviation of position for $y$ component} $$

$$ \text{$\sigma_{x'}$ - standart deviation of velocity for $x$ component} $$

$$ \text{$\sigma_{y'}$ - standart deviation of velocity for $y$ component} $$

Since we know about $(14)$ we can define $\sigma_{x}$, $\sigma_{y}$, $\sigma_{x'}$ and $\sigma_{y'}$ as:

$$ \sigma_{x} = \sigma_{x''} \frac{\Delta t^2}{2} \tag{35}$$

$$ \sigma_{y} = \sigma_{y''} \frac{\Delta t^2}{2} \tag{36}$$

$$ \sigma_{x'} = \sigma_{x''} \Delta t \tag{37}$$

$$ \sigma_{y'} = \sigma_{y''} \Delta t \tag{38}$$

$$\text{, where $\sigma_{x''}$ and $\sigma_{y''}$ - standart deviation of acceleration (tuned values)} $$

And now process noise covariance matrix $Q$ could be defined as:

$$ Q = \begin{bmatrix} (\sigma_{x''} \frac{\Delta t^2}{2})^2 & 0 & \sigma_{x''} \frac{\Delta t^2}{2} \sigma_{x''} \Delta t & 0 \\
0 & (\sigma_{y''} \frac{\Delta t^2}{2})^2 & 0 & \sigma_{y''} \frac{\Delta t^2}{2} \sigma_{y''} \Delta t \\
\sigma_{x''} \frac{\Delta t^2}{2} \sigma_{x''} \Delta t & 0 & (\sigma_{x''} \Delta t)^2 & 0 \\
0 & \sigma_{y''} \frac{\Delta t^2}{2} \sigma_{y''} \Delta t & 0 & (\sigma_{y''} \Delta t)^2 \end{bmatrix} = $$

$$ = \begin{bmatrix} (\sigma_{x''} \frac{\Delta t^2}{2})^2 & 0 & (\sigma_{x''})^2 \frac{\Delta t^2}{2} \Delta t & 0 \\
0 & (\sigma_{y''} \frac{\Delta t^2}{2})^2 & 0 & (\sigma_{y''})^2 \frac{\Delta t^2}{2} \Delta t \\
(\sigma_{x''})^2 \frac{\Delta t^2}{2} \Delta t & 0 & (\sigma_{x''} \Delta t)^2 & 0 \\
0 & (\sigma_{y''})^2 \frac{\Delta t^2}{2}\Delta t & 0 & (\sigma_{y''} \Delta t)^2 \end{bmatrix} = \text{| Knowing that $x''$ and $y''$ - acceleration|} = $$ 
$$ = \begin{bmatrix} (\frac{\Delta t^2}{2})^2 & 0 & \frac{\Delta t^2}{2} \Delta t & 0 \\
0 & (\frac{\Delta t^2}{2})^2 & 0 & \frac{\Delta t^2}{2} \Delta t \\
\frac{\Delta t^2}{2} \Delta t & 0 & \Delta t^2 & 0 \\
0 & \Delta t \frac{\Delta t^2}{2} & 0 & \Delta t^2 \end{bmatrix} \sigma^2_{a''}$$

$$ = \begin{bmatrix} \frac{\Delta t^4}{4} & 0 & \frac{\Delta t^3}{2} & 0 \\
0 & \frac{\Delta t^4}{4} & 0 & \frac{\Delta t^3}{2} \\
\frac{\Delta t^3}{2} & 0 & \Delta t^2 & 0 \\
0 & \frac{\Delta t^3}{2} & 0 & \Delta t^2 \end{bmatrix} \sigma^2_{a''} \tag{39}$$

Covariance of measurement noise $R$ is matrix of size $2 \times 2$ (since there are two components - $x$ and $y$) and it is defined as variance of the measurement noise:

$$R = \begin{matrix}
\begin{matrix}& x & y\end{matrix} \\
\begin{matrix}x \\
y \end{matrix}
  \begin{bmatrix}\sigma^2_{x} & 0 \\
  0 & \sigma^2_{y} \end{bmatrix}
 \\\\
\end{matrix} = \begin{bmatrix}\sigma^2_{x} & 0 \\
  0 & \sigma^2_{y} \end{bmatrix} \tag{40}$$

@todo: rust code / plots
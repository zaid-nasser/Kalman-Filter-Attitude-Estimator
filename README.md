# Kalman-Filter-Attitude-Estimator
An Implementation of two algorithms for Estimating the Attitude of a rigid body from IMU sensor data depending on Kalman filter.

## Description:
Each method folder contains the following files:
- method_function.m: contains the main function of the algorithm. It takes as parameters the algorithm's inputs and performs one iteration producing the attitude according to the representation of the algorithm.
- estimate_code.m: Use the method_function.m function to estimate attitude on IMU sensor data and visualize the result. It can be adjusted to load the data you need for attitude estimation.
- a2A.m: function takes a 3d vector as input and returns the skew-symmetric matrix of this vector.
- Omega_Big.m: function takes as an input a 3d vector $`x`$ and returns the matrix:


$$
  \begin{pmatrix}
  0 & -x_3 & x_2 & x_1 \\
  x_3 & 0 & -x_1 & x_2 \\
  -x_2 & x_1 & 0 & x_3 \\
  -x_1 & -x_2 & -x_3 & 0
  \end{pmatrix}
$$


- TETA_BIG.m: function takes as an input a quaternion $`q`$ and returns the matrix:

$$
  \left(\begin{matrix}
  q_4 & -q_3 & q_2  \\
  q_3 & q_4 & -q_1  \\
  -q_2 & q_1 & q_4  \\
  -q_1 & -q_2 & -q_3 
  \end{matrix}\right)
$$

## the first method
This method uses unit quaternions as a representation of the attitude of the rigid body. 

The state vector $\(\mathbf{x}\)$ is defined as the vector (imaginary) part of the error quaternion $\( q_{\text{error}} \)$. The error quaternion $\( q_{\text{error}} \)$ is obtained by multiplying the real quaternion $\( q_{\text{real}} \)$ by the conjugate of the estimated quaternion $\( q_{\text{estimated}} \)$:

$$q_{\text{error}} = q_{\text{real}} \cdot q_{\text{estimated}}^{-1}$$

Thus, the state vector $\( \mathbf{x} \)$ is given by:

$$\mathbf{x} = \text{Im}(q_{\text{error}}) = \text{Im}(q_{\text{real}} \cdot q_{\text{estimated}}^{-1})$$

The state equation of the system is defined as follows:

$$x_k = F_k \cdot x_{k-1} + G_{k} \cdot w_k$$

Where:

$$
\begin{equation}
F_k = \exp\left(\Delta t \cdot \begin{pmatrix}
0 & -\omega_3 & \omega_2 \\\
\omega_3 & 0 & -\omega_1 \\\
-\omega_2 & \omega_1 & 0
\end{pmatrix}\right)
\end{equation}
$$


$$ G_k = \frac{-\Delta t}{2} \cdot I_{3\times 3}
$$

$` \Delta t`$ is the time difference between each two consecutive readings and $` \omega_k =(\omega_1,\omega_2,\omega_3)`$ is the gyroscope reading at time step number $k$. $` w_k`$ is the state model noise with Gaussian distribution $`\sim \mathcal{N}(0,\,Q_k)\,`$.

The measurement equation of the system is defined as follows:

$$
z_k = H_k \cdot x_k + v_k
$$

Where:

$$
H_k = \begin{bmatrix}
2\cdot \begin{bmatrix} C_n^{\text{ }b-} \cdot f_k^{\text{ }n} \times \end{bmatrix} \\
2\cdot \begin{bmatrix} C_n^{\text{ }b-} \cdot M_k^{\text{ }n} \times \end{bmatrix}
\end{bmatrix}
$$

$`C_n^{\text{ }b-}`$ is the Rotation matrix from the NED frame to the body frame calculated from our gyroscope reading. $`f_k^{\text{ }n}`$ the gravity vector in the NED frame. $`M_k^{\text{ }n}`$ the magnetic field vector in the NED frame. Finally, the notation $` \begin{bmatrix} x \times \end{bmatrix}`$ is the skew-symmetric matrix of the vector $`x`$. $`v_k`$ is the measurments noise with Gaussian distribution $`\sim \mathcal{N}(0,\,R_k)\,`$.

### Algorithm steps:

- Estimating the quaternion -prior estimation- using Gyroscope reading (Prediction stage):

  $`
 \bar{q}_k^{-} = \exp(\frac{\Delta t}{2} \cdot \Omega(\omega_{k-1})) \cdot \bar{q}_{k-1}^+
 `$


 

- Estimating covariance matrix of the state vector:

  $` P_k^{-} = F_{k}\cdot P_{k-1}^+ \cdot F_{k}^T + G_k \cdot Q_k \cdot G_k^T `$

- Calculating Kalman Gain:

  $` K_k = P_k^-\cdot H_k^T(H_k \cdot P_k^-\cdot H_k^T + R_k )^{-1} `$

- Updating the State vector after the measurement:

  $` x_k^+ =K_k(z_k - \hat{z}_k) =K_k(z_k - \begin{bmatrix}
  C_n^{\text{ }b-} \cdot f_k^{\text{ }n} \\
    C_n^{\text{ }b-} \cdot M_k^{\text{ }n}
\end{bmatrix})`$

- Updating the covariance matrix of the state vector:

  $` P_k^+ = (I - K_k\cdot H_k)P_k^- (I - K_k\cdot H_k)^T + K_k\cdot R_k \cdot K_k^T `$

- Updating the quaternion -posterior estimation-:
  
  $`\bar{q}_k^{+} = \bar{q}_k^{-} + \begin{bmatrix} rel(\bar{q}_k^{-}) \cdot I_{3\times 3} \\ -Im(\bar{q}_k^{-})^T\end{bmatrix} \cdot x_k^+`$

## The second method
This method aims to separate the external acceleration applied to the rigid body from the accelerometer sensor reading. Consequently, the only component left in the accelerometer reading is the gravity field vector. The prediction stage of this method uses the gyroscope readings as the system model like the first method. But in this method, we use Euler angles instead of unit quaternions to represent the attitude and we only care about the roll and pitch angles. The yaw angle will not be estimated since the magnetometer sensor is not involved in this method.

The state vector is $`x_k = [g_k^b a_k^b]`$. Where $`g_k^b`$ is the gravity field vector at time step $`k`$ expressed in the body frame and $`a_k^b`$ is the external acceleration applied to the rigid body at time step $`k`$ expressed in the body frame. 

The state equation of the system is as follows:

$$
x_k = F_k \cdot x_{k-1} + G_k \cdot w_k
$$

Where:

$$
F_k = \begin{pmatrix} 
\exp (-\Delta t [\omega_k \times]) & 0_{3\times 3} \\
0_{3\times 3} & c_a \cdot I_{3\times 3} 
\end{pmatrix}
$$

$$
G_k = \begin{pmatrix} 
 -\Delta t [g_{k-1}^b \times]) & 0_{3\times 3} \\
0_{3\times 3} &  I_{3\times 3} 
\end{pmatrix}
$$

$` \Delta t`$ is the time difference between each two consecutive readings and $` \omega_k =(\omega_1,\omega_2,\omega_3)`$ is the gyroscope reading at time step number $k$. $` w_k`$ is the state model noise with Gaussian distribution $`\sim \mathcal{N}(0,\,Q_k)\,`$. $`c_a`$ is a coefficient between 0 and 1 used to model the external acceleration of the rigid body depending on the Gauss-Markov model.

The measurement equation of the system is defined as follows:

$$
z_k = H_k \cdot x_k + v_k
$$

Where:

$$
H_k = [-I_{3 \times 3} -I_{3 \times 3}]
$$

$`v_k`$ is the measurments noise with Gaussian distribution $`\sim \mathcal{N}(0,\,R_k)\,`$.

The Steps of this algorithm is straight forward according to the Kalman filter steps (prediction and correction). Thus I will only illustrate how to obtain the Euler angles from the state vector.

Roll angle:


$$
\Phi = \text{atan2}(x{_2},x{_3})
$$

Pitch angle:


$$
\Theta = -\text{asin}\frac{-x_1}{g}
$$

For a detailed explanation of the method look at the reference:
Fourati, Hassen, and Djamel Eddine Chouaib Belkhiat. Multisensor attitude estimation: fundamental concepts and applications. CRC Press, 2016.

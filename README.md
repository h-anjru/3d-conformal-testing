# 3d-conformal-testing
A series of functions purpose built for testing and comparing various methods of performing the 3D conformal coordinate transformation (i.e. absolute orientation).

## Overview
The 3D conformal coordinate transformation has applications in photogrammetry, lidar, robotics, and geodesy. In these applications, one must express the coordinates of a set of 3D points known in one coordinate system relative to another coordiante system. For example, data collected by a lidar scanner are initially expressed in a coordinate system which is relative to the scanner, but the desired result is to express these data in terms of a world coordinate system. One means of finding the relationship between two coordinate systems is to measure the coordinates of a set of points common to both the coordinate systems. Once the relationship between the two coordinate systems—the transformation—is solved, the whole of the point data in one system can be transformed to the other.

The 3D conformal coordinate transformation is commonly expressed as seven parameters: A uniform scale factor, three sequential rotation angles (one about each of the three axes of the target coordinate system), and three translations (one along each of the three axes of the target coordinate system). This transformation preserves the shape of the data being transformed, i.e., no skew, non-uniform scaling, perspective warp, etc. is introduced.

There are a number of means of solving for the transformation, each with their advantages and shortcomings. The purpose of this experiment is to explore these pros and cons among a number of solution methods for the purpose of providing guidance to researchers and practitioners who wish to examine the efficacy of any software packages which perform a 3D conformal transformation, or to develop their own implementations of the transformation.

## Example usage
After a dataset is generated, the transformation can be solved using the available methods:
```
% Generate the dataset
[arb, con, hgt, noise] = generate3DPoints(5, 0.1);

% Solve for transformation using Horn closed-form method
hgt_horn = hornConf3D(arb + noise, con);

% Solve via Lassiter method (direct linear estimation of rotation only)
[hgt_las, jac_las, K_las] = lasConf3D(arb + noise, con);

% Solve via direct linear estimation of all parameters
[hgt_las2, jac_las2, K_las2] = lasConf3D_2(arb + noise, con);
```
This basic workflow can be nested inside a `for` loop, and outputs stored approriately, for the purposes of Monte Carlo simulations or other tests.

## Horn closed-form unit quaternion method
`hornConf3d()` performs BKP Horn's algorithm for a closed-form solution to the absolute orientation problem (i.e., a 3D conformal transformation)[^1]. 

For a set of $n$ points $\mathbf{r}_i = \langle x_i,y_i,z_i \rangle \ \forall i \in (1,...,n)$ whose coordinates are known in both a left and right coordinate system, the algorithm first solves for the best rotation between the systems before solving for scale and translation. The algorithm, as well as some notes about its implementation in MATLAB, are summerized below. (I reserve derivation and detailed explanation of each step in this algorithm to Horn's concise manuscript.)

### 1. Translate to origin
Find the centroids of the set of points in each coordinate system $\bar{\mathbf{r}_l}$ and $\bar{\mathbf{r}_r}$. Subtract each from its respective set. This yields two sets of translated point sets, denoted by a prime: $\mathbf{r}'_{l,i}$ and $\mathbf{r}'_{r,i}$.

### 2. Matrix of sums of products
There is a $3 \times 3$ matrix whose elements are the sums of products of coordinates in the left system multiplied coordinates in the right system. As notated by Horn, this matrix is
```math
M=\begin{bmatrix}
S_{xx} & S_{xy} & S_{xz} \\
S_{yx} & S_{yy} & S_{yz} \\
S_{zx} & S_{zy} & S_{zz}
\end{bmatrix}
```

where

$S_{xx}= \sum_ {i=1} ^{n} = x' _{l,i} x'_ {r,i}$
$S_{xy}= \sum_ {i=1} ^{n} = x' _{l,i} y'_ {r,i}$
etc.

If the points as column vectors are adjoined as $3 \times n$ matrices $M_l$ and $M_r$ then $M=M_l M_r^T$.

### 3. Real symmetric matrix $N$
There is a $4 \times 4$ symmetric matrix $N$ which can be found using the elements of matrix $M$:

```math
N=\begin{bmatrix}
S_{xx} + S_{yy} + S_{zz} & S_{yz} - S_{zy} & S_{zx} - S_{xz} & S_{xy} - S_{yx} \\
S_{yz} - S_{zy} & S_{xx} - S_{yy} - S_{zz} & S_{xy} + S_{yx} & S_{zx} + S_{xz} \\
S_{zx} - S_{xz} & S_{xy} + S_{yx} & -S_{xx} + S_{yy} - S_{zz} & S_{yz} + S_{zy} \\
S_{xy} - S_{yx} & S_{zx} + S_{xz} & S_{yz} + S_{zy} & -S_{xx} - S_{yy} + S_{zz}
\end{bmatrix}
```

### 4. Finding the rotation as a unit quaternion
In the manuscript it is shown that there is a unit quaternion $\dot{q}$ which describes the rotation between the left and right coordinate systems. It is also shown that the unit quaternion $\dot{q}$ which maximizes the product $\dot{q}N\dot{q}$ is the eiginvector which corresponds to the most positive eigenvalue of natrix $N$. This eignvector is found using `[V, D] = eig(N)` and converted to a rotation matrix using `quat2rotm()`.

### 5. Finding scale
Horn presents three formulas for finding the scale of the transformation. In this code I implement the symmetric scale formula:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$s=\sum_{i=1}^{n} ||\mathbf{r}'_{r,i}|| \ \ / \ \sum_{i=1}^{n} ||\mathbf{r}'_{l,i}||$

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Using the methods shown by Horn (and others), scale can be determined without knowing the rotation, and the solution for rotation is not affected by the scale. 

### 6. Finding the translation
The translation is the difference between the centroid of the right coordinate system and the scaled and rotated centroid of the left coordinate system.

### Discussion of the Horn method
There is a clear advtantage in this method in how simple it is to implement in code, its speed when compared to iterative methods, and the stability of being a closed-form solution. To add, though not mentioned in the summary above, a stochastic model can be added to the algorithm by using weighted means to find the centroids, weighting the norms when finding scale, and weighting the sums of products in matrix $M$. These weights can implement estimates of uncertainty from the measurements in both coordinate systems, which is an advantage over the stochastic model implemented in a typical least squares adjustment.

This solution has two disadvantages, however. First, there is no manner to look for outliers in the data sets, and there is no point in the closed-form solution which would indicate that potential outlier points have been included in the solution. Second, evaluating the quality of the solution after completion is difficult. While residuals can be found as the differences between the measured points in one or both of the coordiante systems versus the transformed points, conventional error propagation cannot be performed in the absence of a coefficient or Jacobian matrix. While the residuals can provide an estimate of the uncertainties of the observations, there is, to my knowledge, no means of ascertaining estimates of the uncertrinaties of the transformation parameters.

[^1]: Berthold K. P. Horn, "Closed-form solution of absolute orientation using unit quaternions," J. Opt. Soc. Am. A 4, 629-642 (1987). https://doi.org/10.1364/JOSAA.4.000629
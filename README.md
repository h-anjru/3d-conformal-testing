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
horn_hgt = hornConf3D(arb + noise, con);

% Solve via Lassiter method (direct linear estimation of rotation only)
[las_hgt, las_jac, las_K] = lasConf3D(arb + noise, con);

% Solve via direct linear estimation of all parameters
[las2_hgt, las2_jac, las2_K] = lasConf3D_2(arb + noise, con);
```
This basic workflow can be nested inside a `for` loop, and outputs stored approriately, for the purposes of Monte Carlo simulations or other tests.

##
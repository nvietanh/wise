# WISE: Wasserstein Inverse-covariance Shrinkage Estimator, MATLAB implementation

This folder contains the source code for the implementation, and all the experiments in the paper.

## Installation guide

Run the install.m script to add the folder '/src' to the working directory of MATLAB. One can use the savepath 

command after running the install.m script to save the path permanently.

## Utilisation guide

The package comes with 2 main functions: wise and wise_structure. The test.m file contains a simple example.


1. The wise command

The wise command returns the Wasserstein inverse covariance estimate by calculating the analytical solution. The user can specify the tolerance for the eigenvalues and for the accuracy of the bisection algorithm.

```matlab
est = wise(S\_hat, rho);
est = wise(S\_hat, rho, 1e-8, 1e-5);
```

2. The wise_structure command

The wise_structure command calculates the Wasserstein inverse-covariance estimate with the sparsity structure 
imposed by the matrix $E$ using the successive quadratic approximation method.

```matlab
options = wise\_structure\_settings('iter\_limit', 100, 'sigma', 1e-4, 'delta\_tol', 1e-6, 'gradient\_tol', 1e-4);
est = wise\_structure(S_hat, rho, E, options);
```
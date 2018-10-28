# WISE: Wasserstein Inverse covariance Shrinkage Estimator, MATLAB implementation

This folder contains the source code for the implementation, and all the experiments in the paper.

## Installation guide

Run the install.m script to add the folder '/src' to the working directory of MATLAB. One can use the savepath 

command after running the install.m script to save the path permanently.

## Utilisation guide

The package comes with 2 main functions: wise and wise_structure. The test.m file contains a simple example.


### The wise command

The wise command returns the Wasserstein inverse covariance estimate by calculating the analytical solution. The user can specify the tolerance for the eigenvalues and for the accuracy of the bisection algorithm.

```matlab
est = wise(S_hat, rho);
est = wise(S_hat, rho, 1e-8, 1e-5);
```

The function returns a matlab structure:

1. est.value: the precision matrix estimate
2. est.x: the eigenvalues of the estimator
3. est.eigbasis: the eigenbasis of the estimator
4. est.gamma: the optimal dual variable
5. est.min_gamma, est.max_gamma: the range of gamma used for bisection
6. est.msg: output message
    

### The wise_structure command

The wise_structure command calculates the Wasserstein inverse-covariance estimate with the sparsity structure 
imposed by the matrix $E$ using the successive quadratic approximation method.

```matlab
options = wise_structure_settings('iter_limit', 100, 'sigma', 1e-4, 'delta_tol', 1e-6, 'gradient_tol', 1e-4);
est = wise_structure(S_hat, rho, E, options);
```


The function return a matlab structure:
1. est.value: the precision matrix estimate with sparsity constraints
2. est.msg: output message
# WISE: Wasserstein Inverse-covariance Shrinkage Estimator, Synthetic experiment

This folder contains the source code for the for synthetic experiment.

Notice that it may take a long time to run the synthetic data due to the fine-scale search grid. The result data can be downloaded from the author's main [website](https://www.vietanhnguyen.net). 

## Run to get the sweet spot curve

Run file synthetic_run.m to calculate the necessary data to plot the sweet spot curve (Figure 2 and 3). We choose a very fine vector for the alpha, beta and rho vector to have a smooth curve.

## Run to get the learning curve

Run file learning_curve_run.m to calculate the data for the learning curve. We choose a very fine vector for rho.
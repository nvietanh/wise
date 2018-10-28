# WISE: Wasserstein Inverse covariance Shrinkage Estimator, Synthetic experiment

This folder contains the source code for the for synthetic experiment.

Notice that it may take a long time to run the synthetic data due to the fine-scale search grid. The result data can be found in the result folder along with the script to plot Figure 2, 3 and 4 in the paper. 

## Run to get the sweet spot curve

Run file synthetic_run.m to calculate the necessary data to plot the sweet spot curve (Figure 2 and 3). We choose a very fine vector for the alpha, beta and rho vector to have a smooth curve.

## Run to get the learning curve

Run file learning_curve_run.m to calculate the data for the learning curve. We choose a very fine vector for rho.

## Solution data file

In the folder result, we have included some .mat file containing the numerical results of the experiment at low resolution. Please be advised that this data is not comprehensive due to the file size restriction of github. In order to have a full smooth curve as in the paper, please contact the authors to have the full numerical data file.
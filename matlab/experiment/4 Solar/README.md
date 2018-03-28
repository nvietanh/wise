# WISE: Wasserstein Inverse covariance Shrinkage Estimator, Solar irradiation experiment

This folder contains the source code for the for minimum variance portfolio optimization experiment.

1. Run file solar_run.m to execute the experiment. 

2. Folder 'data' contains the input data for the experiment. Since the original data for the whole Switzerland is only available on SwissMeteo website, we include only the input data which is used for the experiment. The data contains the following variables:
	- kfold: the data used for running the experiment, including 13 cells corresponding to 13 years of data. 			-- kfold{i}.training is the sample covariance matrix estimated for year i.
		-- kfold{i}.test is the sample covariance matrix estimated using the rest of the data. 
	
	- E: matrix specifies the sparsity structure

	- SAI: the Swiss Average Irradiation over the period of 13 years. One can use the command `imagesc(SAI)` to visualize this matrix.

3. Folder 'result' contains the .mat result files.
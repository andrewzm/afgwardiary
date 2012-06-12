If you have the generated MAT files Results_without_cov_estimation.mat and Results_with_cov_estimation.mat, for result plotting simply run Analyze_results.m. If you do not have the kernel result mat files the program will stop prematurely.

If you only have Test_data.mat follow the following steps use VBIDE.m to carry out the inference with and without precision estimation as follows
VBIDE('Results_without_cov_estimation','n');
VBIDE('Results_with_cov_estimation','y'). This may take a few hours.

If you wish to generate results using a nonparametric intensity function estimator run Kernelestimation.m (this uses some variables from Results_with_cov_estimation so should be run after the previous statements). Warning: This will take a long time and generate Kernel_result1.mat - Kernel_result29.mat for different kernel configurations.




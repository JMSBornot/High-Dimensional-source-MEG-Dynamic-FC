# High-Dimensional-Source-MEG-Dynamic-FC

These scripts and functions are a template. Code is designed to run using a parallel cluster in Matlab (parfor). Analysis are run for all the subjects, individually, by modifying particular codes. In general, some setting variables should be modified to reflect your directory path, accordingly, and other information related to particular analyses.

NOTE: Our tools require that users setup the OHBA software library accordingly to the authors' recommendations (see https://ohba-analysis.github.io/osl-docs/).

Keep tuned! We will make available a friendly GUI in a future research work!

CODE **** CPSpipe01_preprocessing.m - data preprocessing ****

CODE **** run_dataset.m ****

This small script starts the parallel cluster and runs the initial MATLAB parallel (parfor) code in a HPC cluster. Exactly, it runs these two scripts:
- CPSpipe02_blockwise_FCandSTAT.m
- CPSpipe03_blockwise_assemble.m

CODE **** CPSpipe02_blockwise_FCandSTAT.m - blockwise functional connectivity (FC) and statistical (Wilcoxon sign-rank) statistical analysis ****

This function performs bivariate FC analysis using an imaginary coherence index. Instead of directly estimating the 16403 x 16403 connectivity matrix, it performs this computation blockwise. Simultaneously, it performs statistical analysis of the post-stimulus estimated dynamic FC values versus corresponding baseline samples. It also apply thresholding to select outstanding features according to predetermined suprathreshold p-values in the range from 10^-9 to 10^-5.

CODE **** CPSpipe03_blockwise_assemble.m - assemble the blockwise FC indices from the block subspace to the full-space of 16403 x 16403 FC indices  ****

This function only performs the assembling of the highly-sparse selected suprathreshold connections.

CODE **** scADSpain_preproc_parfor.m - data preprocessing ****

1. Read the ELEKTA-Neuromag "fif" file and read it using SPM toobox, which generates an object containing the data.
2. Data is downsample to 200 Hz using an SPM function: "spm_eeg_downsample".
3. Individual mri image is read and meshes are extracted using SPM tools.
4. Anatomical and MEG spaces are coregistered using the semi-automatic tool "spm_coreg", which was programmed by the authors with the purpose of providing a more flexible coregistration.
5. The SPM MEG data object is saved with the fields updated for the coregistered case.

NOTE: Some variables have to be setup accordingly to read your data, for all the subjects

CODE **** spm_coreg.m - semi-automatic MRI-MEG coregistration tool ****

NOTE: It is used only inside "scADSpain_preproc_parfor" as part of the coregistration process.

CODE **** scADSpain_invsol_parfor.m - estimate the source activity ****

1. Reads the SPM12 object containing participants data.
2. Pass-band data filtering using SPM12 function "spm_eeg_filter".
3. Runs the Bayesian minimum norm method as implemented in SPM12 (COH) using the function "spm_eeg_invert", for 2s segmented data with applied Hanning tapering per segment.
4. FFT analysis per segment.
5. FFT coefficients saved to hard-disk for post-hoc FC analysis

CODE **** scADSpain_sourceFC_parfor.m - Run the source FC analysis ****

1. FC analysis per frequency or frequency band (uncomment code as suits your needs)
2. Estimate imaginary coherence (iCOH) and the EIC method, other methods can be inserted in the code as needed.
3. Save the outcome to hard-disk

CODE **** scADSpain_clustperm_indsel_parfor.m - Run the cluster-permutation analysis ****

This script is controlled by several flag variables to control the running of particular cell code. The flag variables are:

1. fcompute_block - to estimate the sparse supra-threshold FC indices for each FC block and frequency band.
2. itype - currently take values from 1 to 7 to compute above analysis for different methods (Wilcoxon, Spearman, standard correlation, etc.) and using different psychological tests (MMSE, IRM, DRM).
3. fcompute_stat - assemble the FC indices for all the blocks, separately for each frequency band, and save the outcome.
4. fvisualize_stat - it is settled to true in the code because we always show some part of the outcome at the end as an inspection of the results. This last part use the function "plot_bsurfconn_spm" to provide a visualization of the cortico-cortical FC map in the SPM12 template surface.

NOTE: Users have to specify the particular thresholds for the upper and lower distribution tails, correspondingly to the rank-sum, or correlation analyses, accordingly with its particular dataset and analysis.

CODE **** plot_bsurfconn_spm.m - Provide visualization of FC map for spm template surface ****

NOTE: It is used only inside "scADSpain_clustperm_indsel_parfor".

CODE **** plot_bsurfconn.m - Provide visualization of FC map for generic surfaces ****

NOTE: It is used only inside "plot_bsurfconn_spm".

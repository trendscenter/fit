# **FIT Updates**

## FITv2.0e (Jan 21, 2020):
- Toolbox is updated to include three-way Parallel group ICA+ICA. Please see "Three-way parallel group independent component analysis: Fusion of spatial and spatiotemporal magnetic resonance imaging data", S. Qi, R. F. Silva et al., HBM, 2022, Volume 43, Issue 4, March 2022:1280-1294.
- Two-way parallel group ICA + ICA is now added. Please see "Parallel group ICA+ICA: Joint estimation of linked functional network variability and structural covariation with application to schizophrenia", S. Qi, J.Chen et al., HBM, 2019, Volume 40, Issue 13, September 2019:3795-3809.
- aNy-way fusion based on Independent Component Analysis is now added. Please see "aNy-way Independent Component Analysis", K. Duan, R. F. Silva, V. D. Calhoun, and J. Liu, Annu Int Conf IEEE Eng Med Biol Soc. 2020 Jul; 2020: 1770–1774.
- Neural net fusion (now Deep Fusion) is added in the parallel ICA toolbox. Please see "Reading the (functional) writing on the (structural) wall: Multimodal fusion of brain structure and function via a deep neural network based translation approach reveals novel impairments in schizophrenia", S. Plis, F. Amin et al., NeuroImage, Volume 181, 1 November 2018, Pages 734-747.

***

## FITv2.0d (Aug 29, 2017):
- Transposed IVA (tIVA) and Multi-set Canonical Correlation Analysis (MCCA) toolboxes are added. Please see paper for more information. T. Adali, Y. Levin-Schwartz, and V. D. Calhoun, "Multi-modal data fusion using source separation: Two effective models based on ICA and IVA and their properties", Proc IEEE Inst Electr Electron Eng. 2015 Sep 1; 103(9): 1478–1493.
- Three way parallel ICA option is added in the parallel ICA toolbox. Please see "A three-way parallel ICA approach to analyze links among genetics, brain structure and brain function", V. Vergara, A. Ulloa, V. D. Calhoun, D. Boutte, J. Chen, J. Liu", NeuroImage, 2014 Sep; 98:386-94.
- PCA-CCA based dimensionality estimation is now added in the fusion toolbox. Please see "Sample-Poor Estimation of Order and Common Signal Subspace with Application to Fusion of Medical Imaging Data", Y. Levin-Schwartz, Y. Song et al., Neuroimage, PMC 2017 Jul 1.
- Reference based parallel ICA is now included. Please see "Multi-Reference Parallel ICA: A Semi-blind Multivariate Approach", J. Chen, V. D. Calhoun, and J. Liu, Conf Proc IEEE Eng Med Biol Soc. 2014; 2014: 6659–6662.

***

## FITv2.0c (May 10, 2012):
- ICASSO plugin is now added in the Fusion ICA toolbox (Joint ICA, CCA + joint ICA and Parallel ICA).
- Added dimensionality estimation for SNPs based on paper J. Chen, V. D. Calhoun, J. Liu, "ICA Order Selection Based on Consistency: Application to Genotype Data" (under review EMBS 2012).
- Added option to do permutation test in parallel ICA utilities menu.
- Option is now provided to do two sample t-test on loading coefficients in the parallel ICA toolbox.

***

## FITv2.0b (Mar 20, 2009):
- Combi ICA algorithm developed by Petr Tichavsky, Zbynek Koldovsky, Arie Yeredor, German Gomez-Herrero and Eran Doron is now added to the FIT toolbox.
- Added an option to run analysis only on pair wise and individual combinations of features. To run this, set STACK_ALL_FEATURES to 0 in ica_fuse_defaults.m and set OPTIMIZE_FEATURES to 'Yes'. Please note that with this option you will be able to use only utilities like optimal features and plotting histograms.
- Analysis code is cleaned to handle the memory issues.
- Added a multi bar plot to view optimal feature results when large number of feature combinations are run.
- Removed the limitation to use MATLAB statistics toolbox when components are sorted based on two sample t-test on mixing coefficients. Sorting components based ontwo sample Kolmogorov-Smirnov test is removed for now.
- Removed the limitation to use at most two features when components are sorted based on spatial divergence.

***

## FITv2.0a (Feb 8, 2008):
- Parallel ICA method implemented by [**Jean Liu**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/fit/publications/2007_Liu_SNPfMRI_final.pdf) is included to examine the shared information across modalities.
- Main GUI is changed to include the parallel ICA method.

***

## FITv1.2a (Jan 15, 2008):
- Sorting GUI now contains an option to use feature histograms when spatial divergence sorting criteria is selected.
- Option is provided in display GUI to create EEG-fMRI movie when fusion is done using EEG and fMRI modalities.
- HTML help manual is provided to the run the analysis quickly.
- Latest SPM updates for volume functions are installed.
- Flip parameter for analyze images compatible with SPM2 is stored infusion information file. If the flip parameter is changed during display a warning message is printed.
- Convention for plotting right and left text is changed to make it the same as in SPM.

***

## FITv1.1a (April 25, 2007):
- Option is provided to specify mask for each feature.
- Now back reconstructs group components.
- Magnified view of each joint component is provided when you click image orERP signal.
- EEG data is down sampled while saving the data to the disk.
- Now reads Nifti-data.
- Now allows fMRI features to have different dimensions.
- Batch script is also updated.

***

## FITv1.0_beta: (Aug 2, 2006):
- New GUI is provided which implements joint ICA method. GUI contains user interface controls to do pre-analysis, analysis and display steps.
- Option is provided to sort joint ICA components.
- Batch script is provided to complement the user interface.

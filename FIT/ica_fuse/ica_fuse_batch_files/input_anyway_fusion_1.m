%% INPUT FILE FOR RUNNING Anyway ICA FUSION %%
% After entering the parameters, call batch script using ica_fuse_anyway_fusion_batch(inputFile); at the MATLAB command prompt.
% The input parameters are defined below. Please use full path for files or directory

%% Output directory
outputDir = 'C:\Users\srrac\gsu\srinivas\results\test_fusion\results';

%% Output files prefix
% Specify a value that generates a valid variable name.
% Output prefix must not contain Empty space, \, /, ?, :, *, ", <, > or any
% other character that doesn't generate a valid variable name.
prefix = 'Aod';

%% Mask file
% 1. Default mask - Specify empty value like [].
%
% 2. When specifying maskFiles enter in a cell array. For each feature there
% is a mask.
% fMRI - full file path
% Gene - indices like '1:367'
% sMRI - full file path
%
% NOTE:
% To replicate the same mask over features you can use repmat for example
% maskFile = repmat({'I:\Fusion_Data\fmri_fmri\SZ\s271(201)\targets_ME.img'}, 2, 1);
maskFile = [];

%% Feature names 
% Enter feature names in a cell array of length 2
featureNames = {'fMRI', 'GENE'};

%% Modality selection.
% This will be a vector of length equal to 2.
% Please see variable PARALLEL_ICA_MODALITIES in ica_fuse_defaults.m for available
% modalities.
modality = {'fmri', 'gene'};

%% Input files is a cell array of length number of features. Each cell
% contains a character array of file names.
load 'C:\Users\srrac\gsu\srinivas\data\Fusion_Data\fmri_erp' erp_files fmri_files
files = {fmri_files, erp_files};

%% Number of components. Length of the vector is equal to the number of features.
numComp = [5, 5];

%% Scale components to z-scores. Options are 'none' and 'z-scores'
scale_components = 'z-scores';

%% Algorithm options
algorithm_opts.lambda_iva = 0.2;


%% display
display_results.format = 'html';
display_results.image_values = 'positive and negative';
display_results.anatomical_plane = 'axial';
display_results.slices_in_mm = -40:4:72;
display_results.convert_to_z = 'yes';
display_results.threshold = 1.5;
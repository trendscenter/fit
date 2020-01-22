%% INPUT FILE FOR RUNNING PARALLEL ICA FUSION %%
% After entering the parameters, call batch script using ica_fuse_parallelICA_batch_file(inputFile); at the MATLAB command prompt.
% The input parameters are defined below. Please use full path for files or directory

%% Output directory
outputDir = '/na/homes/test1';

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

%% Group names
% Enter group names in a cell array
groupNames = {'Healthy', 'Schizophrenia'};

%% Feature names 
% Enter feature names in a cell array of length 2
featureNames = {'fMRI', 'GENE'};

%% Modality selection.
% This will be a vector of length equal to 2.
% Please see variable PARALLEL_ICA_MODALITIES in ica_fuse_defaults.m for available
% modalities.
modality = {'fmri', 'gene'};


%% Data selection. Options are 1 and 2
% 1 - Specify directory and file patterns
% 2 - Paste filenames
dataSelectionOption = 2;

%% Data Selection (Option 1)

% Answer for file pattern.
% 1 means the file pattern will be same between groups.
% 0 means the file pattern will be different between groups
answerFilePattern = 1;

% Answer for directory
% 1 means all the data is organized in one folder
% 0 means data is organized in different folders
answerDir = 0;

% Groups file pattern
% If answerFilePattern is 1 then only group1_file_pattern variable will be
% read.
group1_file_pattern = {'con_0008.img', '*.asc'}; % group 1 file pattern
group2_file_pattern = {'con_0008.img', '*.asc'}; % group 2 file pattern

% Specify the root folder for groups. The data will be selected as follows:
%
% 1. Data with the specified file pattern will be searched first in the root folder and each file will be treated
% as an individual for that group.
% 2. If the files are not in root folder and each subject for that group
% has a sub-folder, then individuals for that group will be selected based on the matching file pattern. First it will
% search subfolder followed by sub-subfolders. If the data is missing error
% will be printed to the MATLAB command prompt.
%
% NOTE:
% 1. After data is selected, the files selected for each group's task
% will be printed to a text file. You can review this text file.
% 2. Keep one file per subject for each group's task so that the toolbox selects the data correctly without reviewing the text file.

group1_dir = '/na/homes/test1/fusion_data/fmri_gene/Healthy';
group2_dir = '/na/homes/test1/fusion_data/fmri_gene/SZ';

%% Data Selection (option 2)
% Files is a cell array of dimensions no. of groups by features. Each cell
% contains a character array of file names.
files = {'F:\Fusion_Data\fmri_fmri\Healthy\tar*.img', 'F:\Fusion_Data\fmri_fmri\Healthy\*.asc';
    'F:\Fusion_Data\fmri_fmri\SZ\tar*.img', 'F:\Fusion_Data\fmri_fmri\SZ\*.asc'};


%% Number of components. Length of the vector is equal to the number of features.
numComp = [5, 7];

%% Type of ICA analysis:
% Options are 'average' and 'icasso'
type_ica = 'icasso';

%% Number of times you want to run ICA
num_ica_runs = 15;

%% Type of parallel ICA
% Options are AA, AS and AA-ref
% 1. AA - Correlation between modalities is computed using loading
% coefficients.
% 2. AS - Correlation between modalities is computed using loading
% coefficients of first modality with the source of the second modality.
% 3. AA-ref - Parallel ICA with SNP reference file as constraint. Specify
% SNP reference file under snp_ref_file variable.
type_parallel_ica = 'AA';

%% Type of PCA.
% Options are Reference and Standard
% Note that if you have selected AA as the type for parallel ICA you will
% have the option to select reference or standard for type of PCA. If you
% have selected AS  as the type for parallel ICA, standard option will be
% automatically selected.
type_pca = 'reference';

%% Reference vector to distinguish groups if you have selected AA as the type of parallel ICA.
% In the example data there are 63 subjects, 1 and -1 are used as flags for
% Healthy and patients. Leave it as empty if you are not sure about this.
% NOTE: The length of this vector must equal the total number of subjects.
reference = [ones(1, 43), -ones(1, 20)];

%% SNP reference file name for using AA-ref algorithm
snp_ref_file = '';


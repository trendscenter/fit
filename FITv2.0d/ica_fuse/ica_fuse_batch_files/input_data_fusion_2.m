%%%%%%% INPUT FILE FOR RUNNING FUSION USING A BATCH SCRIPT %%%%%%%
% The input parameters are defined below. Please use full path for files or
% directory

%% Output directory
outputDir = 'E:\test_fusion\batch_results\';

%% Prefix
% Output files prefix. Specify a value that generates a valid variable name.
% Output prefix must not contain Empty space, \, /, ?, :, *, ", <, > or any
% other character that doesn't generate a valid variable name.
prefix = 'Aod';

%% Mask file
% 1. Default mask - Specify empty value like [].
%
% 2. When specifying maskFiles enter in a cell array. For each feature there
% is a mask. 
% fMRI - full file path
% sMRI - full file path
% EEG - indices like '50:500'
% Note:
% To replicate the same mask over features use repmat for example
% maskFile = repmat({'I:\Fusion_Data\fmri_fmri\SZ\s271(201)\targets_ME.img'}, 2, 1);
maskFile = [];

%% Normalization. 1 is default
normalize = 1;

%% Group names
groupNames = {'Healthy', 'Schizophrenia'};

%% Feature names
featureNames = {'Auditory Oddball', 'Sternberg'};

%% Modalities
% Modality selection. This will be a vector of length equal to length of featureNames.
% See variable AVAILABLE_MODALITIES in ica_fuse_defaults for available
% modalities.
modality = {'fmri', 'fmri'};


%% Data selection. Options are 1 and 2
% 1 - Specify directory and file patterns
% 2 - Paste filenames
dataSelectionOption = 2;

%% Data Selection (Option 1)

% answer for file pattern. 
% 1 means the file pattern will be same between groups.
% 0 means the file pattern will be different between groups
answerFilePattern = 1;

% answer for directory
% 1 means all the data is organized in one folder
% 0 means data is organized in different folders
answerDir = 0;

% group1 file pattern
% if answerFilePattern is 1 then only group1_file_pattern variable will be
% read.
group1_file_pattern = {'tar*.img', 're*.img'}; % group 1 file pattern
group2_file_pattern = {'tar*.img', 're*.img'}; % group 2 file pattern

% Specify the root folder for groups. The data will be selected as follows:
%
% 1. Data with the specified file pattern will be searched first in the root folder and each image will be treated
% as an individual for that group. 
% 2. If the images are not in root folder and each subject for that group
% has a sub-folder, then individuals for that group will be selected based on the matching file pattern. First it will
% search subfolder followed by sub-subfolders. If the data is missing error
% will be printed to the Matlab command prompt.
%
% NOTE: 
% 1. After data is selected, the files selected for each group's task
% will be printed to a text file. You can review this text file.
% 2. Keep one image per subject for each group's task so that the toolbox selects the data correctly without reviewing the text file. 

group1_dir = 'F:\Fusion_Data\fmri_fmri\Healthy';
group2_dir = 'F:\Fusion_Data\fmri_fmri\SZ';


%% Data Selection (option 2)
% Files is a cell array of dimensions no. of groups by features. Each cell
% contains a character array of file names.
files = {'F:\Fusion_Data\fmri_fmri\Healthy\tar*.img', 'F:\Fusion_Data\fmri_fmri\Healthy\re*.img';
    'F:\Fusion_Data\fmri_fmri\SZ\tar*.img', 'F:\Fusion_Data\fmri_fmri\SZ\re*.img'};

%% Scale components to z-scores
% Options are 'No' and 'Yes'
z_scores = 'Yes';

%% Number of components
numComp = 10;

%% Type of PCA.
% Options are standard, reference, cca, mcca and mccar.
% NOTE: 
% 1. If you selected reference or mccar as the type of data reduction then enter the
% groups information in the reference variable.
% 2. If you have selected cca or mcca as type of data reduction, specify cca options in variable cca_opts.
type_pca = 'standard';

%% Reference vector to distinguish groups if you selected reference as the type of PCA.
% In the example data there are 30 subjects, 1 and -1 are used as flags for
% Healthy and patients. Leave it as empty if you are not sure about this.
% NOTE: The length of this vector must equal the total number of subjects.
reference = [ones(1, 15), -ones(1, 15)];

%% CCA options
% These options will be used only if you use cca as data reduction
% strategy. Length of principal components must equal the no. of features. These
% numbers must be >= numComp.
cca_opts.numPC = [20, 20];

%% Type of ICA analysis:
% Options are 'average' and 'icasso'
type_ica = 'icasso';

%% Number of times you want to run ICA
num_ica_runs = 10;

%% ICA Algorithm
% Available algorithms are Infomax, FastICA, ...
algorithm = 2;

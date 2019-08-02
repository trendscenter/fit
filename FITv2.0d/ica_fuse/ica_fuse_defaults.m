function ica_fuse_defaults
%% DATA FUSION DEFAULTS
%

%% COLOR DEFAULTS
global FIG_BG_COLOR; % Figure background color
global FIG_FG_COLOR; % Figure foreground color
global UI_BG_COLOR; % background color for user interface controls except pushbutton
global BUTTON_BG_COLOR; % background color for pushbutton
global UI_FG_COLOR; % Foreground color for controls except pushbutton
global BUTTON_FG_COLOR; % Foreground color for pushbutton
global AX_COLOR; % axes color
global LEGEND_COLOR; % LEGEND Text COLOR for Matlab 7.0
global LEGEND_EDGE_COLOR; % Legend edge color
global HELP_FG_COLOR;

%% FONT DEFAULTS
global UI_FONT_NAME; % font name
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size

%% SCREEN SIZE
global SCREENSIZE;
global WSCREEN;
global MIN_SCREEN_DIM_IN_PIXELS;
global PERCENT_SCREEN_OCUPIED;

% VARS FOR READING DATA
global AVAILABLE_MODALITIES;

%% MAT FILES
global FUSION_INFO_MAT_FILE;

%% Optimization of features
global OPTIMIZE_FEATURES;

%% Subject file
global FUSION_SUBJECT_MAT_FILE;

%% Data reduction file name
global DATA_REDUCTION_FILE;

%% ICA file
global ICA_FILE;

%% Back reconstruction file
global BACK_RECONSTRUCT_FILE;

%% Scale components file
global SCALE_COMP_FILE;

%%%%%%% TEXT FILES %%%%%%%
global SELECTED_DATA_TXTFILE;

%% Filters
global MRI_DATA_FILTER;
global EEG_DATA_FILTER;
global BEH_DATA_FILTER;

%% DETREND NUMBER
global DETREND_NUMBER;

global FWHM_VALUE;

%% Image conventions
global JOINT_COMPONENT_NAMING;

%% Display defaults
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGE_VALUES;
global IMAGES_PER_FIGURE;
global ANATOMICAL_PLANE;
global ANATOMICAL_FILE;
global USE_DEFAULT_SLICES;
global SLICESINMM;
global COLORMAPFILE;
global DIVERGENCE_PARAMETERS;

%% CRITICAL PARAMETERS
global EEG_DATA_INDICES;
global NUM_RUNS_ICA;
global STANDARDIZE_SUBJECTS;

% Variable for changing sign of components during scaling
global FLIP_SIGN_COMPONENTS;

%% Histogram defaults
global Z_THRESHOLD_HISTOGRAM;
global NUM_BINS;

%% Remove mean before doing PCA
global RM_PCA;

%% PARALLEL ICA Fusion Defaults
global PARALLEL_ICA_INFO_MAT_FILE;
global PARALLEL_ICA_SEL_DATA_TXT_FILE;

global PARALLEL_ICA_PCA_FILE;
global PARALLEL_ICA_ICA_FILE;

% Parallel ICA modalities
global PARALLEL_ICA_MODALITIES;

% Gene data filter
global GENE_DATA_FILTER;

% Parallel ICA component naming
global PARALLEL_ICA_COMPONENT_NAMING;

% Variable for flipping analyze images
global FLIP_ANALYZE_IM;

% SNP Z Threshold
global SNP_Z_THRESHOLD;

%%%%%%% End for PARALLEL ICA Fusion Defaults %%%%%%%%%%

%% Open display window
global OPEN_DISPLAY_WINDOW;

%% Option to stack all features
global STACK_ALL_FEATURES;

%% Option to enforce MAT file version.
global ENFORCE_MAT_FILE_VERSION;

%% Number of bar elements per axes (Optimal features)
global NUM_BAR_ELEMENTS_PER_AXES;

%% Talairach defaults
global TALAIRACHDIST;
global TALAIRACHTHRESHOLD;

%% Default MCCA cost func
global MCCA_COST_FUNC;

%% Default MCCAR Penalty parameter
global MCCAR_LAMBDA;

%% Color Defaults
FIG_BG_COLOR = [0 0 0]; % Figure back ground color
FIG_FG_COLOR = [1 1 1]; % Figure font color
UI_BG_COLOR = [.2 .2 .2]; % User Interface controls background color except push button
UI_FG_COLOR = [1 1 1]; % Foreground color for UI except push button
AX_COLOR = [0, 0, 0]; % Axes color
LEGEND_COLOR = [1 1 1]; % Legend color for figures on Matlab 7

%% Button background and font colors
BUTTON_BG_COLOR = [.2 .2 .2]; %(default is grey)
BUTTON_FG_COLOR = [1 1 1]; %(default is white)
HELP_FG_COLOR = [1 1 0]; % default is yellow

%% Font Size
UI_FONT_NAME = 'times';
UI_FONT_UNITS = 'points';
UI_FONT_SIZE = 12;

%% Screen size
SCREENSIZE = get(0, 'ScreenSize');
WSCREEN = [SCREENSIZE(3)/1280, SCREENSIZE(4)/960, SCREENSIZE(3)/1280, SCREENSIZE(4)/960];
MIN_SCREEN_DIM_IN_PIXELS = min([SCREENSIZE(3) SCREENSIZE(4)]);
PERCENT_SCREEN_OCUPIED = .9;

%% MAT Files
% SUBJECT FILE
FUSION_SUBJECT_MAT_FILE = '_dataInfo';

% FUSION MAT FILE
FUSION_INFO_MAT_FILE = '_ica_fusion';

% DATA REDUCTION
DATA_REDUCTION_FILE = '_pca_';

% ICA FILE
ICA_FILE = '_ica_';

% Back Reconstruction file
BACK_RECONSTRUCT_FILE = '_br_';

% Scale component file
SCALE_COMP_FILE = '_sc_';

% TEXT Files
SELECTED_DATA_TXTFILE = '_selected_data';


%% AVAILABLE MODALITIES
AVAILABLE_MODALITIES = {'fmri', 'eeg', 'smri', 'gene', 'behavioral'};

%% Files Filter. Same variable used for writing joint components.
% MRI data filter. Options are *.img and *.nii
MRI_DATA_FILTER = '*.img'; % MRI DATA FILTER

% Format for writing EEG data is *.asc (Ascii)
EEG_DATA_FILTER = '*.asc'; % EEG DATA FILTER

% Format for selecting behavioral data 
BEH_DATA_FILTER = '*.asc'; 

%% OPTIMIZATION OF FEATURES
% OPTIONS are 'no' and 'yes'
OPTIMIZE_FEATURES = 'Yes';

FWHM_VALUE = [5.2, 5.2, 5.2];

DETREND_NUMBER = 0; % DETREND NUMBER

JOINT_COMPONENT_NAMING = '_joint_comp_ica_';

%% Display Defaults
CONVERT_TO_Z = 'Yes';
Z_THRESHOLD = 1.5;
IMAGE_VALUES = 'Positive and Negative';
IMAGES_PER_FIGURE = '1';
% Anatomical plane like 'Axial', 'Sagittal', 'Coronal'
ANATOMICAL_PLANE = 'Axial';
% Anatomical file
ANATOMICAL_FILE = which('ch2bet.nii'); %which('nsingle_subj_T1_2_2_5.nii'); 
% Anatomical SLICES in MM
USE_DEFAULT_SLICES = 1;
SLICESINMM = '-40:4:72';
COLORMAPFILE = which('ica_fuse_colors.mat');

%% Default divergence name and number
% Options for name are: 'kl', 'j', 'alpha', 'renyi'
% KL and J divergence doesn't need number
DIVERGENCE_PARAMETERS = {'renyi', 2};

%% CRITICAL PARAMETERS
EEG_DATA_INDICES = (50:500); % EEG data indices
NUM_RUNS_ICA = 1; % Number of ica runs

%%%%%%%%%%%%% Defaults specific to type of fusion %%%%%
STANDARDIZE_SUBJECTS = struct;

% Defaults specific to fmri-fmri fusion
STANDARDIZE_SUBJECTS.fmri_fmri = 'no';

% Defaults specific to eeg-fmri fusion
STANDARDIZE_SUBJECTS.eeg_fmri = 'yes';

% Defaults specific to eeg-eeg fusion
STANDARDIZE_SUBJECTS.eeg_eeg = 'yes';

% Defaults specific to smri-smri fusion
STANDARDIZE_SUBJECTS.smri_smri = 'no';

% Defaults specific to fmri-smri fusion
STANDARDIZE_SUBJECTS.fmri_smri = 'no';

% Defaults specific to eeg-smri fusion
STANDARDIZE_SUBJECTS.eeg_smri = 'yes';

%% Z Threshold
Z_THRESHOLD_HISTOGRAM = 1.5;

%% Number of bins
NUM_BINS = 50;

%% Variable for changing sign of components during scaling
% 1 - Change sign of components w.r.t mean of data
% 0 - Don't change sign of components.
FLIP_SIGN_COMPONENTS = 1;

%% Remove mean while calculating eigen values during PCA step
% 1 - Remove mean
% 0 - Don't remove mean
RM_PCA = 0;

%% Variable for flipping analyze images
FLIP_ANALYZE_IM = 1;

%% Parallel ICA defaults
PARALLEL_ICA_MODALITIES = {'fmri', 'gene', 'smri', 'eeg', 'behavioral'};

% Parallel ICA fusion info MAT file
PARALLEL_ICA_INFO_MAT_FILE = '_para_ica_fusion';

% Parallel ICA selected data text file
PARALLEL_ICA_SEL_DATA_TXT_FILE = '_para_ica_sel_data';

% DATA REDUCTION
PARALLEL_ICA_PCA_FILE = 'para_ica_pca';

% ICA FILE
PARALLEL_ICA_ICA_FILE = 'para_ica_ica';

% Gene data filter
GENE_DATA_FILTER = '*.asc';

% Component namings
PARALLEL_ICA_COMPONENT_NAMING = 'para_ica_comp_';

% Z-threshold for SNPS
SNP_Z_THRESHOLD = 2.5;

%% Open display window
OPEN_DISPLAY_WINDOW = 1;

%% STACK ALL FEATURES
% Options are 0 and 1
% 0 - Don't stack all features.
% 1 - Stack all features.
%
% NOTE: If the number of features is less than 3, all features are stacked
% together whereas for features greater than 3, component images will not be written if you set this variable to 1. However, you will be
% still able to determine optimal features and plot histograms.
STACK_ALL_FEATURES = 1;

%% Enforce MAT files versioning for MATLAB versions greater than 6.5. Use
% the correct option. For more help on version compatibility, please check
% MATLAB save command options for MAT files.
ENFORCE_MAT_FILE_VERSION = '-v6';


%% Number of bar elements per axes. This is used while plotting results of
% optimal features.
NUM_BAR_ELEMENTS_PER_AXES = 10;


%% Talairach defaults %
% Distance between contiguous voxels in mm
TALAIRACHDIST = 4;
% Threshold. This pulls out positive and negative areas by applying data
% > abs(threshold) and -data > abs(threshold)
TALAIRACHTHRESHOLD = 3.5;


%% DEFAULT MCCA COST FUNCTION (N >=3 modalities)
% Options are ssqcor, genvar, maxvar
MCCA_COST_FUNC = 'genvar';

%% MCCA Reference Penalty PARAMETERS
MCCAR_LAMBDA = 0.98;

function [inputText] = ica_fuse_define_parameters_paraICA(inputFile)
% The input parameters are as follows:
% To error check each parameter put an additional field errorCheck that
% defines what criteria must be checked
%
% Note: read_from_input_file field states that field can be read from input file
% if provided

ica_fuse_defaults;
global NUM_RUNS_ICA;

if ~exist('inputFile', 'var')
    inputFile = [];
end

% Parameter 1
numParameters = 1;

% Output prefix
inputText(numParameters).promptString = 'Enter Name(Prefix) Of Output Files';
inputText(numParameters).answerString = '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'prefix'; % tag creates the field in input structure
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 1;
inputText(numParameters).executeCallback = 1;

numParameters = numParameters + 1;

% Data selection
inputText(numParameters).promptString = 'Have You Selected The Data Files?';
inputText(numParameters).answerString = 'Select';
inputText(numParameters).uiType = 'pushbutton';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'dataInfo';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;

numParameters = numParameters + 1;

% Mask Selection
inputText(numParameters).promptString =  'What Mask Do You Want To Use?';
inputText(numParameters).answerString =  char('Default Mask', 'Select Mask');
inputText(numParameters).tag = 'maskFile';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 1;
inputText(numParameters).executeCallback = 1;


numParameters = numParameters + 1;

% Mask Selection
inputText(numParameters).promptString =  'Select Type Of Pre-processing The Data?';
inputText(numParameters).answerString =  char('None', 'Z-scores');
inputText(numParameters).tag = 'preproc_type';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;


numParameters = numParameters + 1;


% MDL estimation
inputText(numParameters).promptString =  'Do You Want To Estimate The Number Of Independent Components?';
inputText(numParameters).answerString =  char('No', 'Yes');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'estimate_components';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;


numParameters = numParameters + 1;

% Number of Independent Components
inputText(numParameters).promptString =  'Number of PC for features';
inputText(numParameters).answerString =  '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'numComp';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;
% 
% numParameters = numParameters + 1;
% 
% % Number of Independent Components
% inputText(numParameters).promptString =  'Number of PC for feature 2';
% inputText(numParameters).answerString =  '';
% inputText(numParameters).uiType = 'edit';
% inputText(numParameters).answerType = 'numeric';
% inputText(numParameters).tag = 'modality2_numComp';
% inputText(numParameters).enable = 'inactive';
% inputText(numParameters).value = 1;
% inputText(numParameters).read_from_file = 1;
% inputText(numParameters).executeCallback = 1;

numParameters = numParameters + 1;


% Type of parallel ICA
inputText(numParameters).promptString = 'Select type of parallel ICA';
inputText(numParameters).answerString = char('AA', 'AS', 'AA-ref');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'type_parallel_ica';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 1;
inputText(numParameters).executeCallback = 1;


numParameters = numParameters + 1;


% Type of PCA
inputText(numParameters).promptString = 'Select type of PCA';
inputText(numParameters).answerString = char('Reference', 'Standard');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'type_pca';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 1;
inputText(numParameters).executeCallback = 1;

numParameters = numParameters + 1;

% Type of ICA
inputText(numParameters).promptString = 'Select Type Of ICA';
inputText(numParameters).answerString = char('Average', 'ICASSO');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'type_ica';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;

numParameters = numParameters + 1;

% No. of ICA runs
inputText(numParameters).promptString = 'Number of times ICA will run';
inputText(numParameters).answerString = num2str(NUM_RUNS_ICA);
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'num_ica_runs';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;


% Read parameters from file and set answer string
if ~isempty(inputFile)
    [inputText] = ica_fuse_checkInputFile(inputFile, inputText);
end
% end for checking input file
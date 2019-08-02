function [inputText] = ica_fuse_define_parameters_mcca(inputFile)
% The input parameters are as follows:
% To error check each parameter put an additional field errorCheck that
% defines what criteria must be checked
%
% Note: read_from_input_file field states that field can be read from input file
% if provided

ica_fuse_defaults;
global NUM_RUNS_ICA;
global DATA_REDUCTION_TYPE;
global ICA_ALGORITHM_NAME;

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

% Normalization type
inputText(numParameters).promptString =  'How do you want to normalize the data?';
inputText(numParameters).answerString =  char('Default', 'Norm 2', 'Std', 'None');
inputText(numParameters).tag = 'normalize';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 1;
inputText(numParameters).executeCallback = 1;

numParameters = numParameters + 1;


% MDL estimation
inputText(numParameters).promptString =  'Do You Want To Estimate The Number Of Components?';
inputText(numParameters).answerString =  char('No', 'Yes');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'estimate_components';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;

numParameters = numParameters + 1;

% Scale components
inputText(numParameters).promptString =  'How Do You Want To Scale Components?';
inputText(numParameters).answerString =  char('Data-Units(eg. EEG-mV)', 'Z-scores');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'z_scores';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 1;
inputText(numParameters).executeCallback = 1;


numParameters = numParameters + 1;

% Number of PC1
inputText(numParameters).promptString =  'Enter No. Of Principal Components For Features In A Vector';
inputText(numParameters).answerString =  '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'numPC';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;


% numParameters = numParameters + 1;
%
% % Number of PC2
% inputText(numParameters).promptString =  'Enter No. Of Principal Components For Feature 2';
% inputText(numParameters).answerString =  '';
% inputText(numParameters).uiType = 'edit';
% inputText(numParameters).answerType = 'numeric';
% inputText(numParameters).tag = 'numPC2';
% inputText(numParameters).enable = 'off';
% inputText(numParameters).value = 1;
% inputText(numParameters).read_from_file = 0;
% inputText(numParameters).executeCallback = 1;


numParameters = numParameters + 1;

% Number of Independent Components
inputText(numParameters).promptString =  'Number of canonical variates';
inputText(numParameters).answerString =  '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'numComp';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 1;
inputText(numParameters).executeCallback = 1;


% Read parameters from file and set answer string
if ~isempty(inputFile)
    
    inputText = ica_fuse_checkInputFile(inputFile, inputText);
    
end
% end for checking input file
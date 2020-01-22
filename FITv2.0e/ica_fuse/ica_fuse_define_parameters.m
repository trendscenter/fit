function [inputText] = ica_fuse_define_parameters(inputFile)
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

% Number of Independent Components
inputText(numParameters).promptString =  'Number of Independent Components';
inputText(numParameters).answerString =  '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'numComp';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 1;
inputText(numParameters).executeCallback = 1;

numParameters = numParameters + 1;


dataRedAlgo = 1;
dataRedList = char('Standard', 'Reference', 'CCA', 'MCCA', 'MCCAR');
if (~isempty(DATA_REDUCTION_TYPE))
    if (~isnumeric(DATA_REDUCTION_TYPE))
        dataRedAlgo = strmatch(lower(DATA_REDUCTION_TYPE), lower(dataRedList), 'exact');
    end
    if (isempty(dataRedAlgo))
        dataRedAlgo = 1;
    end
end

% Type of PCA
inputText(numParameters).promptString = 'Select Type Of Data Reduction';
inputText(numParameters).answerString = dataRedList;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'type_pca';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = dataRedAlgo;
inputText(numParameters).read_from_file = 0;
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

% Number of ICA runs
inputText(numParameters).promptString = 'Number of times ICA will run';
inputText(numParameters).answerString = num2str(NUM_RUNS_ICA);
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'num_ica_runs';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;

numParameters = numParameters + 1;

icaAlgoValue = 1;
algorithmList = ica_fuse_icaAlgorithm;
if (~isempty(ICA_ALGORITHM_NAME))
    if (~isnumeric(ICA_ALGORITHM_NAME))
        icaAlgoValue = strmatch(lower(ICA_ALGORITHM_NAME), lower(algorithmList), 'exact');
    end
    if (isempty(icaAlgoValue))
        icaAlgoValue = 1;
    end
end

% ICA Algorithms
inputText(numParameters).promptString = 'Which ICA/IVA Algorithm Do You Want To Use?';
inputText(numParameters).answerString = algorithmList;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'algorithm';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = icaAlgoValue;
inputText(numParameters).read_from_file = 1;
inputText(numParameters).executeCallback = 1;

% Read parameters from file and set answer string
if ~isempty(inputFile)
    
    [inputText] = ica_fuse_checkInputFile(inputFile, inputText);
    
    % Type of PCA
    keywd = 'type_pca';
    pca_ind = strmatch(keywd, cellstr(char(inputText.tag)), 'exact');
    inputText(pca_ind).read_from_file = 1;
    
    try
        inputData = ica_fuse_read_variables(inputFile, keywd, {'character'});
        tempVal = getfield(inputData, keywd);
        tempInd = strmatch(lower(tempVal), lower(cellstr(inputText(pca_ind).answerString)), 'exact');
        if ~isempty(tempInd)
            inputText(pca_ind).value = tempInd;
        end
    catch
    end
    
end
% end for checking input file
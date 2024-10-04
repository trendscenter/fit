function inputText = ica_fuse_define_parameters_cm_ica
% The input parameters are as follows:
% To error check each parameter put an additional field errorCheck that
% defines what criteria must be checked
%
% Note: read_from_input_file field states that field can be read from input file
% if provided


ica_fuse_defaults;
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
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;

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
inputText(numParameters).executeCallback = 0;

numParameters = numParameters + 1;

% Scale components
inputText(numParameters).promptString =  'Number of PC in the Data Reduction First Step?';
inputText(numParameters).answerString =  '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'numPC1';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;


numParameters = numParameters + 1;

% Scale components
inputText(numParameters).promptString =  'Number of PC/IC in the Data Reduction Second Step?';
inputText(numParameters).answerString =  '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'numPC2';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;

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
inputText(numParameters).promptString = 'Which ICA Algorithm Do You Want To Use?';
inputText(numParameters).answerString = algorithmList;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'algorithm';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = icaAlgoValue;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;

% Read parameters from file and set answer string
if ~isempty(inputFile)
    
    [inputText] = ica_fuse_checkInputFile(inputFile, inputText);
    
end
% end for checking input file
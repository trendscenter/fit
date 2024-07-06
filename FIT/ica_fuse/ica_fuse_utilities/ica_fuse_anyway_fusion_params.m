function inputText = ica_fuse_anyway_fusion_params
% The input parameters are as follows:
% To error check each parameter put an additional field errorCheck that
% defines what criteria must be checked
%
% Note: read_from_input_file field states that field can be read from input file
% if provided


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
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;


numParameters = numParameters + 1;

% Number of Independent Components
inputText(numParameters).promptString =  'Number of Principal/Independent Components (One per feature)';
inputText(numParameters).answerString =  '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'numComp';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;

numParameters = numParameters + 1;

% Number of Independent Components
inputText(numParameters).promptString =  'How Do You Want To Scale Components?';
inputText(numParameters).answerString =  char('Z-scores', 'None');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'scale_components';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;


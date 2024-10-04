function ica_fuse_setup_joint_cmica(fusionFile)
%% Select data for cmICA
%


ica_fuse_defaults;


if (~exist('fusionFile', 'var') || isempty(fusionFile))
    % select output directory
    outputDir = ica_fuse_selectEntry('typeEntity', 'directory', 'title', 'Select output directory for joint cmICA analysis');
else
    outputDir = fileparts(fusionFile);
end

if (isempty(outputDir))
    error('Output directory is not selected for Joint cmICA analysis');
end

cd(outputDir);

inputText = ica_fuse_define_parameters_cm_ica;

[InputHandle] = ica_fuse_plot_controls_fig(inputText, 'JointCmICA', 'on', 'Done', 'Cancel', 1);


% Output prefix callback
%set(findobj(InputHandle, 'tag', 'prefix'), 'callback', {@prefixCallback, InputHandle});

% Set callback for data sets
set(findobj(InputHandle, 'tag', 'dataInfo'), 'callback', {@listCallback, InputHandle});

set(findobj(InputHandle, 'tag', 'Done'), 'callback', {@parseInputs, InputHandle});

set(findobj(InputHandle, 'tag', 'Done'), 'userdata', outputDir);

for nH = 1:length(inputText)
    helpHandle = findobj(InputHandle, 'tag', ['help_', inputText(nH).tag]);
    set(helpHandle, 'callback', {@showHelpDialog, InputHandle});
end


function listCallback(hObject, event_data, handles)
%% Feature list

ica_fuse_defaults;
global UI_FONT_SIZE;

inputData = get(handles, 'userdata');


featureList = '';
try
    featureList = char(inputData.feature_name);
catch
end


InputHandle = ica_fuse_getGraphics('Enter Features', 'normal', 'features_sel', 'on');
set(InputHandle, 'menubar', 'none');
controlWidth = 0.52;
promptHeight = 0.05;
promptWidth = controlWidth;
listboxHeight = controlWidth; listboxWidth = controlWidth;
xOffset = 0.02; yOffset = promptHeight; yPos = 0.92;
okWidth = 0.12; okHeight = promptHeight;

dropDownWidth = 0.4;

%% Design dropdown box
promptPos = [0.25 - 0.5*dropDownWidth, yPos - 0.5*yOffset, 0.4, promptHeight];

yPos = yPos - promptPos(4) - yOffset;

%% Features text and listbox
promptPos = [0.5 - 0.5*controlWidth, yPos - 0.5*yOffset, promptWidth, promptHeight];

textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Features For Joint cmICA', 'tag', ...
    'prompt_features', 'fontsize', UI_FONT_SIZE - 1);
ica_fuse_wrapStaticText(textH);
listboxXOrigin = promptPos(1);
listboxYOrigin = promptPos(2) - yOffset - listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', featureList, 'tag', ...
    'features_listbox', 'fontsize', UI_FONT_SIZE - 1, 'min', 0, 'max', 1, 'value', 1);

addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_features_button', 'fontsize',...
    UI_FONT_SIZE - 1, 'callback', {@addFeatures, InputHandle, handles});
ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_features_button', 'fontsize',...
    UI_FONT_SIZE - 1, 'callback', {@removeFeatures, InputHandle, handles});


promptPos = listboxPos;

%% Add cancel, save and run buttons
cancelPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
cancelPos(2) = cancelPos(2) - 0.5*cancelPos(4);
ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', cancelPos, 'string', 'Cancel', 'tag', 'cancel_button', 'fontsize',...
    UI_FONT_SIZE - 1, 'callback', 'delete(gcbf);');

okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'create_button', 'fontsize',...
    UI_FONT_SIZE - 1, 'callback', {@featuresChk, InputHandle, handles});


function dataSelectCallback(hObject, handles)
%% Data select callback
%

% Parameter 1
numParameters = 1;

% Feature name
inputText(numParameters).promptString = 'Enter Feature Name';
inputText(numParameters).answerString = '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'feature_name'; % tag creates the field in input structure
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;

numParameters = numParameters + 1;

% Modality selection
inputText(numParameters).promptString = 'Select Modality. 1st should be fMRI type';
inputText(numParameters).answerString = char('fMRI', 'DTI');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'modality';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;

numParameters = numParameters + 1;

% Scale components
inputText(numParameters).promptString =  'Have You Selected Data';
inputText(numParameters).answerString =  'Select';
inputText(numParameters).uiType = 'pushbutton';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'files';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;


numParameters = numParameters + 1;

% Scale components
inputText(numParameters).promptString =  '*Select Mask';
inputText(numParameters).answerString =  'Mask';
inputText(numParameters).uiType = 'pushbutton';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'mask';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;


paramsH = ica_fuse_plot_controls_fig(inputText, 'DataSelect', 'on', 'Done', 'Cancel', 1, -0.2);

tempH = ica_fuse_uicontrol('parent', paramsH, 'units', 'normalized', 'style', 'text', ...
    'position', [ 0.0100    0.7800    0.9900    0.1000], ...
    'string', '*Note: 1st feature mask may need to be (NIfTI) nii. 2nd whole brain mask should match your DTI indeces ', ...
    'tag', ['prompt', 'dummy'], 'panel_tag', ['panel_text', 'dymmy2'], 'fontsize', 11);
tempH = ica_fuse_uicontrol('parent', paramsH, 'units', 'normalized', 'style', 'text', ...
    'position', [ 0.0100    0.700    0.9900    0.1000], ...
    'string', 'eg., desc-wmstreamline-wholebrain_res-5x5x5_mask.nii', ...
    'tag', ['prompt', 'dummy'], 'panel_tag', ['panel_text', 'dymmy2'], 'fontsize', 11);

for nH = 1:length(inputText)
    helpHandle = findobj(paramsH, 'tag', ['help_', inputText(nH).tag]);
    set(helpHandle, 'callback', {@showHelpDialog, paramsH});
end


% Set callback for data sets
set(findobj(paramsH, 'tag', 'files'), 'callback', {@filesCallback, paramsH});

% Set callback for done
set(findobj(paramsH, 'tag', 'Done'), 'callback', {@getInputs, paramsH, handles});

% Set callback for mask
set(findobj(paramsH, 'tag', 'mask'), 'callback', {@maskCallback, paramsH});

drawnow;

try
    waitfor(paramsH);
catch
end

try
    delete(paramsH);
catch
end

function getInputs(hObject, event_data, handles, InputH)
%% Get inputs
%

featureH = findobj(handles, 'tag', 'feature_name');
featureName = get(featureH, 'string');
featureName = strtrim(deblank(featureName));

if (isempty(featureName))
    error('Feature name is not selected');
end

filesH = findobj(handles, 'tag', 'files');
files = get(filesH, 'userdata');

if (isempty(files))
    error('Files are not selected');
end

inputData = get(InputH, 'userdata');

try
    inds = strmatch(lower(featureName), lower(cellstr(char(inputData.feature_name))), 'exact');
    if (~isempty(inds))
        error(['Feature ', featureName , ' already exists. Use a different name']);
    end
catch
end

modalityH = findobj(handles, 'tag', 'modality');
modalities = cellstr(get(modalityH, 'string'));
modalityVal = get(modalityH, 'value');
modalitySel = modalities{modalityVal};

maskH = findobj(handles, 'tag', 'mask');
mask = get(maskH, 'userdata');

if (isempty(inputData))
    % Store input parameters
    inputData.modality = modalitySel;
    inputData.feature_name = featureName;
    inputData.files = files;
    inputData.mask = mask;
else
    %     inds = strmatch(lower(featureName), lower(cellstr(char(inputData.feature_name))), 'exact');
    %     if (~isempty(inds))
    %         error(['Feature ', featureName , ' already exists. Use a different name']);
    %     end
    % Append data
    featureLen = length(inputData) + 1;
    inputData(featureLen).modality = modalitySel;
    inputData(featureLen).feature_name = featureName;
    inputData(featureLen).files = files;
    inputData(featureLen).mask = mask;
end

set(InputH, 'userdata', inputData);

try
    delete(handles);
catch
end


function filesCallback(hObject, event_data, handles)

featureH = findobj(handles, 'tag', 'feature_name');
featureName = get(featureH, 'string');
featureName = strtrim(deblank(featureName));

if (isempty(featureName))
    error('Feature name is not selected');
end

titleStr = ['Feature ', featureName];

files = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', 'title', ['Select all subjects for ', titleStr], 'filter', '*.nii;*.dot');

if (isempty(files))
    error([titleStr, ' files are not selected']);
end

set(hObject, 'userdata', files);


function maskCallback(hObject, event_data, handles)
%% Get inputs
%

featureH = findobj(handles, 'tag', 'feature_name');
featureName = get(featureH, 'string');
featureName = strtrim(deblank(featureName));

if (isempty(featureName))
    error('Feature name is not selected');
end

titleStr = ['Feature ', featureName];
mask = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'title', ['Select Mask for ', titleStr], 'filter', '*.nii');

if (isempty(mask))
    error([titleStr, ' mask is not selected']);
end

set(hObject, 'userdata', mask);

function featuresChk(hObject, event_data, figH, handles)
%% Add features
%

ud = get(handles, 'userdata');

if (isempty(ud))
    error('Features data is not selected');
end

selectH = findobj(handles, 'tag', 'dataInfo');
set(selectH, 'foregroundcolor', [0, 1, 0]);


try
    delete(figH);
catch
end

function addFeatures(hObject, event_data, figH, handles)
%% Add features


listH = findobj(figH, 'tag', 'features_listbox');

addH = findobj(figH, 'tag', 'add_features_button');

dataSelectCallback(addH, handles);

drawnow;

ud = get(handles, 'userdata');
set(listH, 'value', 1);
set(listH, 'string', cellstr(char(ud.feature_name)));

function removeFeatures(hObject, event_data, figH, handles)
%% Remove features
%

inputData = get(handles, 'userdata');
listH = findobj(figH, 'tag', 'features_listbox');
val = get(listH, 'value');

if (~isempty(val))
    check = ica_fuse_questionDialog('title', 'Remove Feature', 'textbody', 'Do you want to remove the feature from the list?');
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(inputData.feature_name));
    inputData(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(handles, 'userdata', inputData);
catch
end

function parseInputs(hObject, event_data, handles)
%% Parse INputs


inputData = get(handles, 'userdata');
outputDir = get(hObject, 'userdata');

if (isempty(inputData))
    error('Features data is not selected');
end

output_prefix = strtrim(get(findobj(handles, 'tag', 'prefix'), 'string'));

numPC1 = str2num(strtrim(get(findobj(handles, 'tag', 'numPC1'), 'string')));
numPC2 = str2num(strtrim(get(findobj(handles, 'tag', 'numPC2'), 'string')));
algoList = cellstr(get(findobj(handles, 'tag', 'algorithm'), 'string'));
algoVal = get(findobj(handles, 'tag', 'algorithm'), 'value');

cmICAInfo.outputDir = outputDir;
cmICAInfo.output_prefix = output_prefix;
cmICAInfo.dataInfo = inputData;
cmICAInfo.numPC1 = numPC1;
cmICAInfo.numComp = numPC2;
cmICAInfo.algorithm = algoList{algoVal};

fileName = fullfile(outputDir, [output_prefix, '_joint_cmica_info.mat']);

disp(['Saving session information in file ', fileName, ' ...']);
save(fileName, 'cmICAInfo');
delete(handles);
disp('Done');
fprintf('\n');


function showHelpDialog(hObject, event_data, handles)
% Show help dialog box

inputTag = get(hObject, 'tag');

inputTag = strrep(inputTag, 'help_', '');

switch (inputTag)
    
    case 'prefix'
        
        titleStr = 'Prefix';
        D(1).string = 'Enter a valid variable name as application data and files will be stored with this name.';
        
    case {'files', 'dataInfo'}
        
        titleStr = 'Data';
        D(1).string = 'Select the data for the joint connectivity matrix ICA analysis. After the data selection, a drop down box will appear with options as ''Yes'' and ''No''. No means data can be selected again.';
        
        
    case 'numPC1'
        
        titleStr = 'No. of principal components in the first PCA step';
        D(1).string = 'Enter the no. of principal components you want to extract from the data.';
        
    case 'numPC2'
        
        titleStr = 'No. of independent components';
        D(1).string = 'Enter the no. of indepenent components you want to extract from the data for each modality.';
        
        
    case 'algorithm'
        
        titleStr = 'ICA Algorithm';
        D(1).string = 'Currently, 15 ICA algorithms are available in the toolbox like Infomax, Fast ICA, Erica, Simbec, Evd, Jade Opac, Amuse, SDD ICA, CCICA, Combi, EBM, ERBM, IVA-G, IVA-GGD adn pmjlICA.';
        
    case 'feature_name'
        
        titleStr = 'Feature Name';
        D(1).string = 'Enter feature name';
        
    case 'modality'
        
        titleStr = 'Select modality';
        D(1).string = 'Options are fMRI and DTI';
        
    case 'mask'
        
        titleStr = 'Mask';
        D(1).string = 'Mask used for the analysis. For fMRI modality, mask is optional whereas for DTI you need to select mask for the probability tracts';
        
        
        
    otherwise
        
        titleStr = 'Unknown';
        D(1).string = 'Unknown option specified';
        
end
% End for checking

msgStr = str2mat(D.string);
ica_fuse_dialogBox('title', titleStr, 'textBody', msgStr, 'textType', 'large');

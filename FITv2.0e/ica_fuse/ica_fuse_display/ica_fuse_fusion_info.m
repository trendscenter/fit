function ica_fuse_fusion_info(fusionFile)
% Display analysis information
% Plot Parameter, Analysis, Output files

ica_fuse_defaults;
global FUSION_INFO_MAT_FILE;


if ~exist('fusionFile', 'var')
    fusionFile = ica_fuse_selectEntry('title', 'Select joint ICA fusion information file', ...
        'typeEntity', 'file', 'typeSelection', 'Single', 'filter', ['*', FUSION_INFO_MAT_FILE, '*.mat']);
end

if isempty(fusionFile)
    error(['Fusion information file is not selected']);
end

drawnow;

% Load fusion file
load(fusionFile);

% Make sure the selected file is a valid fusion parameter file
if ~exist('fusionInfo', 'var')
    error(['Selected file: ', fusionFile, ' is not a valid fusion parameter file']);
end


if ~isfield(fusionInfo, 'run_analysis')
    error('Please run the analysis to view the fusion analysis information');
end

[outputDir, fileN, extn] = fileparts(fusionFile);
if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

%% Check this for backward compatibility
if ~isfield(fusionInfo.run_analysis, 'all_comb')
    % Assuming all the features are stacked
    stackAllFeatures = 1;
    optimalFeatures = 0;
    if (length(fusionInfo.run_analysis.pcaFiles) > 1)
        optimalFeatures = 1;
    end
    % Get all the combinations
    fusionInfo.run_analysis.all_comb = ica_fuse_get_combinations(fusionInfo.run_analysis.numFeatures, stackAllFeatures, optimalFeatures);
end

fusionInfo.run_analysis.good_inds = ica_fuse_good_cells(fusionInfo.run_analysis.all_comb);

figTag = 'Fusion Info';

handleVec = findobj('tag', figTag);
for nn = 1:length(handleVec)
    delete(handleVec(nn));
end

inputHandle = ica_fuse_getGraphics(figTag, 'normal', figTag, 'on');
set(inputHandle, 'menubar', 'none');

% help on analysis info
fitHelpTitle = uimenu('parent', inputHandle, 'label', 'FIT-Help');
fitHelpMenu = uimenu(fitHelpTitle, 'label', 'Analysis Info', 'callback', ...
    'ica_fuse_openHTMLHelpFile(''fit_analysis_info.htm'');');

xOffset = 0.025; yOffset = 0.03;

% Button width
buttonWidth = (1 - xOffset*4) / 3; buttonHeight = 0.05;

% Button 1: Parameter Info
buttonPos = [xOffset, yOffset, buttonWidth, buttonHeight];
paramH = ica_fuse_uicontrol('parent', inputHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', buttonPos, 'string', 'Parameter Info', 'tag', 'parameter_info', 'callback', ...
    {@paramInfoCallback, inputHandle});

% Button 2: Analysis Info
buttonPos(1) = buttonPos(1) + buttonPos(3) + xOffset;
analysisH = ica_fuse_uicontrol('parent', inputHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', buttonPos, 'string', 'Analysis Info', 'tag', 'analysis_info', 'callback', ...
    {@analysisInfoCallback, inputHandle});

% Button 3: Output Files
buttonPos(1) = buttonPos(1) + buttonPos(3) + xOffset;
outputFilesH = ica_fuse_uicontrol('parent', inputHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', buttonPos, 'string', 'Output Files Info', 'tag', 'output_files_info', 'callback', ...
    {@outputFilesCallback, inputHandle});


%%% Plot listbox
% Listbox width and height
listboxWidth = 1 - 2*xOffset;
listboxOrigin = buttonPos(2) + buttonPos(4) + yOffset;
listboxHeight = 1 - yOffset - listboxOrigin;
listboxPos = [xOffset, listboxOrigin, listboxWidth, listboxHeight];

listboxH = ica_fuse_uicontrol('parent', inputHandle, 'units', 'normalized', 'style', ...
    'listbox', 'position', listboxPos, 'string', [], 'tag', 'analysis_info_listbox');

% Apply conditions for dialog box differently for different platforms
if ispc
    set(listboxH, 'enable', 'inactive');
else
    set(listboxH, 'enable', 'on');
end

set(listboxH, 'min', 0, 'max', 2);

% make no selection
set(listboxH, 'value', []);

set(inputHandle, 'userdata', fusionInfo);

paramInfoCallback(paramH, [], inputHandle);


%%%%%%%%%%%%% FUNCTION CALLBACKS %%%%%%%%%%%%%
function paramInfoCallback(hObject, event_data, handles)
%% Parameter information callback
%

set(handles, 'pointer', 'watch');

listboxHandle = findobj(handles, 'tag', 'analysis_info_listbox');

%% Fusion information
fusionInfo = get(handles, 'userdata');

dataInfo = fusionInfo.run_analysis.dataInfo;

%% Mask string
if isempty(fusionInfo.setup_analysis.maskFile)
    maskStr{1} = ['Mask: Default mask created from data'];
else
    maskStr = cell(length(dataInfo(1).feature), 1);
    % Loop over features
    for nn = 1:length(dataInfo(1).feature)
        maskStr{nn} = ['Mask For Feature ', dataInfo(1).feature(nn).name, ': ', fusionInfo.setup_analysis.maskFile{nn}];
    end
    % End loop over features
end

maskStr = str2mat(maskStr);

%% Normalization method
inputText = ica_fuse_define_parameters;
matchIndex = strmatch('normalize', lower(cellstr(str2mat(inputText.tag))), 'exact');
str = inputText(matchIndex).answerString;
clear inputText;
str = deblank(str(fusionInfo.run_analysis.normalize, :));

normalizationStr = ['Normalization Method: ', str];

clear str;

%% Scaling method
z_scores = 0;
if isfield(fusionInfo.run_analysis, 'z_scores')
    z_scores = fusionInfo.run_analysis.z_scores;
end

if z_scores
    str = 'Z scores';
else
    str = 'Percent signal change';
end

scaledMethodStr = ['Scaling Method: ', str];

%% Type of PCA
str = 'Standard';
try
    str = fusionInfo.run_analysis.type_pca;
catch
end

typePCAStr = ['Type of PCA: ', upper(str)];

clear str;

%% ICA Algorithm str
icaAlgo = ica_fuse_icaAlgorithm;
str = deblank(icaAlgo(fusionInfo.run_analysis.algorithm, :));
ICAAlgoStr = ['ICA Algorithm: ', str];
clear str;

icaType = 'Average';
try
    icaType = upper(fusionInfo.run_analysis.type_ica);
catch
end

if (strcmpi(icaType, 'average'))
    icaType = 'Average over runs';
end

icaTypeStr = ['ICA Analysis Type: ', icaType];

% %% String for listbox
listString = {['Output Prefix: ', fusionInfo.run_analysis.prefix]; ...
    ['Number Of Groups: ', num2str(fusionInfo.run_analysis.numGroups)]; ...
    ['Number Of Features: ', num2str(fusionInfo.run_analysis.numFeatures)];
    ['Group Names: ', ica_fuse_formatStr(str2mat(dataInfo.name))]; ...
    ['Feature Names: ', ica_fuse_formatStr(str2mat(dataInfo(1).feature.name))]; ...
    ['Feature Modality: ', ica_fuse_formatStr(str2mat(dataInfo(1).feature.modality))]; ...
    ['Number Of Subjects: ', ica_fuse_formatStr(num2str(fusionInfo.run_analysis.numSubjects(:)))]; ...
    ['Number Of Components: ', num2str(fusionInfo.run_analysis.numComp)]; ...
    maskStr; normalizationStr; scaledMethodStr; typePCAStr; ICAAlgoStr; icaTypeStr};

listString = cellstr(str2mat(listString));
%% Do text wrap
listString = textwrap(listboxHandle, listString);
set(listboxHandle, 'string', listString);

set(handles, 'pointer', 'arrow');

function analysisInfoCallback(hObject, event_data, handles)
%% Analysis information callback
%

set(handles, 'pointer', 'watch');

listboxHandle = findobj(handles, 'tag', 'analysis_info_listbox');

% fusion information
fusionInfo = get(handles, 'userdata');

good_inds = fusionInfo.run_analysis.good_inds;

%% Normalize string
normalizeStr = ['Normalization parameters: ', ica_fuse_formatStr(num2str(fusionInfo.run_analysis.featureNormPara(:)))];

%% PCA files string
combinationNames = str2mat(fusionInfo.run_analysis.pcaFiles(good_inds).combinationName);
formatStr = ica_fuse_formatStr(combinationNames);
pcaFilesStr = ['Principal component analysis was run on combination/combinations: ', formatStr];

%% ICA files string
icaFilesStr = ['Independent component analysis was run on combination/combinations: ', formatStr];

%% Backreconstruction files string
brFilesStr = ['Components are back-reconstructed for combination/combinations: ', formatStr];

listString = {normalizeStr; pcaFilesStr; icaFilesStr; brFilesStr};

%% Scale components files
if (good_inds(1))
    listString{end + 1} = ['Components are scaled for combination: ', fusionInfo.run_analysis.scaleCompFiles(1).combinationName];
end

listString = textwrap(listboxHandle, listString);
set(listboxHandle, 'string', listString);

set(handles, 'pointer', 'arrow');

function outputFilesCallback(hObject, event_data, handles)
%% Output files information callback
%

set(handles, 'pointer', 'watch');

listboxHandle = findobj(handles, 'tag', 'analysis_info_listbox');

% fusion information
fusionInfo = get(handles, 'userdata');

% Good cells
good_inds = find(fusionInfo.run_analysis.good_inds ~= 0);

%% PCA MAT files
listString(1).string = 'PCA MAT file information: ';
% Loop over good indices
for nn = good_inds
    listString(length(listString) + 1).string = ['Combination ', fusionInfo.run_analysis.pcaFiles(nn).combinationName, ': ', ...
        deblank(fusionInfo.run_analysis.pcaFiles(nn).name(1, :))];
end
% End loop over good indices
listString(length(listString) + 1).string = '';

%% ICA MAT files
listString(length(listString) + 1).string = 'ICA MAT file information: ';
% Loop over good indices
for nn = good_inds
    listString(length(listString) + 1).string = ['Combination ', fusionInfo.run_analysis.icaFiles(nn).combinationName, ': ', ...
        deblank(fusionInfo.run_analysis.icaFiles(nn).name(1, :))];
end
% End loop over good indices
listString(length(listString) + 1).string = '';

%% Back reconstruction files
listString(length(listString) + 1).string = 'Back reconstruction MAT file information: ';
% Loop over good indices
for nn = good_inds
    listString(length(listString) + 1).string = ['Combination ', ...
        fusionInfo.run_analysis.backReconstructFiles(nn).combinationName, ': ', ...
        deblank(fusionInfo.run_analysis.backReconstructFiles(nn).name(1, :))];
end
% End loop over good indices
listString(length(listString) + 1).string = '';

if (good_inds(1) == 1)
    % If the first feature combination is run
    listString(length(listString) + 1).string = 'Scaling components MAT file information: ';
    listString(length(listString) + 1).string = ['Combination ', ...
        fusionInfo.run_analysis.scaleCompFiles(1).combinationName, ': ', ...
        deblank(fusionInfo.run_analysis.scaleCompFiles(1).name(1, :))];
    listString(length(listString) + 1).string = '';
    listString(length(listString) + 1).string = 'Joint ICA components output file information: ';
    
    % Loop over all output files
    for nn = 1:length(fusionInfo.run_analysis.outputFiles)
        listString(length(listString) + 1).string = ['Feature ', ...
            fusionInfo.run_analysis.outputFiles(nn).feature_name, ': ', ...
            deblank(fusionInfo.run_analysis.outputFiles(nn).name(1, :))];
    end
    % End loop over output files
    listString(length(listString) + 1).string = '';
end

% Set string for listbox
listStr = str2mat(listString.string);
clear listString;

listString = textwrap(listboxHandle, cellstr(listStr));
set(listboxHandle, 'string', listString);

set(handles, 'pointer', 'arrow');

function ica_fuse_groupsNSortingCriteria
% GUI for getting sorting criteria and groups

% Load defaults
ica_fuse_defaults;

global FUSION_INFO_MAT_FILE;
global UI_FONT_NAME;
global UI_FONT_SIZE;
global AX_COLOR;
global COLORMAPFILE;
global FIG_FG_COLOR;
global HISTOGRAM_THRESHOLD;

fontSize = 0.01;
fusionFile = ica_fuse_selectEntry('title', 'Select fusion information file for plotting histograms', ...
    'typeEntity', 'file', 'typeSelection', 'Single', 'filter', ['*', FUSION_INFO_MAT_FILE, '*.mat']);
drawnow;

load(fusionFile);

if ~exist('fusionInfo', 'var')
    error(['Selected file: ', fusionFile, ' is not a valid fusion parameter file']);
end

if ~isfield(fusionInfo, 'run_analysis')
    error('Please run the analysis inorder to plot histograms of individual subjects');
end

[outputDir, fileN, extn] = fileparts(fusionFile);
if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

%%%%%%% Get the required parameters from run_analysis field %%%%%%%%%%%

% Get the information from the fusion file
numComp = fusionInfo.run_analysis.numComp;

% Output prefix
output_prefix = fusionInfo.run_analysis.prefix;

% Get the feature names
featureNames = str2mat(fusionInfo.run_analysis.dataInfo(1).feature.name);

% Get the groups names
groupNames = str2mat(fusionInfo.run_analysis.dataInfo.name);

% Get the ouput files
outputFiles = fusionInfo.run_analysis.outputFiles;

% Get the dataInfo
dataInfo = fusionInfo.run_analysis.dataInfo;

% Get the PCA files
pcaFiles = fusionInfo.run_analysis.pcaFiles;

% Get the ICA Files
icaFiles = fusionInfo.run_analysis.icaFiles;

% Back reconstruct files
backReconstructFiles = fusionInfo.run_analysis.backReconstructFiles;

% Number of groups
numGroups = fusionInfo.run_analysis.numGroups;

% Number of features
numFeatures = fusionInfo.run_analysis.numFeatures;

% voxel indices
maskIndices = fusionInfo.run_analysis.mask_ind;

% Number of subjects
numSubjects = fusionInfo.run_analysis.numSubjects;

featureDataLength = fusionInfo.run_analysis.featureDataLength;
% Feature data length
dataLength = featureDataLength(1).Length;

normalize = fusionInfo.run_analysis.normalize;

clear fusionInfo;


displayTag = 'OptimalFeatures';

% Get divergence name and value from defaults
[div_name, div_value] = ica_fuse_getDivergencePara;

if ~isempty(div_value)
    divStr = ['(', upper(div_name), '(', num2str(div_value), '))'];
else
    divStr = ['(', upper(div_name), ')'];
end

%%%%%% Delete any previous figures%%%%%
displayH = findobj('tag', displayTag);

for ii = 1:length(displayH)
    delete(displayH(ii));
end
%%%%%% end for deleting previous figures of display GUI %%%%

[displayHandle] = ica_fuse_getGraphics('Optimal Features', 'normal', displayTag, 'on');
set(displayHandle, 'menubar', 'none');

numPara = 1;
sortingCriteria(numPara).displayString = 'Two sample t-test on mixing coefficients between groups';
sortingCriteria(numPara).string = 'ttest2';

numPara = numPara + 1;
sortingCriteria(numPara).displayString =  ['Spatial divergence ', divStr, ' between groups'];
sortingCriteria(numPara).string = 'divergence';

clear numPara;


choiceString = str2mat('Select sorting criteria ...', str2mat(sortingCriteria.displayString));

% offsets
xOffset = 0.04; yOffset = 0.035;

popupTextHeight = 0.05; popupTextWidth = 1 - 2*xOffset;

yPos = 0.95;
popupTextPos = [xOffset, yPos - yOffset - popupTextHeight, popupTextWidth, popupTextHeight];

% Plot popup
sortTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupTextPos, 'String', choiceString, 'fontsize', UI_FONT_SIZE - 1);

listboxWidth = 0.5;

%%%%%%%%%%%%%% Groups Listbox %%%%%%%%%%%%%%%%%%%%

% Plot textbox
textBoxWidth = listboxWidth; textboxHeight = 0.05;
textboxPos = [0.5 - 0.5*listboxWidth, popupTextPos(2) - popupTextPos(4) - 2*yOffset, textBoxWidth, textboxHeight];

% Plot groups text
groupsTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Select Groups', 'fontsize', UI_FONT_SIZE - 1);

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
listboxPos(4) = 0.4;
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);

if numGroups >= 2
    groupsListVal = [1:2];
end

% Plot groups listbox
groupsListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', groupNames, 'value', ...
    groupsListVal, 'min', 0, 'max', 2, 'tag', 'selGroups', 'fontsize', UI_FONT_SIZE - 1, 'TooltipString', ...
    'Select groups ...');

%%%%%%%%%%%%%%% End for plotting groups listbox and text %%%%%%%%%%%%


groupListPos = get(groupsListH, 'position');

buttonHeight = 0.05;
buttonWidth = 0.2;

donePos = [0.5 - 0.5*buttonWidth, yOffset + 0.5*buttonHeight, buttonWidth, buttonHeight];

doneButtonH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', donePos, 'string', 'Calculate', 'tag', 'done_button', 'callback', ...
    {@doneCallback, displayHandle}, 'fontsize', UI_FONT_SIZE - 1);

%%%%%%%%%%%%%%%%%%%%% Object Callbacks %%%%%%%%%%%%%%%%%%%%%

function doneCallback(hObject, event_data, handles)
% Done callback


selGroupsHandle = findobj(handles, 'tag', 'selGroups'); % Groups handle

selGroups = get(selGroupsHandle, 'value');

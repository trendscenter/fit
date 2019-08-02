function ica_fuse_plotHistogram
% Plot cross-task histograms for the selected groups.
% Option is provided to plot histograms for both features and group
% components
%
% Input:
% 1. fusionFile - full file path of fusion MAT file


% Load defaults
ica_fuse_defaults;

global FUSION_INFO_MAT_FILE;
global UI_FONT_NAME;
global UI_FONT_SIZE;
global AX_COLOR;
global COLORMAPFILE;
global FIG_FG_COLOR;
global Z_THRESHOLD_HISTOGRAM;

fontSize = 0.01;
% Select fusion file
%if ~exist('fusionFile', 'var')
fusionFile = ica_fuse_selectEntry('title', 'Select fusion information file for plotting histograms', ...
    'typeEntity', 'file', 'typeSelection', 'Single', 'filter', ['*', FUSION_INFO_MAT_FILE, '*.mat']);
drawnow;
%end

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

% Get the feature names
%featureNames = str2mat(fusionInfo.run_analysis.dataInfo(1).feature.name);

% Get the groups names
groupNames = str2mat(fusionInfo.run_analysis.dataInfo.name);

pcaFiles = fusionInfo.run_analysis.pcaFiles;

% Get the dataInfo
dataInfo = fusionInfo.run_analysis.dataInfo;

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
featureNormPara = fusionInfo.run_analysis.featureNormPara;


%% Check this for backward compatibility
if (isfield(fusionInfo.run_analysis, 'all_comb'))
    all_comb = fusionInfo.run_analysis.all_comb;
else
    % Assuming all the features are stacked
    stackAllFeatures = 1;
    optimalFeatures = 0;
    if (length(pcaFiles) > 1)
        optimalFeatures = 1;
    end
    % Get all the combinations
    all_comb = ica_fuse_get_combinations(numFeatures, stackAllFeatures, optimalFeatures);
end

% good cells
good_inds = find(ica_fuse_good_cells(all_comb) ~= 0);

if (~isfield(fusionInfo.run_analysis, 'dims'))
    [dims, voxels] = ica_fuse_getFeatureDIM(mask_ind);
else
    dims = fusionInfo.run_analysis.dims;
    voxels = fusionInfo.run_analysis.voxels;
end

if isfield(fusionInfo.run_analysis, 'allCombNames')
    combinationsListStr = fusionInfo.run_analysis.allCombNames;
    combinationsListStr = combinationsListStr(good_inds);
else
    combinationsListStr = repmat({''}, length(good_inds), 1);
    % Loop over good cells
    for nC = 1:length(combinationsListStr)
        load(fullfile(outputDir, pcaFiles(good_inds(nC)).name), 'combinationName');
        combinationsListStr{nC} = combinationName;
        clear combinationName;
    end
    % End loop over good cells
end

clear fusionInfo;

listSize = [300 400];

threshold = Z_THRESHOLD_HISTOGRAM;

compStr = num2str((1:numComp)');

selectedCombNumber = 1;
if (length(good_inds) > 1)             
    titleStr = 'Select the feature combination';    
    
    %% Open list dialog box
    selectedCombNumber = ica_fuse_listdlg('promptstring', titleStr, 'liststring', combinationsListStr, 'windowstyle', 'normal', 'selectionmode', 'single', 'movegui', ...
        'center', 'title_fig', titleStr);    
    
    if isempty(selectedCombNumber)
        error('Please select the feature combination');
    end        
end

selectedCombStr = combinationsListStr{selectedCombNumber};            

% Absolute index 
selectedCombNumber = good_inds(selectedCombNumber);

% Feature names
featureNames = str2mat(strread(selectedCombStr, '%s', 'delimiter', '&'));

%%%%%%%%%%% Plot Histogram GUI %%%%%%%%

displayTag = 'histogramGUI';


%%%%%% Delete any previous figures%%%%%
displayH = findobj('tag', displayTag);

for ii = 1:length(displayH)
    delete(displayH(ii));
end
%%%%%% end for deleting previous figures of display GUI %%%%

[displayHandle] = ica_fuse_getGraphics('Histogram GUI', 'normal', displayTag, 'on');
set(displayHandle, 'menubar', 'none');


fitHelpTitle = uimenu('parent', displayHandle, 'label', 'FIT-Help');
fitHelpMenu = uimenu(fitHelpTitle, 'label', 'Histogram Plot', 'callback', ...
    'ica_fuse_openHTMLHelpFile(''fit_histogram_plot.htm'');');

% offsets
xOffset = 0.04; yOffset = 0.035;

popupTextHeight = 0.05; popupTextWidth = 0.52;

yPos = 0.95;
popupTextPos = [xOffset, yPos - yOffset - popupTextHeight, popupTextWidth, popupTextHeight];

popupTextString = {'What Histograms Do You Want To View?'};
% Plot popup
histTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', popupTextPos, 'String', popupTextString, 'fontsize', UI_FONT_SIZE - 1);

[popupTextString, newPos] = textwrap(histTextH, popupTextString);

popupTextPos(4) = newPos(4);
popupTextPos(2) = yPos - yOffset - popupTextPos(4);

% Set the wrapped text and the new position
set(histTextH, 'string', popupTextString);
set(histTextH, 'position', popupTextPos);


%%%% Plot Popup control
popupPos = popupTextPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = 0.32;
popupPos(4) = 0.05;
popupPos(2) = popupPos(2) + 0.5*popupTextPos(4) - 0.5*popupPos(4);
%popupString = str2mat('Features', 'Group Components');

popupString = 'Features';

histPopupH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupPos, 'String', popupString, 'fontsize', UI_FONT_SIZE - 1, 'tag', ...
    'histogram_criteria');

listboxWidth = 0.35;

%%%%%%%%%%%%%% Groups Listbox %%%%%%%%%%%%%%%%%%%%

% Plot textbox
textBoxWidth = listboxWidth; textboxHeight = 0.05;
textboxPos = [xOffset, popupPos(2) - popupPos(4) - 2*yOffset, textBoxWidth, textboxHeight];

% Plot groups text
groupsTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Groups', 'fontsize', UI_FONT_SIZE - 1);

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
listboxPos(4) = 0.22;
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);

if (numGroups >= 2)
    groupsListVal = (1:2);
else
    groupsListVal = 1;
end

% Plot groups listbox
groupsListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', groupNames, 'value', ...
    groupsListVal, 'min', 0, 'max', 2, 'tag', 'selGroups', 'fontsize', UI_FONT_SIZE - 1, 'TooltipString', ...
    'Select groups ...');

%%%%%%%%%%%%%%% End for plotting groups listbox and text %%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% Features Text and Listbox %%%%%%%%%%%%%%%%%%%%%%
textboxPos(1) = 1 - textboxPos(3) - 2*xOffset;
textboxPos(3) = listboxWidth;

% Plot Feature text
featureTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Feature');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);


if (size(featureNames, 1) >= 2)
    featuresListVal = (1:2);
else
    featuresListVal = 1;
end


% Plot component listbox
featureListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', featureNames, 'value', ...
    featuresListVal, 'min', 0, 'max', 2, 'tag', 'selFeatures', 'fontsize', ...
    UI_FONT_SIZE - 1, 'TooltipString', 'Select features ...');


groupListPos = get(groupsListH, 'position');


textboxPos(2) = groupListPos(2) - 2*yOffset;
textboxPos(1) = groupListPos(1);

% Plot Component text
compTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Component');


listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);

listboxPos(1) = textboxPos(1);

% Plot component listbox
compListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', compStr, 'value', ...
    1, 'min', 0, 'max', 1, 'tag', 'selComponent', 'fontsize', ...
    UI_FONT_SIZE - 1, 'TooltipString', 'Select a component ...');

featuresListPos = get(featureListH, 'position');
thresholdTextPos = textboxPos;
thresholdTextPos(1) = featuresListPos(1);

thresholdText = {'Z-Threshold'};
thresholdTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', thresholdTextPos, 'String', thresholdText);


[thresholdText, newPos] = textwrap(thresholdTextH, thresholdText);

thresholdTextPos(4) = newPos(4);
thresholdTextPos(2) = textboxPos(2) - thresholdTextPos(4);

% Set the wrapped text and the new position
set(thresholdTextH, 'string', thresholdText);
set(thresholdTextH, 'position', thresholdTextPos);

thresholdPos = thresholdTextPos;
thresholdPos(4) = 0.05;
thresholdPos(2) = thresholdPos(2) - yOffset - thresholdPos(4);
thresholdPos(3) = 0.2;
thresholdPos(1) = thresholdTextPos(1) + 0.5*thresholdTextPos(3) - 0.5*thresholdPos(3);
thresholdH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', thresholdPos, 'String', num2str(Z_THRESHOLD_HISTOGRAM), 'tag', 'histogram_threshold');


%%%% Plot display button %%%%
buttonHeight = 0.05;
buttonWidth = 0.2;
% displayPos(3) = buttonWidth;
% displayPos(4) = buttonHeight;

displayPos = [1 - xOffset - buttonWidth, yOffset, buttonWidth, buttonHeight];

doneButtonH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', displayPos, 'string', 'Done', 'tag', 'done_button', 'callback', ...
    {@doneCallback, displayHandle}, 'fontsize', UI_FONT_SIZE - 1);


%%% Set the figure data %%%%
figureData.numGroups = numGroups;
figureData.numFeatures = numFeatures;
figureData.numSubjects = numSubjects;
figureData.dataInfo = dataInfo;
figureData.groupNames = groupNames;
figureData.numComp = numComp;
figureData.backReconstructFiles = backReconstructFiles;
figureData.maskIndices = maskIndices;
figureData.dataLength = dataLength;
figureData.featureDataLength = featureDataLength;
figureData.threshold = threshold;
figureData.fusionFile = fusionFile;
figureData.normalize = normalize;
figureData.featureNormPara = featureNormPara;
figureData.outputDir = outputDir;
figureData.all_comb = all_comb;
figureData.dims = dims;
figureData.voxels = voxels;
figureData.selectedCombNumber = selectedCombNumber;
figureData.selectedFeatures = featureNames;


set(displayHandle, 'userdata', figureData);



%%%%%%%%%%%%%%%%%%%%% Object Callbacks %%%%%%%%%%%%%%%%%%%%%

function doneCallback(hObject, event_data, handles)
% Done callback


try
    
    selGroupsHandle = findobj(handles, 'tag', 'selGroups'); % Groups handle
    selFeaturesHandle = findobj(handles, 'tag', 'selFeatures'); % Features handle
    selCompHandle = findobj(handles, 'tag', 'selComponent'); % Component handle
    thresholdTextH = findobj(handles, 'tag', 'histogram_threshold'); % threshold text
    
    popupH = findobj(handles, 'tag', 'histogram_criteria');
    
    popupString = get(popupH, 'string');
    popupVal = get(popupH, 'value');
    
    histogram_criteria = lower(deblank(popupString(popupVal, :)));
    
    
    % Selected groups features and component
    selGroups = get(selGroupsHandle, 'value');
    selFeatures = get(selFeaturesHandle, 'value');
    selComp = get(selCompHandle, 'value');
    
    if isempty(selGroups)
        error('Groups are not selected');
    end
    
    if isempty(selFeatures)
        error('Features are not selected');
    end
    
    
    % Selected groups
    if length(selGroups) > 2
        disp('Atmost two groups can be selected ...');
        selGroups = selGroups(1:2);
    end
    set(selGroupsHandle, 'value', selGroups);
    
    % Selected features
    if length(selFeatures) > 2
        disp('Atmost two features can be selected ...');
        selFeatures = selFeatures(1:2);
    end
    set(selFeaturesHandle, 'value', selFeatures);
    
    figureData = get(handles, 'userdata');
    
    fusionFile = figureData.fusionFile;
    
    threshold = str2num(get(thresholdTextH, 'string'));
    
    if isempty(threshold)
        error('Please enter a valid number for Z-Threshold');
    end
    
    set(handles, 'pointer', 'watch');
    
    %%% Form histParameters structure
    histParameters.numGroups = figureData.numGroups;
    histParameters.numSubjects = figureData.numSubjects;
    histParameters.selGroupsVal = selGroups;
    histParameters.groupNames = figureData.groupNames;
    histParameters.dataInfo = figureData.dataInfo;
    histParameters.backReconstructFiles = figureData.backReconstructFiles;
    histParameters.mask_ind = figureData.maskIndices;
    histParameters.featureDataLength = figureData.featureDataLength;
    histParameters.normalize = figureData.normalize;
    histParameters.featureNormPara = figureData.featureNormPara;
    histParameters.outputDir = figureData.outputDir;
    histParameters.all_comb = figureData.all_comb;
    histParameters.dims = figureData.dims;
    histParameters.voxels = figureData.voxels;
    histParameters.selectedCombNumber = figureData.selectedCombNumber;
    histParameters.selectedFeatures = figureData.selectedFeatures;
    
    ica_fuse_display_histogram(histParameters, selComp, selFeatures, threshold, histogram_criteria);
    
    set(handles, 'pointer', 'arrow');
    
catch
    set(handles, 'pointer', 'arrow');
    ica_fuse_errorDialog(lasterr, 'Histogram Plot', 'modal');
    ica_fuse_displayErrorMsg;
end
function sortResults = ica_fuse_sortingGUI(fusionFile)
%% Sorting GUI:
% Sorts joint components based on the sorting criteria. Presently sorting
% criteria available are 'Two sample t-test on mixing coefficients between
% groups' 'Divergence between groups'
%
% Inputs:
% 1. fusionFile - Fusion MAT file
%
% Outputs:
% 1. sortResults - Sort results data structure
%

% Load defaults
ica_fuse_defaults;

global FUSION_INFO_MAT_FILE;
global ANATOMICAL_FILE;

% Display Defaults
global IMAGE_VALUES;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global Z_THRESHOLD_HISTOGRAM;
global IMAGES_PER_FIGURE;
global ANATOMICAL_PLANE;
global UI_FONT_NAME;
global UI_FONT_SIZE;


% Select fusion file
if ~exist('fusionFile', 'var')
    fusionFile = ica_fuse_selectEntry('title', 'Select fusion information file for sorting components', ...
        'typeEntity', 'file', 'typeSelection', 'Single', 'filter', ['*', FUSION_INFO_MAT_FILE, '*.mat']);
end

% Load fusion file
load(fusionFile);

% Make sure the selected file is a valid fusion parameter file
if ~exist('fusionInfo', 'var')
    error(['Selected file: ', fusionFile, ' is not a valid fusion parameter file']);
end

[outputDir, fileN, extn] = fileparts(fusionFile);
if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

if ~isfield(fusionInfo, 'run_analysis')
    error('Please run the analysis to view the results');
end

sortResults = struct;

%%%%%%% Get the required parameters from run_analysis field %%%%%%%%%%%
% Output directory
sortParameters.outputDir = outputDir;

% Get the information from the fusion file
numComp = fusionInfo.run_analysis.numComp;

% Number of components
sortParameters.numComp = numComp;

% Output prefix
sortParameters.output_prefix = fusionInfo.run_analysis.prefix;

% Get the feature names
sortParameters.featureNames = str2mat(fusionInfo.run_analysis.dataInfo(1).feature.name);

% Get the groups names
sortParameters.groupNames = str2mat(fusionInfo.run_analysis.dataInfo.name);

% Get the ouput files
sortParameters.outputFiles = fusionInfo.run_analysis.outputFiles;

% Get the dataInfo
sortParameters.dataInfo = fusionInfo.run_analysis.dataInfo;

% Get the PCA files
sortParameters.pcaFiles = fusionInfo.run_analysis.pcaFiles;

% Get the ICA Files
sortParameters.icaFiles = fusionInfo.run_analysis.icaFiles;

% Back reconstruct files
sortParameters.backReconstructFiles = fusionInfo.run_analysis.backReconstructFiles;

% Scale component file
sortParameters.scaleCompFiles = fusionInfo.run_analysis.scaleCompFiles;

% Number of groups
sortParameters.numGroups = fusionInfo.run_analysis.numGroups;

% Number of features
sortParameters.numFeatures = fusionInfo.run_analysis.numFeatures;

% voxel indices
sortParameters.mask_ind = fusionInfo.run_analysis.mask_ind;
[sortParameters.mask_ind] = ica_fuse_form_maskInd(sortParameters.mask_ind, sortParameters.dataInfo);

% feature data length
sortParameters.featureDataLength = fusionInfo.run_analysis.featureDataLength;

% Normalization parameters
sortParameters.featureNormPara = fusionInfo.run_analysis.featureNormPara;

% Number of subjects
numSubjects = fusionInfo.run_analysis.numSubjects;

if length(numSubjects) ~= length(sortParameters.dataInfo)
    % Loop over groups
    numSubjects = repmat(numSubjects, 1, sortParameters.numGroups);
end

sortParameters.numSubjects = numSubjects;

sortParameters.normalize = fusionInfo.run_analysis.normalize;

%% Check this for backward compatibility
if (isfield(fusionInfo.run_analysis, 'all_comb'))
    sortParameters.all_comb = fusionInfo.run_analysis.all_comb;
else
    % Assuming all the features are stacked
    stackAllFeatures = 1;
    optimalFeatures = 0;
    if (length(sortParameters.pcaFiles) > 1)
        optimalFeatures = 1;
    end
    % Get all the combinations
    sortParameters.all_comb = ica_fuse_get_combinations(sortParameters.numFeatures, stackAllFeatures, optimalFeatures);
end

if (~isfield(fusionInfo.run_analysis, 'dims'))
    [dims, voxels] = ica_fuse_getFeatureDIM(sortParameters.mask_ind);
else
    dims = fusionInfo.run_analysis.dims;
    voxels = fusionInfo.run_analysis.voxels;
end

sortParameters.dims = dims;
sortParameters.voxels = voxels;

%%%%%%%%%%%%% Form prompt string for sorting criteria %%%%%%%%
% Check the necessary M files
% check_ttest2 = which('ttest2.m');
% check_kstest2 = which('kstest2.m');
%
% if isempty(check_ttest2)
%     error('Need statistics toolbox to sort components');
% end

%%%%% Need atleast two groups to sort the components %%%%%%
if sortParameters.numGroups < 2

    error('Need atleast two groups to sort components');

end
%%%%%%% end for error checking regarding number of groups %%%%

% Get divergence name and value from defaults
[div_name, div_value] = ica_fuse_getDivergencePara;

if ~isempty(div_value)
    divStr = ['(', upper(div_name), '(', num2str(div_value), '))'];
else
    divStr = ['(', upper(div_name), ')'];
end

% Draw a figure with two listboxes, a check box, Convert to z scores,
% Threshold, Slices, Anatomical Plane
displayTag = 'SortingGUI_Fusion';

%%%%%% Delete any previous figures of display GUI %%%%%
displayH = findobj('tag', displayTag);

for ii = 1:length(displayH)
    delete(displayH(ii));
end
%%%%%% end for deleting previous figures of display GUI %%%%

[displayHandle] = ica_fuse_getGraphics('Sorting GUI', 'normal', displayTag, 'off');
set(displayHandle, 'menubar', 'none');


numPara = 1;
sortingCriteria(numPara).displayString = 'ttest2 on mixing coeff';
sortingCriteria(numPara).string = 'ttest2';

numPara = numPara + 1;
sortingCriteria(numPara).displayString =  ['Spatial div ', divStr];
sortingCriteria(numPara).string = 'divergence';

clear numPara;


choiceString = str2mat(sortingCriteria.displayString);


% Plot popups for sorting criteria and histogram type
% plot listbox for selecting groups
% Z-threshold on component maps

% offsets
xOffset = 0.04; yOffset = 0.04;
popupTextHeight = 0.05; popupTextWidth = 1 - 2*xOffset; yPos = 0.95;

popupTextWidth = 0.45;
popupTextPos = [xOffset, yPos - yOffset - popupTextHeight, popupTextWidth, popupTextHeight];

%%% Sorting criteria text %%%%

% Plot Text
sortTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', popupTextPos, 'String', 'Select sorting criteria', 'fontsize', UI_FONT_SIZE - 1, ...
    'horizontalalignment', 'center');

[sortTextH] = ica_fuse_wrapStaticText(sortTextH);

popupTextPos = get(sortTextH, 'position');

popupWidth = 0.4;
popupPos = popupTextPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = popupWidth;

% Plot popup
sortPopupH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupPos, 'String', choiceString, 'fontsize', UI_FONT_SIZE - 1, 'userdata', sortingCriteria, ...
    'tag', 'sorting_criteria', 'callback', {@sortPopupCallback, displayHandle});

%%%% End for plotting sorting criteria and text %%%%%


%%%%%% Histogram Criteria %%%%%%

popupTextPos = [xOffset, popupTextPos(2) - yOffset - popupTextPos(4), popupTextWidth, popupTextHeight];

% Plot Text
histTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', popupTextPos, 'String', 'Select type of histogram', 'fontsize', UI_FONT_SIZE - 1, ...
    'horizontalalignment', 'center');

[histTextH] = ica_fuse_wrapStaticText(histTextH);

popupTextPos = get(histTextH, 'position');

popupPos = popupTextPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = popupWidth;

choiceString = str2mat('Feature', 'Component');

% Plot popup
histPopupH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupPos, 'String', choiceString, 'fontsize', UI_FONT_SIZE - 1, 'tag', 'histogram_criteria', ...
    'callback', {@histPopupCallback, displayHandle});

%%%%%% End for plotting histogram criteria %%%%%%


editTextHeight = 0.05;
editTextWidth = 0.45;
editTextPos = popupTextPos;
editTextPos = [xOffset, editTextPos(2) - yOffset - editTextPos(4), editTextWidth, editTextHeight];

%% Select z-threshold for feature histograms
zTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', editTextPos, 'String', 'Select Z-threshold', 'fontsize', UI_FONT_SIZE - 1, ...
    'horizontalalignment', 'center');

[zTextH] = ica_fuse_wrapStaticText(zTextH);

editTextPos = get(zTextH, 'position');

editPos = editTextPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.2;

% Plot text
zEditH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', num2str(Z_THRESHOLD_HISTOGRAM), 'fontsize', UI_FONT_SIZE - 1, 'tag', ...
    'threshold');

listboxWidth = 0.4;

%%%%%%%%%%%%%% Groups Listbox %%%%%%%%%%%%%%%%%%%%

% Plot textbox
textBoxWidth = listboxWidth; textboxHeight = 0.05;
textboxPos = [xOffset, editPos(2) - editPos(4) - yOffset, textBoxWidth, textboxHeight];

% Plot groups text
groupsTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Groups', 'fontsize', UI_FONT_SIZE - 1);

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
listboxPos(4) = 0.16;
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);

if sortParameters.numGroups >= 2
    groupsListVal = [1:2];
else
    groupsListVal = 1;
end

% Plot groups listbox
groupsListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', sortParameters.groupNames, 'value', ...
    groupsListVal, 'min', 0, 'max', 2, 'tag', 'selGroups', 'fontsize', UI_FONT_SIZE - 1, 'TooltipString', ...
    'Select two groups ...');

%%%%%%%%%%%%%%% End for plotting groups listbox and text %%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% Features Text and Listbox %%%%%%%%%%%%%%%%%%%%%%
textboxPos(1) = 1 - xOffset - listboxWidth;
textboxPos(3) = listboxWidth;

% Plot Feature text
featureTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Feature');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);

toolTipStr = 'Select features ...';

% Plot feature listbox
featureListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', sortParameters.featureNames, 'value', ...
    [1:size(sortParameters.featureNames, 1)], 'min', 0, 'max', 2, 'tag', 'selFeature', 'fontsize', UI_FONT_SIZE - 1, ...
    'TooltipString', toolTipStr);

% Textbox width
clear textboxPos;

% plot display button
buttonWidth = 0.12; buttonHeight = 0.05;
donePos = [0.5 - 0.5*buttonWidth, listboxPos(2) - 2*yOffset - buttonHeight, buttonWidth, buttonHeight];
doneButtonH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', donePos, 'string', 'Done', 'tag', 'display_button', 'callback', ...
    {@doneCallback, displayHandle}, 'fontsize', UI_FONT_SIZE - 1);

set(displayHandle, 'userdata', sortParameters);

sortPopupH = findobj(displayHandle, 'tag', 'sorting_criteria');
sortPopupCallback(sortPopupH, [], displayHandle);

% Make the graphics visible after plotting all controls
set(displayHandle, 'visible', 'on');

figure(displayHandle);

try
    set(displayHandle, 'visible', 'on');
    waitfor(displayHandle);
catch
end


if isappdata(0, 'sortCompData')
    sortResults = getappdata(0, 'sortCompData');
    rmappdata(0, 'sortCompData');
else
    error('Sorting GUI figure was quit');
end

%%%%%%%%%%%%%%%%%%%%%%% End for plotting features listbox and text %%%%%%%

function displayDefaultsCallback(hObject, event_data, handles)

try

    controlPara = get(hObject, 'userdata');
    sortParameters = get(handles, 'userdata');

    [controlPara, sortParameters] = ica_fuse_get_display_parameters(sortParameters, controlPara, 'on');

    set(hObject, 'userdata', controlPara);
    set(handles, 'userdata', sortParameters);

catch

    disp(lasterr);

end


function doneCallback(hObject, event_data, handles)

try

    set(handles, 'pointer', 'watch');

    % Component listbox
    %groupListH = findobj(handles, 'tag', 'selGroups');

    sortParameters = get(handles, 'userdata');

    % Number of groups and features
    %numGroups = sortParameters.numGroups;

    %numFeatures = sortParameters.numFeatures;

    % Number of components
    %numComp = sortParameters.numComp;

    % Output directory
    outputDir = sortParameters.outputDir;

    cd(outputDir);

    % Output prefix
    %output_prefix = sortParameters.output_prefix;

    % data info
    %dataInfo = sortParameters.dataInfo;
    % Output files
    %outputFiles = sortParameters.outputFiles;

    % Pca files
    %pcaFiles = sortParameters.pcaFiles;

    % ica files
    %icaFiles = sortParameters.icaFiles;

    % Back reconstruct files
    %back_reconstruct_files = sortParameters.backReconstructFiles;

    % Voxels
    %mask_ind = sortParameters.mask_ind;

    %sortingCriteria = sortParameters.sortingCriteria;

    % Get the display fields

    [selGroupNames, selGroupsVal] = ica_fuse_get_value_uicontrol(handles, 'selGroups');
    if length(selGroupsVal) > 2
        error('Number of groups that can be selected is 2');
    end

    if length(selGroupsVal) == 1
        error('Cannot sort the components as there is only one group');
    end
    selGroupNames = selGroupNames(selGroupsVal, :);
    sortParameters.selGroupNames = selGroupNames;
    sortParameters.selGroupsVal = selGroupsVal;

    selFeatureH = findobj(handles, 'tag', 'selFeature');

    % Selected Features
    [selectedFeature, selectedFeatureVal] = ica_fuse_get_value_uicontrol(handles, 'selFeature');


    % Type of histogram
    histogramH = findobj(handles, 'tag', 'histogram_criteria');
    histPopVal = get(histogramH, 'value');
    histPopStr = get(histogramH, 'string');
    histogramCriteria = lower(deblank(histPopStr(histPopVal, :)));

    % Threshold
    thresholdTextH = findobj(handles, 'tag', 'threshold');
    threshold = str2num(get(thresholdTextH, 'string'));

    % Check the selected threshold
    if strcmpi(histogramCriteria, 'feature')
        if isempty(threshold)
            error('Please enter a valid number for threshold');
        end
    end

    selectedFeature = deblank(selectedFeature(selectedFeatureVal, :));
    sortParameters.selectedFeature = selectedFeature;
    sortParameters.selectedFeatureVal = selectedFeatureVal;
    sortParameters.z_threshold = threshold;
    sortParameters.histogramCriteria = histogramCriteria;

    delete(handles);

    % Sort components
    sortResults = ica_fuse_sort_components(sortParameters);


    % Store sorting results in application data
    setappdata(0, 'sortCompData', sortResults);

catch

    delete(handles);

    ica_fuse_displayErrorMsg;

end

function sortPopupCallback(hObject, event_data, handles)
% Sort popup callback
% hObject: Popup object
% handles: Sorting GUI figure

sortingCriteria = '';

getString = get(hObject, 'string');
getValue = get(hObject, 'value');
objUserdata = get(hObject, 'userdata');

sortParameters = get(handles, 'userdata');

sortingCriteria = objUserdata(getValue).string;
sortParameters.sortingCriteria = sortingCriteria;

selFeatureH = findobj(handles, 'tag', 'selFeature');

selFeatureString = get(selFeatureH, 'string');
selFeatureVal = get(selFeatureH, 'value');

selGroupH = findobj(handles, 'tag', 'selGroups');

selGroupString = get(selGroupH, 'string');
selGroupVal = get(selGroupH, 'value');


histogramH = findobj(handles, 'tag', 'histogram_criteria');
z_thresholdH = findobj(handles, 'tag', 'threshold');

if strcmpi(sortingCriteria, 'ttest2')
    set(selFeatureH, 'value', (1:size(selFeatureString, 1)));
    set(selFeatureH, 'enable', 'off');
    set(histogramH, 'enable', 'off');
    set(z_thresholdH, 'enable', 'off');
else
    set(selFeatureH, 'enable', 'on');
    set(histogramH, 'enable', 'on');
    set(z_thresholdH, 'enable', 'on');
end

if strcmpi(get(histogramH, 'enable'), 'on')
    histPopupCallback(histogramH, [], handles);
end

set(handles, 'userdata', sortParameters);


function histPopupCallback(hObject, event_data, handles)
% Enable z-threshold text box only if histogram of features is selected


z_thresholdH = findobj(handles, 'tag', 'threshold');

% Get the selected option
getVal = get(hObject, 'value');
getStr = get(hObject, 'string');
selectedStr = deblank(getStr(getVal, :));

% Enable threshold text only for features
if strcmpi(selectedStr, 'feature')
    set(z_thresholdH, 'enable', 'on');
else
    set(z_thresholdH, 'enable', 'off');
end

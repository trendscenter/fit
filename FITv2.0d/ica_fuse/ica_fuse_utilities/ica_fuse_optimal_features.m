function ica_fuse_optimal_features(fusionFile, selGroups, sortingCriteria)
%% Function to get optimal features.
% Loads back reconstruction for the respective combinations and displays
% a bar graph ranking the features.
%
% Inputs:
% fusionFile - full file path


% if isempty(which('ttest2.m'))
%     error('Need statistics toolbox to rank features');
% end

% Load fusion defaults
ica_fuse_defaults;

% Declare variables
global FUSION_INFO_MAT_FILE;
global FIG_FG_COLOR;
global FIG_FG_COLOR;
global UI_FONT_NAME;
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size
global Z_THRESHOLD_HISTOGRAM;


if ~exist('fusionFile', 'var')
    fusionFile = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        ['*', FUSION_INFO_MAT_FILE, '*.mat'], 'title', 'Select joint ICA fusion information file for ranking features');
    drawnow;
end

outputDir = fileparts(fusionFile);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

load(fusionFile);

if ~exist('fusionInfo', 'var')
    error(['Selected file: ', fusionFile, ' is not a valid fusion parameter file']);
end

if ~isfield(fusionInfo, 'run_analysis')
    error('Please run the analysis inorder to rank features');
end

[div_name, div_value] = ica_fuse_getDivergencePara;

%% Get the necessary parameters
figureData.outputDir = outputDir;
figureData.dataInfo = fusionInfo.run_analysis.dataInfo;
figureData.numComp = fusionInfo.run_analysis.numComp;
figureData.output_prefix = fusionInfo.run_analysis.prefix;
figureData.featureNames = str2mat(fusionInfo.run_analysis.dataInfo(1).feature.name);
figureData.groupNames = str2mat(fusionInfo.run_analysis.dataInfo.name);
figureData.pcaFiles = fusionInfo.run_analysis.pcaFiles;
figureData.icaFiles = fusionInfo.run_analysis.icaFiles;
figureData.backReconstructFiles = fusionInfo.run_analysis.backReconstructFiles;
figureData.numGroups = fusionInfo.run_analysis.numGroups;
figureData.numFeatures = fusionInfo.run_analysis.numFeatures;
figureData.mask_ind = fusionInfo.run_analysis.mask_ind;
[figureData.mask_ind] = ica_fuse_form_maskInd(figureData.mask_ind, figureData.dataInfo);
figureData.featureDataLength = fusionInfo.run_analysis.featureDataLength;
figureData.numGroups = fusionInfo.run_analysis.numGroups;
figureData.numFeatures = fusionInfo.run_analysis.numFeatures;
figureData.normalize = fusionInfo.run_analysis.normalize;
figureData.featureNormPara = fusionInfo.run_analysis.featureNormPara;
figureData.div_name = div_name;
figureData.div_value = div_value;

%% Check this for backward compatibility
if (isfield(fusionInfo.run_analysis, 'all_comb'))
    figureData.all_comb = fusionInfo.run_analysis.all_comb;
else
    % Assuming all the features are stacked
    stackAllFeatures = 1;
    optimalFeatures = 0;
    if (length(figureData.pcaFiles) > 1)
        optimalFeatures = 1;
    end
    % Get all the combinations
    figureData.all_comb = ica_fuse_get_combinations(figureData.numFeatures, stackAllFeatures, optimalFeatures);
end

if (~isfield(fusionInfo.run_analysis, 'dims'))
    [dims, voxels] = ica_fuse_getFeatureDIM(figureData.mask_ind);
else
    dims = fusionInfo.run_analysis.dims;
    voxels = fusionInfo.run_analysis.voxels;
end

figureData.dims = dims;
figureData.voxels = voxels;

%dataInfo = figureData.dataInfo;
numGroups = figureData.numGroups;
numFeatures = figureData.numFeatures;
back_reconstruct_files = figureData.backReconstructFiles;
groupNames = str2mat(figureData.dataInfo.name);

if length(fusionInfo.run_analysis.numSubjects) ~= length(figureData.dataInfo)
    numSubjects = fusionInfo.run_analysis.numSubjects;
    numSubjects = repmat(numSubjects, 1, length(figureData.dataInfo));
else
    numSubjects = fusionInfo.run_analysis.numSubjects;
end

clear fusionInfo;

figureData.numSubjects = numSubjects;

if numGroups == 1
    error('Need atleast two groups to do optimization');
end

if (length(back_reconstruct_files) == 1) && (numFeatures > 1)
    error('Error:OptimalFeatures', ['Cannot run optimal features as only joint ICA is run. ', ...
        ' \nYou need to set OPTIMIZE_FEATURES variable in ica_fuse_defaults.m as ''yes'' \nand run analysis again.']);
end

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
%%%%%% end for deleting previous figures  %%%%

[displayHandle] = ica_fuse_getGraphics('Optimal Features', 'normal', displayTag, 'on');
set(displayHandle, 'menubar', 'none');

fitHelpTitle = uimenu('parent', displayHandle, 'label', 'FIT-Help');
fitHelpMenu = uimenu(fitHelpTitle, 'label', 'Optimal Features', 'callback', ...
    'ica_fuse_openHTMLHelpFile(''fit_optimal_features.htm'');');

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
    'tag', 'sorting_criteria');

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
    'position', popupPos, 'String', choiceString, 'fontsize', UI_FONT_SIZE - 1, 'tag', 'histogram_type', ...
    'callback', {@histPopupCallback, displayHandle});

%%%%%% End for plotting histogram criteria %%%%%%


listboxWidth = 0.4;

%%%%%%%%%%%%%% Groups Listbox %%%%%%%%%%%%%%%%%%%%

% Plot textbox
textBoxWidth = listboxWidth; textboxHeight = 0.05;
textboxPos = [xOffset, popupTextPos(2) - popupTextPos(4) - 2*yOffset, textBoxWidth, textboxHeight];

% Plot groups text
groupsTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Select Groups', 'fontsize', UI_FONT_SIZE - 1, ...
    'horizontalalignment', 'center');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
listboxPos(4) = 0.3;
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


%%%%%% Z-Threshold %%%%%%
textboxPos = textboxPos;
textboxPos(2) = listboxPos(2) + 0.5*listboxPos(4) + 0.5*textboxPos(4);
textboxPos(1) = popupPos(1);
textboxPos(3) = 0.3;

% Plot Text
zTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Z-Threshold', 'fontsize', UI_FONT_SIZE - 1, ...
    'horizontalalignment', 'center');

editPos = textboxPos;

editWidth = 0.2;
editHeight = 0.05;
editPos(1) = editPos(1) + 0.5*editPos(3) - 0.5*editWidth;
editPos(2) = editPos(2) - editPos(4) - 0.5*editHeight;
editPos(3) = editWidth;
editPos(4) = editHeight;

% Plot edit
thresholdH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', num2str(Z_THRESHOLD_HISTOGRAM), 'fontsize', UI_FONT_SIZE - 1, 'tag', ...
    'threshold');

%%%%%% End for plotting Z-Threshold %%%%%%


groupListPos = get(groupsListH, 'position');

buttonHeight = 0.05;
buttonWidth = 0.2;

donePos = [0.75 - 0.5*buttonWidth, yOffset + 0.5*buttonHeight, buttonWidth, buttonHeight];

doneButtonH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', donePos, 'string', 'Calculate', 'tag', 'done_button', 'callback', ...
    {@doneCallback, displayHandle}, 'fontsize', UI_FONT_SIZE - 1);

cancelPos = [0.25 - 0.5*buttonWidth, yOffset + 0.5*buttonHeight, buttonWidth, buttonHeight];

cancelButtonH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', cancelPos, 'string', 'Close', 'tag', 'cancel', 'callback', ...
    {@cancelCallback, displayHandle}, 'fontsize', UI_FONT_SIZE - 1);

set(displayHandle, 'userdata', figureData);



%%%%% Done Callback %%%%%
function doneCallback(hObject, event_data, handles)

%global OPTIMIZATION_DATA;
ica_fuse_defaults;
global FIG_FG_COLOR;
global FIG_FG_COLOR;
global UI_FONT_NAME;
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size
global NUM_BAR_ELEMENTS_PER_AXES;

try
    set(handles, 'pointer', 'watch');
    % Selected groups
    selGroupsHandle = findobj(handles, 'tag', 'selGroups');
    selGroups = get(selGroupsHandle, 'value');
    groupNames = get(selGroupsHandle, 'string');

    if (length(selGroups) < 2)
        error('Need atleast two groups to sort the components');
    else
        selGroups = selGroups(1:2);
    end

    % Sorting criteria
    sortingCriteriaH = findobj(handles, 'tag', 'sorting_criteria');
    sortingVal = get(sortingCriteriaH, 'value');

    sortingData = get(sortingCriteriaH, 'userdata');
    sortingType = deblank(sortingData(sortingVal).string);


    % Type of histogram
    histogramH = findobj(handles, 'tag', 'histogram_type');
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

    figureData = get(handles, 'userdata');

    figureData.histogramCriteria = histogramCriteria;
    figureData.sortingCriteria = sortingType;
    figureData.selGroupsVal = selGroups;
    figureData.groupNames  = groupNames;
    figureData.z_threshold = threshold;

    div_name = figureData.div_name;   

    % Rank features
    [optimization_data, legendString, sorted_ind] = ica_fuse_rank_features(figureData);

    % Output file
    outFile = fullfile(figureData.outputDir, [figureData.output_prefix, '_optimal_feature_results.mat']);

    %     XTickStr = cell(1, length(sorted_ind));
    %     for nn = 1:length(sorted_ind)
    %         XTickStr{nn} = num2str(sorted_ind(nn));
    %     end

    XTickStr = cellstr(num2str(sorted_ind(:)));

    yLabelStr = [upper(div_name), ' divergence'];


    % Plot bar graph
    titleFig = 'Optimal Feature Results';
    selectedgroupNames = cellstr(deblank(groupNames(selGroups, :)));
    axesTitle = ['Spatial divergence between groups ', selectedgroupNames{1}, ' & ', selectedgroupNames{2}];

    if (length(optimization_data) < NUM_BAR_ELEMENTS_PER_AXES)
        %% Single plot
        figHandle = ica_fuse_getGraphics(titleFig, 'normal', titleFig, 'on');
        axes_h = axes('parent', figHandle, 'units', 'normalized', 'position', [0.12 0.12 0.8 0.8]);
        bar_h = bar([optimization_data.divergence]');
        ylabel(yLabelStr, 'parent', axes_h);
        xlabel('Features', 'parent', axes_h);
        title(axesTitle, 'parent', axes_h);
        set(axes_h, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
        set(axes_h, 'XTickLabel', XTickStr);
        set(axes_h, 'fontname', UI_FONT_NAME, 'fontunits', UI_FONT_UNITS, 'fontsize', UI_FONT_SIZE);
        figure(figHandle);
        ica_fuse_legend(str2mat(legendString));

    else
        %% Multi bar plot
        ica_fuse_multi_bar_plot([optimization_data.divergence]', 'XTickLabel', XTickStr, 'xlabel', 'Features', 'ylabel', yLabelStr, ...
            'title', titleFig, 'axesTitle', axesTitle, 'legend', legendString);

    end

    assignin('base', 'OPTIMIZATION_DATA', optimization_data);
    ica_fuse_save(outFile, 'optimization_data', 'selectedgroupNames');

    fprintf('Divergence values are stored in workspace variable OPTIMIZATION_DATA and\nalso saved in file %s', outFile);
    fprintf('\n');
    set(handles, 'pointer', 'arrow');

catch
    set(handles, 'pointer', 'arrow');
    ica_fuse_errorDialog(lasterr, 'Optimal Features Error', 'modal');
    ica_fuse_displayErrorMsg;
    delete(handles);
end

function cancelCallback(hObject, event_data, handles)

delete(handles);


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
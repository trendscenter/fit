function ica_fuse_create_movie_erp_fmri(fusionFile, dispParameters)
% Create ERP-fMRI movie using the information from the fusion information
% file
%
% Inputs:
% 1. FusionFile - full file path
% 2. dispParameters - Display parameters structure
% For example:
% a. dispParameters.convert_to_z = 1; % converts to z-scores
% b. dispParameters.z_threshold = 1.0; % threshold
% c. Image Values: 1 - Positive and negative, 2 - Positive, 3 - Abs, 4 -
% Negative
% dispParameters.image_values = 1;
% d. Anatomical View: Axial, sagittal, coronal
% dispParameters.anatomical_plane = 'axial';

% Load defaults
ica_fuse_defaults;

global FUSION_INFO_MAT_FILE; % Fusion information file
global JOINT_COMPONENT_NAMING; % JOINT COMPONENT IMAGE FILE NAMING

% Figure color
global FIG_BG_COLOR; % Figure background color
global FIG_FG_COLOR; % Figure foreground color

% FONT DEFAULTS
global UI_FONT_NAME; % font name
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size

if ~exist('fusionFile', 'var')
    fusionFile = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        ['*', FUSION_INFO_MAT_FILE, '*.mat'], 'title', 'Select fusion information file ...');
    drawnow;
end

load(fusionFile);

if ~exist('fusionInfo', 'var')
    error(['Selected file: ', fusionFile, ' is not a valid fusion info file']);
end

[outputDir, fileName, extn] = fileparts(fusionFile);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

% Remove the field run_analysis from fusionInfo variable
if ~isfield(fusionInfo, 'run_analysis')
    if ~isfield(fusionInfo.run_analysis, 'outputFiles')
        error(['Please run the analysis in order to view the ERP-fMRI movie']);
    end
end


% Output files prefix
output_prefix = fusionInfo.run_analysis.prefix;

% Number of groups
numGroups = fusionInfo.run_analysis.numGroups;

% Number of features
numFeatures = fusionInfo.run_analysis.numFeatures;

% Data Info
dataInfo = fusionInfo.run_analysis.dataInfo;

% Output files
outputFiles = fusionInfo.run_analysis.outputFiles;

% Modalities
modalities = lower(str2mat(dataInfo(1).feature.modality));

% Feature Names
featureNames = str2mat(dataInfo(1).feature.name);

% Group names
groupNames = str2mat(dataInfo.name);

% Number of components
numComp = fusionInfo.run_analysis.numComp;

% Mask indices
mask_ind = fusionInfo.run_analysis.mask_ind;

[dims, voxels] = ica_fuse_getFeatureDIM(mask_ind);

% Number of subjects
numSubjects = fusionInfo.run_analysis.numSubjects;

checkEEG = strmatch('eeg', modalities, 'exact');
checkfMRI = strmatch('fmri', modalities, 'exact');

if isempty(checkEEG)
    error(['EEG modality doesn''t exist in order to create ERP-fMRI movie']);
end

if isempty(checkfMRI)
    error(['fMRI modality doesn''t exist in order to create ERP-fMRI movie']);
end

%%%%%%%%%%% GUI for getting selected fMRI and EEG features%%%%%%%%

displayTag = 'erp_fmri_movie_GUI';


%%%%%% Delete any previous figures%%%%%
displayH = findobj('tag', displayTag);

for ii = 1:length(displayH)
    delete(displayH(ii));
end
%%%%%% end for deleting previous figures of display GUI %%%%

[displayHandle] = ica_fuse_getGraphics('ERP-fMRI movie GUI', 'normal', displayTag, 'off');
set(displayHandle, 'menubar', 'none');


fitHelpTitle = uimenu('parent', displayHandle, 'label', 'FIT-Help');
fitHelpMenu = uimenu(fitHelpTitle, 'label', 'Create ERP-fMRI Movie', 'callback', ...
    'ica_fuse_openHTMLHelpFile(''fit_fmri_eeg_fusion.htm'');');


% offsets
xOffset = 0.04; yOffset = 0.035;

outputFileTextHeight = 0.05; outputFileTextWidth = 0.52;

yPos = 0.95;
outputFilePos = [xOffset, yPos - yOffset - outputFileTextHeight, outputFileTextWidth, outputFileTextHeight];

outputTextString = {'Enter output file name to save movie: '};
% Plot static text
outputFileH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', outputFilePos, 'String', outputTextString, 'fontsize', UI_FONT_SIZE - 1);

[outputTextString, newPos] = textwrap(outputFileH, outputTextString);

outputFilePos(4) = newPos(4);
outputFilePos(2) = yPos - yOffset - outputFilePos(4);

% Set the wrapped text and the new position
set(outputFileH, 'string', outputTextString);
set(outputFileH, 'position', outputFilePos);


%%%% Plot edit control
editPos = outputFilePos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.32;
editPos(4) = 0.05;
editPos(2) = editPos(2) + 0.5*outputFilePos(4) - 0.5*editPos(4);

editString = 'erp_fmri';

editH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', editString, 'fontsize', UI_FONT_SIZE - 1, 'tag', ...
    'output_file');

listboxWidth = 0.35;

%%%%%%%%%%%%%% fmri Listbox %%%%%%%%%%%%%%%%%%%%

% Plot textbox
textBoxWidth = listboxWidth; textboxHeight = 0.05;
textboxPos = [xOffset, editPos(2) - editPos(4) - 2*yOffset, textBoxWidth, textboxHeight];

% Plot fmri text
fmriTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'fMRI', 'fontsize', UI_FONT_SIZE - 1);

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
listboxPos(4) = 0.17;
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);

% Plot fmri listbox
fmriListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', str2mat(featureNames(checkfMRI, :)), 'value', ...
    1, 'min', 0, 'max', 1, 'tag', 'selfMRI', 'fontsize', UI_FONT_SIZE - 1, 'TooltipString', ...
    'Select fmri feature ...');

%%%%%%%%%%%%%%% End for plotting groups listbox and text %%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% EEG Listbox %%%%%%%%%%%%%%%%%%%%%%
textboxPos(1) = 1 - textboxPos(3) - 2*xOffset;
textboxPos(3) = listboxWidth;

% Plot Feature text
eegTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'EEG');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);

% Plot eeg listbox
eegListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', str2mat(featureNames(checkEEG, :)), 'value', ...
    1, 'min', 0, 'max', 1, 'tag', 'selEEG', 'fontsize', ...
    UI_FONT_SIZE - 1, 'TooltipString', 'Select eeg feature ...');



%%%%% Plot Groups listbox %%%%%
fmriListPos = get(fmriListH, 'position');


textboxPos(2) = fmriListPos(2) - 2*yOffset;
textboxPos(1) = fmriListPos(1);

% Plot groups text
groupTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Select a group');


listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);
listboxPos(1) = textboxPos(1);

% Plot groups listbox
groupListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', groupNames, 'value', 1, 'min', 0, 'max', 1, 'tag', ...
    'selGroups', 'fontsize', UI_FONT_SIZE - 1, 'TooltipString', 'Select a group ...');


%%%%%%%%%%%%%%% Plot Step Size %%%%%%%%%%%%%%%
% Plot text for step size
tboxPos(1) = listboxPos(1);
tboxPos(2) = listboxPos(2) - 0.05 - yOffset;
tboxPos(3) = 0.25;
tboxPos(4) = 0.05;

tbH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', tboxPos, 'String', 'Select Step Size:', 'fontsize', UI_FONT_SIZE - 1);
tbH = ica_fuse_wrapStaticText(tbH);
tboxPos = get(tbH, 'position');

editWidth = 0.2;
tboxPos(1) = tboxPos(1) + 0.5*tboxPos(3) - 0.5*editWidth;
tboxPos(3) = editWidth;
tboxPos(2) = tboxPos(2) - 0.5*tboxPos(4) - yOffset;


if voxels < 10
    step = 1;
elseif voxels > 10 && voxels < 100
    step = 10;
elseif voxels > 100 && voxels < 1000
    step = 100;
else
    step = 500;
end

tbEditH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', tboxPos, 'String', num2str(step), 'tag', 'step_size', 'fontsize', UI_FONT_SIZE - 1);

%%%%%%%% End for plotting step size %%%%%%%%%%%%

tboxPos2 = get(tbH, 'position');
tboxPos2(1) = tboxPos2(1) + tboxPos2(3) + xOffset; %1 - textboxPos(3) - 2*xOffset;
tb2H = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', tboxPos2, 'String', 'Select Time Units:', 'fontsize', UI_FONT_SIZE - 1);
tb2H = ica_fuse_wrapStaticText(tb2H);

tboxPos2(1) = tboxPos2(1) + 0.5*tboxPos2(3) - 0.5*editWidth;
tboxPos2(2) = tboxPos2(2) - 0.5*tboxPos2(4) - yOffset;
tboxPos2(3) = editWidth;
ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', tboxPos2, 'String', {'Millisecs', 'Secs'}, 'tag', 'time_units', 'fontsize', UI_FONT_SIZE - 1);

tboxPos3 = get(tb2H, 'position');
tboxPos3(1) = tboxPos3(1) + tboxPos3(3) + xOffset; %1 - textboxPos(3) - 2*xOffset;
tb3H = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', tboxPos3, 'String', 'Time Range (millisecs):', 'fontsize', UI_FONT_SIZE - 1);
tb3H = ica_fuse_wrapStaticText(tb3H);

tboxPos3(1) = tboxPos3(1) + 0.5*tboxPos3(3) - 0.5*editWidth;
tboxPos3(2) = tboxPos3(2) - 0.5*tboxPos3(4) - yOffset;
tboxPos3(3) = editWidth;
ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', tboxPos3, 'String', '', 'tag', 'time_range', 'fontsize', UI_FONT_SIZE - 1);


listboxPos(1) = 1 - textboxPos(3) - 2*xOffset;
compTextPos = textboxPos;
compTextPos(1) = listboxPos(1);
% Plot components text
compTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', compTextPos, 'String', 'Select components');


listboxPos(2) = compTextPos(2) - 0.5*compTextPos(4) - listboxPos(4);
listboxPos(1) = compTextPos(1);

for nComp = 1:numComp
    compStr(nComp).name = num2str(nComp);
end

% Plot groups listbox
compListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', str2mat(compStr.name), 'value', 1, 'min', 0, 'max', 2, 'tag', ...
    'selComp', 'fontsize', UI_FONT_SIZE - 1, 'TooltipString', 'Select components ...');



%%%% Plot display button %%%%
buttonHeight = 0.05;
buttonWidth = 0.2;

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
figureData.featureNames = featureNames;
figureData.numComp = numComp;
figureData.mask_ind = mask_ind;
figureData.fusionFile = fusionFile;
figureData.outputDir = outputDir;
figureData.outputFiles = outputFiles;
figureData.checkfMRI = checkfMRI;
figureData.checkEEG = checkEEG;
figureData.dispParameters = dispParameters;

set(displayHandle, 'userdata', figureData, 'visible', 'on');


function doneCallback(hObject, event_data, handles)

ica_fuse_defaults;

% Anatomical file
global ANATOMICAL_FILE;

% Figure colors
global FIG_BG_COLOR;
global FIG_FG_COLOR;

% FONT DEFAULTS
global UI_FONT_NAME; % font name
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size

% Figure data
figureData = get(handles, 'userdata');

try
    
    outputFiles = figureData.outputFiles;
    numGroups = figureData.numGroups;
    numFeatures = figureData.numFeatures;
    numSubjects = figureData.numSubjects;
    dataInfo = figureData.dataInfo;
    groupNames = figureData.groupNames;
    featureNames = figureData.featureNames;
    numComp = figureData.numComp;
    mask_ind = figureData.mask_ind;
    fusionFile = figureData.fusionFile;
    outputDir = figureData.outputDir;
    checkfMRI = figureData.checkfMRI;
    checkEEG = figureData.checkEEG;
    dispParameters = figureData.dispParameters;
    
    [dims, voxels] = ica_fuse_getFeatureDIM(mask_ind);
    
    % Get step size
    stepH = findobj(handles, 'tag', 'step_size');
    
    try
        step = str2num(get(stepH, 'string'));
    catch
        error('Enter a valid number for step size');
    end
    
    step = round(step);
    
    if step > voxels
        error('Error:StepSize', 'Step size (%s) cannot be greater than the number of voxels (%s)', num2str(step), num2str(voxels));
    end
    
    if step == 0
        error('Please enter a valid step size.');
    end
    
    % End for getting step size
    
    movieFileH = findobj(handles, 'tag', 'output_file');
    movieFileName = get(movieFileH, 'string');
    
    [status, message] = ica_fuse_errorCheck(movieFileName, 'output_prefix');
    if status == 0
        error(message);
    end
    
    movieFileName = fullfile(outputDir, [movieFileName, '.avi']);
    
    
    selfMRIH = findobj(handles, 'tag', 'selfMRI');
    selEEGH = findobj(handles, 'tag', 'selEEG');
    selGroupsH = findobj(handles, 'tag', 'selGroups');
    selCompH = findobj(handles, 'tag', 'selComp');
    
    % Get the values from the handles
    selfMRI = get(selfMRIH, 'value');
    selEEG = get(selEEGH, 'value');
    selGroups = get(selGroupsH, 'value');
    selComp = get(selCompH, 'value');
    
    timeH = findobj(handles, 'tag', 'time_units');
    timeVal = get(timeH, 'value');
    timeStr = cellstr(get(timeH, 'string'));
    
    timeRange = get(findobj(handles, 'tag', 'time_range'), 'string');
    
    disp('Using information from display defaults in display GUI to threshold the spatial map ...');
    
    % Display parameters
    convertToZ = dispParameters.convert_to_z;
    thresh = abs(dispParameters.z_threshold);
    imageValues = dispParameters.image_values;
    anatomicalView = lower(dispParameters.anatomical_plane);
    slices_in_mm = dispParameters.slices_in_mm;
    
    if (min(thresh) == max(thresh))
        thresh = thresh(1);
    end
    
    % Selected fMRI and EEG index
    fmriIndex = checkfMRI(selfMRI);
    erpIndex = checkEEG(selEEG);
    
    fmriFeatureName = deblank(featureNames(fmriIndex, :));
    erpFeatureName = deblank(featureNames(erpIndex, :));
    
    
    delete(handles);
    
    disp(['Selected group is ', deblank(groupNames(selGroups, :))]);
    disp(['Selected features are ', fmriFeatureName, ' and ', erpFeatureName]);
    
    % Get the fmri and ERP files
    if numGroups == 1
        
        erpFiles = str2mat(outputFiles(erpIndex).name);
        fmriFiles = str2mat(outputFiles(fmriIndex).name);
        
    else
        
        erpFiles = str2mat(outputFiles(erpIndex).groupFiles(selGroups).comp.name);
        fmriFiles = str2mat(outputFiles(fmriIndex).groupFiles(selGroups).comp.name);
        
    end
    
    fmriFiles = ica_fuse_fullFile('directory', outputDir, 'files', fmriFiles);
    erpFiles = ica_fuse_fullFile('directory', outputDir, 'files', erpFiles);
    
    erpInputFiles = str2mat(dataInfo(selGroups).feature(erpIndex).files.name);
    fmriInputFiles = str2mat(dataInfo(selGroups).feature(fmriIndex).files.name);
    
    
    %%%%%%%%%%%%%% Interpolate fMRI data %%%%%%%%%%%%%%%%%%%%%%%%
    featureInfo.feature_name = fmriFeatureName;
    featureInfo.groupNames = dataInfo(selGroups).name;
    featureInfo.numSubjects = numSubjects(selGroups);
    
    [fmriData, anatData] = ica_fuse_loadCompData('component_files', fmriFiles, 'slices_in_mm', slices_in_mm, ...
        'anatomical_file', ANATOMICAL_FILE, 'anatomical_view', anatomicalView, 'component_numbers', ...
        [1:numComp], 'interp_message', ['Interpolating components ', ...
        'of feature ', fmriFeatureName], 'interp_title', ['Interp comp of ', fmriFeatureName], ...
        'input_files', fmriInputFiles, 'voxels', voxels, 'feature_info', featureInfo, ...
        'mask_ind', mask_ind(fmriIndex).ind, 'outputDir', outputDir);
    
    size_anat_data = size(anatData);
    
    if length(size_anat_data) == 2
        size_anat_data = [size_anat_data, 1];
    end
    
    for nn = 1:numComp
        if nn == 1
            compMaskInd = (squeeze(fmriData(nn, :, :, :)) ~= 0);
        else
            compMaskInd = compMaskInd & (squeeze(fmriData(nn, :, :, :)) ~= 0);
        end
    end
    
    compMaskInd = find(compMaskInd);
    
    fmriData = reshape(fmriData, numComp, prod(size_anat_data));
    
    fmriData = fmriData';
    
    fmriData = fmriData(compMaskInd, :);
    
    %%%%%%%%%% End for interpolating fMRI data %%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%% Interpolate ERP data %%%%%%%%%%%%%%%%%%%%%%%%
    featureInfo.feature_name = erpFeatureName;
    featureInfo.groupNames = dataInfo(selGroups).name;
    featureInfo.numSubjects = numSubjects(selGroups);
    
    [erpData, dd, meanERPData] = ica_fuse_loadCompData('component_files', erpFiles, 'slices_in_mm', [], ...
        'anatomical_file', ANATOMICAL_FILE, 'anatomical_view', anatomicalView, 'component_numbers', ...
        [1:numComp], 'interp_message', ['Interpolating components ', ...
        'of feature ', erpFeatureName], 'interp_title', ['Interp comp of ', erpFeatureName], ...
        'input_files', erpInputFiles, 'voxels', voxels, 'feature_info', featureInfo, ...
        'mask_ind', mask_ind(erpIndex).ind, 'outputDir', outputDir);
    %%%%%%%%%% End for interpolating ERP data %%%%%%%%%%%%%%%%%%%
    
    
    if size(load(deblank(erpInputFiles(1, :))), 2) == 1
        erpUnits = 'index';
    else
        erpUnits = 'milliseconds';
    end
    
    % ERP indices
    xAxis = squeeze(erpData(1, :, 1));
    
    if ((~strcmpi(erpUnits, 'index')) && strcmpi(deblank(timeStr{timeVal}), 'secs'))
        xAxis = xAxis*1000;
    end
    
    try
        timeRange = str2num(timeRange);
    catch
    end
    
    erpData = squeeze(erpData(:, :, 2)); % Components by time
    
    disp(str2mat(['Components selected to get the weighted fMRI (Voxels by time) component are: '], ...
        num2str(selComp)));
    
    if (length(timeRange) == 2)
        timeInd = (xAxis >= min(timeRange)).*(xAxis <= max(timeRange));
        timeInd = find(timeInd == 1);
    else
        if (length(timeRange) == 1)
            timeInd = find(xAxis >= timeRange(1));
        else
            timeInd = (1:length(xAxis));
        end
    end
    %     timeInd = find(xAxis >= 1);
    if (isempty(timeInd))
        error('Specified time range is invalid. Please select the correct time range.');
    end
    
    xind = timeInd(1:step:end);
    erpData = erpData(selComp, :);
    erpic2 = erpData(:, xind);
    
    % Weighted IC
    weightedIC = fmriData(:, selComp)*abs(erpic2); % Weighted component (Voxels by time)
    
    weightedIC = weightedIC / std(weightedIC(:));
    
    clear fmriData;
    
    fll = (min([erpic2(:);meanERPData(:)]):0.01:max([erpic2(:);meanERPData(:)]));
    
    msgStr = 'Creating movie may take some time. Please don''t open any windows as it will affect creating movie.';
    helpHandle = helpdlg(msgStr, 'Creating movie ...');
    
    pause(3);
    
    try
        delete(helpHandle);
    catch
    end
    
    drawnow;
    
    xOrigin = 0.01; yOrigin = 0.15;
    xOffset = 0.04; yOffset = 0.075;
    axesWidth = 0.46;
    axesHeight = 0.5;
    fontSize = 0.05;
    cmap = ica_fuse_getColormap(imageValues, 1);
    handles = figure('color', FIG_BG_COLOR, 'DefaultAxesColor', FIG_BG_COLOR, 'name', ...
        [fmriFeatureName, '-', erpFeatureName, ' movie'], 'colormap', cmap);
    disp(msgStr);
    % Initialise movie object
    mObject = repmat(struct('cdata', [], 'colormap', []), length(xind), 1);
    for nn = 1:length(xind)
        figure(handles);
        h1 = subplot(1, 2, 1);
        compImg = zeros([1, size_anat_data]);
        compImg(compMaskInd) = (weightedIC(:, nn))';
        [compImg] = ica_fuse_applyDispParameters('image', compImg, 0, imageValues, thresh);
        compImg = reshape(compImg, [1, size_anat_data]);
        %make_composite(compImg, anatData, thresh);
        [funcImg, minICAIm, maxICAIm, minInterval, maxInterval] = ica_fuse_create_composite(anatData, compImg, ...
            imageValues, anatomicalView, size_anat_data, thresh);
        funcImg = squeeze(funcImg);
        ImageAxis = image(funcImg, 'parent', h1, 'CDataMapping', 'scaled');
        set(h1, 'CLIM', [minInterval, 2*maxInterval]);
        axis(h1, 'image');
        
        axes1Pos = [xOrigin, yOrigin, axesWidth, axesHeight];
        set(h1, 'position', axes1Pos);
        set(h1, 'fontname', UI_FONT_NAME, 'fontunits', UI_FONT_UNITS, 'fontsize', UI_FONT_SIZE);
        
        if ~strcmpi(erpUnits, 'index')
            titleStr = [num2str(round(xAxis(xind(nn)))), ' ms'];
        else
            titleStr = [num2str(round(xAxis(xind(nn)))), ' index'];
        end
        
        title(titleStr, 'color', FIG_FG_COLOR, 'parent', h1);
        
        % ERP Plot
        h2 = subplot(1, 2, 2);
        axes2Pos = axes1Pos;
        axes2Pos(1) = axes2Pos(1) + axes2Pos(3) + xOffset;
        %axes2Pos(3) = 0.5;
        set(h2, 'position', axes2Pos);
        set(h2, 'xminorTick', 'on');
        grid on;
        plot(xAxis, meanERPData, 'y.', 'parent', h2); hold on;
        plot(xAxis, erpData', 'parent', h2);
        hold on;
        %axis tight;
        lineH = plot(xAxis(xind(nn)), fll, 'w.', 'parent', h2);
        hold off;
        axis(h2, 'tight');
        set(h2, 'fontname', UI_FONT_NAME, 'fontunits', UI_FONT_UNITS, 'fontsize', UI_FONT_SIZE);
        set(h2, 'YColor', [1 1 1], 'XColor', [1 1 1]);
        
        xlabel(erpUnits, 'parent', h2);
        ylabel('microvolts', 'parent', h2);
        
        mObject(nn) = getframe(handles);
        
        drawnow;
        
    end
    
    drawnow;
    delete(handles);
    disp('Done creating movie');
    movie2avi(mObject, movieFileName, 'fps', 5, 'compression', 'none');
    disp(str2mat([fmriFeatureName, '-', erpFeatureName, ' movie is stored in file: '], movieFileName));
    
catch
    ica_fuse_displayErrorMsg;
end
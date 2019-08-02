function ica_fuse_display(fusionFile)
% Display GUI for fusion

% Load defaults
ica_fuse_defaults;

global FUSION_INFO_MAT_FILE;
global UI_FONT_NAME;
global UI_FONT_SIZE;
global ANATOMICAL_FILE;

% Display Defaults
global IMAGE_VALUES;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGES_PER_FIGURE;
global ANATOMICAL_PLANE;
global SORT_COMP;


% Select fusion file
if ~exist('fusionFile', 'var')
    fusionFile = ica_fuse_selectEntry('title', 'Select fusion information file for displaying components', ...
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

anatomicalView = 'axial';

if ~isfield(fusionInfo, 'run_analysis')
    error('Please run the analysis to view the results');
end

if ~isfield(fusionInfo.run_analysis, 'outputFiles')
    error('No output component images to display');
end

%%%%%%% Get the required parameters from run_analysis field %%%%%%%%%%%

% Output directory
displayParameters.outputDir = outputDir;


% Flip parameter for analyze images
flip_analyze_images = [];
if isfield(fusionInfo.run_analysis, 'flip_analyze_images')
    flip_analyze_images = fusionInfo.run_analysis.flip_analyze_images;
end


displayMsg = 'Opening Display Window ...';

disp(displayMsg);

fprintf('\n');

helpHandle = helpdlg(displayMsg, displayMsg);

displayParameters.flip_analyze_images = flip_analyze_images;

% Get the information from the fusion file
numComp = fusionInfo.run_analysis.numComp;

% Number of components
displayParameters.numComp = numComp;

% Output prefix
displayParameters.output_prefix = fusionInfo.run_analysis.prefix;

% Get the feature names
displayParameters.featureNames = str2mat(fusionInfo.run_analysis.dataInfo(1).feature.name);

% Get the groups names
displayParameters.groupNames = str2mat(fusionInfo.run_analysis.dataInfo.name);

% Get the ouput files
displayParameters.outputFiles = fusionInfo.run_analysis.outputFiles;

% Get the dataInfo
displayParameters.dataInfo = fusionInfo.run_analysis.dataInfo;

% Get the PCA files
displayParameters.pcaFiles = fusionInfo.run_analysis.pcaFiles;

% Get the ICA Files
displayParameters.icaFiles = fusionInfo.run_analysis.icaFiles;

% Back reconstruct files
displayParameters.backReconstructFiles = fusionInfo.run_analysis.backReconstructFiles;

% Number of groups
displayParameters.numGroups = fusionInfo.run_analysis.numGroups;

% Number of features
displayParameters.numFeatures = fusionInfo.run_analysis.numFeatures;

% voxel indices
displayParameters.mask_ind = fusionInfo.run_analysis.mask_ind;
[displayParameters.mask_ind] = ica_fuse_form_maskInd(displayParameters.mask_ind, displayParameters.dataInfo);

% Number of subjects
numSubjects = fusionInfo.run_analysis.numSubjects;


if length(numSubjects) ~= length(displayParameters.dataInfo)
    % Loop over groups
    numSubjects = repmat(numSubjects, 1, length(displayParameters.dataInfo) );
end

displayParameters.numSubjects = numSubjects;

% Full file path of the fusion file
displayParameters.fusionFile = fusionFile;

% Draw a figure with two listboxes, a check box, Convert to z scores,
% Threshold, Slices, Anatomical Plane
displayTag = 'DisplayGUI_Fusion';

%%%%%% Delete any previous figures of display GUI %%%%%
displayH = findobj('tag', displayTag);

for ii = 1:length(displayH)
    delete(displayH(ii));
end
%%%%%% end for deleting previous figures of display GUI %%%%

% Display GUI figure
[displayHandle] = ica_fuse_getGraphics('Display GUI', 'displaygui', displayTag, 'off');
set(displayHandle, 'menubar', 'none');

% Plot Title for the figure
% parameters to plot in a menu

% offsets
xOffset = 0.04; yOffset = 0.03;


%%%%%%%%%%%%% Draw Title here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title color
titleColor = [0 0.9 0.9];
% fonts
titleFont = 13;
titleAxesH = axes('Parent', displayHandle, 'position', [0 0 1 1], 'visible', 'off');
axis(titleAxesH, 'off');
xPos = 0.5; yPos = 0.97;
text(xPos, yPos, 'Display GUI', 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', ...
    'fontsize', titleFont, 'HorizontalAlignment', 'center', 'FontName', UI_FONT_NAME);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot display button
buttonWidth = 0.2; buttonHeight = 0.05;
displayButtonPos = [0.75 yOffset buttonWidth buttonHeight];


displayButtonH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', displayButtonPos, 'string', 'Display', ...
    'tag', 'display_button', 'callback', {@displayCallback, displayHandle});

%%%%%%%%%%%%%% Component Listbox %%%%%%%%%%%%%%%%%%%%

compListWidth = 0.20; listboxWidth = 0.4; yPos = 0.92;

% Plot textbox
textBoxWidth = compListWidth; textboxHeight = 0.05;
textboxPos = [xOffset, yPos - yOffset, textBoxWidth, textboxHeight];

% Plot component text
compTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Comp No.');

% form component string
for nComp = 1:numComp
    compStr(nComp).name = num2str(nComp);
end


listboxPos = textboxPos;
listboxPos(4) = 0.22;
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);
% Plot component listbox
compListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', str2mat(compStr.name), 'value', [1:numComp], ...
    'min', 0, 'max', 2,  'tag', 'selComp');
%%%%%%%%%%%%% End for plotting component listbox and text %%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% Features Text and Listbox %%%%%%%%%%%%%%%%%%%%%%
textboxPos(1) = 1 - xOffset - listboxWidth;
textboxPos(3) = listboxWidth;

% Plot Feature text
featureTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Feature');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
% Plot feature listbox
featureListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', displayParameters.featureNames, 'value', ...
    [1:size(displayParameters.featureNames, 1)], 'min', 0, 'max', 2, 'tag', 'selFeature');

% Textbox width
clear textboxPos;
%%%%%%%%%%%%%%%%%%%%%%% End for plotting features listbox and text %%%%%%%

% Plot drop down box for sorting joint components
textPos = [xOffset, listboxPos(2) - yOffset, 0.6, 0.05];
textToPlot = {'Do You Want To Sort Components?'};

sortTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textPos, 'String', textToPlot);
[textToPlot, newPos] = textwrap(sortTextH, textToPlot);
textPos(4) = newPos(4);
textPos(2) = listboxPos(2) - yOffset - textPos(4);

set(sortTextH, 'string', textToPlot);
set(sortTextH, 'position', textPos);

popupPos = textPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = 0.2;
popupPos(4) = 0.05;
popupPos(2) = textPos(2) + 0.5*textPos(4) - 0.5*popupPos(4);
sortPopupH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupPos, 'String', str2mat('No', 'Yes'), 'tag', ...
    'sort_components', 'callback', {@sortCompCallback, displayHandle});


% plot display defaults
display_defaultsMenu = uimenu('parent', displayHandle, 'label', 'Display Defaults', 'callback', ...
    {@display_defaults_callback, displayHandle});

[controlPara, displayParameters] = ica_fuse_get_display_parameters(displayParameters, [], 'off');

set(compListH, 'userdata', controlPara);

set(displayHandle, 'userdata', displayParameters);

anatomicalPlaneH = findobj(displayHandle, 'tag', 'anatomical_plane');
set(anatomicalPlaneH , 'callback', {@anatomicalPopupCallback, displayHandle});

% Plot Utils menu
utilsMenuH = uimenu('parent', displayHandle, 'label', 'Utilities');
resultsSummaryH = uimenu('parent', utilsMenuH, 'label', 'Results Summary', 'callback', ...
    {@resultsSummaryCallback, displayHandle});
movieMenuH = uimenu('parent', utilsMenuH, 'label', 'Create ERP-fMRI Movie', 'callback', ...
    {@createMovieCallback, displayHandle});


% help on display GUI
fitHelpTitle = uimenu('parent', displayHandle, 'label', 'FIT-Help');
fitHelpMenu = uimenu(fitHelpTitle, 'label', 'Display GUI', 'callback', ...
    'ica_fuse_openHTMLHelpFile(''fit_display_gui.htm'');');

% Make the graphics visible after plotting all controls
set(displayHandle, 'visible', 'on');

try
    delete(helpHandle);
catch
end


%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CALLBACKS %%%%%%%%%%%%%%%%%%%

function displayCallback(hObject, event_data, handles)
% hObject - Display push button
% handles - Display GUI figure
% Purpose: Display the components by features


try
    
    set(handles, 'pointer', 'watch');
    ica_fuse_defaults;
    global ANATOMICAL_FILE;
    
    % Component listbox
    compListH = findobj(handles, 'tag', 'selComp');
    set(compListH, 'enable', 'on');
    %ica_fuse_get_values_display_controls(handles, compListH);
    
    displayParameters = get(handles, 'userdata');
    
    fusionFile = displayParameters.fusionFile;
    
    % Number of groups and features
    numGroups = displayParameters.numGroups;
    
    numFeatures = displayParameters.numFeatures;
    
    % Number of components
    numComp = displayParameters.numComp;
    
    % Output directory
    outputDir = displayParameters.outputDir;
    
    cd(outputDir);
    
    % Output prefix
    output_prefix = displayParameters.output_prefix;
    
    % data info
    dataInfo = displayParameters.dataInfo;
    % Output files
    outputFiles = displayParameters.outputFiles;
    
    % Pca files
    pcaFiles = displayParameters.pcaFiles;
    
    % ica files
    icaFiles = displayParameters.icaFiles;
    
    % Back reconstruct files
    back_reconstruct_files = displayParameters.backReconstructFiles;
    
    % Voxels
    mask_ind = displayParameters.mask_ind;
    
    % Modalities
    modalities = str2mat(dataInfo(1).feature.modality);
    
    % Number of subjects
    numSubjects = displayParameters.numSubjects;
    
    [dims, voxels] = ica_fuse_getFeatureDIM(mask_ind);
    
    % Flip parameter for analyze images
    flip_analyze_images = displayParameters.flip_analyze_images;
    
    % Get the display fields
    
    % Selected components
    [answerString, selectedComp] = ica_fuse_get_value_uicontrol(handles, 'selComp');
    displayParameters.selectedComp = selectedComp;
    
    %%%%% Sort components %%%%
    sortCompH = findobj(handles, 'tag', 'sort_components');
    sortOptionsString = get(sortCompH, 'string');
    sortOptionsVal = get(sortCompH, 'value');
    sortSelString = deblank(sortOptionsString(sortOptionsVal, :));
    
    % Call sorting function
    doSorting(handles);
    
    if strcmpi(sortSelString, 'yes')
        displayParameters = get(handles, 'userdata');
        sortResults = displayParameters.sortResults;
        selectedComp = sortResults.sorted_comp;
        displayParameters.selectedComp = selectedComp;
        sortCompData = sortResults.sortCompData;
        plotType = sortResults.plotType;
    end
    
    %     if isfield(displayParameters, 'sortResults')
    %         sortResults = displayParameters.sortResults;
    %         selectedComp = sortResults.sorted_comp;
    %         displayParameters.selectedComp = selectedComp;
    %         sortCompData = sortResults.sortCompData;
    %         plotType = sortResults.plotType;
    %     end
    
    %%%%% Order components based on their peak of EEG signal if present %%%%%%
    if (length(selectedComp) == numComp) & strcmpi(sortSelString, 'no')
        % Check if eeg data is present
        eegIndices = strmatch('eeg', modalities, 'exact');
        if ~isempty(eegIndices)
            eegIndex = eegIndices(1);
            disp(['Ordering components based on peak of feature ', outputFiles(eegIndex).feature_name, ' ...']);
            fprintf('\n');
            % Load eeg data
            eeg_data = ica_fuse_loadData(outputFiles(eegIndex).name);
            % Get only the yaxis
            eeg_data = squeeze(eeg_data(:, 2, :));
            maxIndices = zeros(1, numComp);
            stdComp = zeros(1, numComp);
            % Loop over number of components
            for nComp = 1:numComp
                [maxValue, maxIndices(nComp)] = max(abs(eeg_data(:, nComp)));
                stdComp(nComp) = std(eeg_data(:, nComp));
            end
            % End loop over number of components
            [maxIndVal, selectedComp] = sort(maxIndices);
            clear eeg_data;
        end
        % end for checking eeg data
        
    end
    %%%%%% end for ordering the components based on their peak %%%%%%
    
    % Selected Features
    [selectedFeature, selectedFeatureVal] = ica_fuse_get_value_uicontrol(handles, 'selFeature');
    selectedFeature = deblank(selectedFeature(selectedFeatureVal, :));
    displayParameters.selectedFeature = selectedFeature;
    displayParameters.selectedFeatureVal = selectedFeatureVal;
    
    % Retrieve display parameters
    convertToZ = displayParameters.convert_to_z;
    imageValues = displayParameters.image_values;
    images_per_figure = displayParameters.images_per_figure;
    threshold = abs(displayParameters.z_threshold);
    anatomical_plane = lower(displayParameters.anatomical_plane);
    slices_in_mm = displayParameters.slices_in_mm;
    
    if (min(threshold) == max(threshold))
        threshold = threshold(1);
    end
    
    
    cm = [];
    text_left_right = [];
    for nOutputFiles = 1:length(selectedFeatureVal)
        % component files
        componentFiles = ica_fuse_fullFile('directory', outputDir, 'files', ...
            outputFiles(selectedFeatureVal(nOutputFiles)).name);
        if isfield(outputFiles(selectedFeatureVal(nOutputFiles)), 'groupFiles')
            groupCompFiles = outputFiles(selectedFeatureVal(nOutputFiles)).groupFiles;
        else
            groupCompFiles = [];
        end
        feature_name = outputFiles(selectedFeatureVal(nOutputFiles)).feature_name;
        
        featureInfo.feature_name = feature_name;
        featureInfo.groupNames = str2mat(dataInfo.name);
        featureInfo.numSubjects = numSubjects;
        %featureInfo.numSubjects = size(str2mat(dataInfo(1).feature(1).files.name), 1);
        
        % Get the input files
        for nInputFiles = 1:length(dataInfo)
            
            if nInputFiles == 1
                feature_input_files = str2mat(dataInfo(nInputFiles).feature(selectedFeatureVal(nOutputFiles)).files.name);
            else
                feature_input_files = str2mat(feature_input_files, ...
                    str2mat(dataInfo(nInputFiles).feature(selectedFeatureVal(nOutputFiles)).files.name));
            end
            
        end
        % end for getting feature files
        
        [compData, anatData, meanData, groupCompData, meanDataLegend, groupCompLegend, text_left_right, ...
            plotType, HInfo] = ica_fuse_loadCompData('component_files', componentFiles, ...
            'group_component_files', groupCompFiles, 'slices_in_mm', slices_in_mm, 'anatomical_file', ...
            ANATOMICAL_FILE, 'anatomical_view', anatomical_plane, 'component_numbers', ...
            selectedComp, 'interp_message', ['Interpolating components ', ...
            'of feature ', feature_name], 'interp_title', ['Interp comp of ', feature_name], ...
            'input_files', feature_input_files, 'voxels', voxels, 'feature_info', featureInfo, ...
            'mask_ind', mask_ind(selectedFeatureVal(nOutputFiles)).ind, 'outputDir', outputDir, ...
            'flip_analyze_images', flip_analyze_images, 'modality', deblank(modalities(selectedFeatureVal(nOutputFiles), :)));
        
        
        if ~isfield(displayParameters, 'text_left_right')
            displayParameters.text_left_right = text_left_right;
        end
        
        clear feature_input_files;
        clear maxICAIm; clear minICAIm;
        maxICAIm = zeros(1, length(selectedComp));
        minICAIm = zeros(1, length(selectedComp));
        minInterval = 0;
        maxInterval = 0;
        
        % For images apply display parameters
        if strcmpi(plotType, 'image')
            structDIM = HInfo.DIM;
            % Apply display parameters
            compData = ica_fuse_applyDispParameters(plotType, compData, convertToZ, imageValues, threshold);
            
            % create composite image
            [compData, minICAIm, maxICAIm, minInterval, maxInterval] = ica_fuse_create_composite(anatData, compData, ...
                imageValues, anatomical_plane, structDIM, threshold);
            colorbarLim = [minInterval, maxInterval];
            % Get the associated colormap
            cm = ica_fuse_getColormap(imageValues, 1);
            clear anatData;
            
        end
        
        % Store the component data in outputData structure
        for nComp = 1:length(selectedComp)
            outputData(nOutputFiles).CompData(nComp).data = squeeze(compData(nComp, :, :));
            compIndex = ica_fuse_returnFileIndex(selectedComp(nComp));
            if ~isempty(groupCompData)
                outputData(nOutputFiles).CompData(nComp).groupCompData = squeeze(groupCompData(:, :, nComp, :));
            else
                outputData(nOutputFiles).CompData(nComp).groupCompData = [];
            end
            if ~isempty(groupCompData)
                outputData(nOutputFiles).CompData(nComp).groupCompLegend = groupCompLegend;
            else
                outputData(nOutputFiles).CompData(nComp).groupCompLegend = {};
            end
            outputData(nOutputFiles).CompData(nComp).axesTitle = ['Comp ', compIndex, ' Feature ', feature_name];
            outputData(nOutputFiles).CompData(nComp).plotType = plotType;
            outputData(nOutputFiles).CompData(nComp).maxICAIm = maxICAIm(nComp);
            outputData(nOutputFiles).CompData(nComp).minICAIm = minICAIm(nComp);
            outputData(nOutputFiles).CompData(nComp).minInterval = minInterval;
            outputData(nOutputFiles).CompData(nComp).maxInterval = maxInterval;
            if isfield(displayParameters, 'text_left_right')
                outputData(nOutputFiles).CompData(nComp).textLeftRight = displayParameters.text_left_right;
            else
                outputData(nOutputFiles).CompData(nComp).textLeftRight = [];
            end
        end
        outputData(nOutputFiles).meanData = meanData;
        outputData(nOutputFiles).meanDataLegend = meanDataLegend;
        clear meanData;
        clear meanDataLegend;
        clear featureInfo;
        clear compData;
        
    end
    
    % Change the data structure such that the components are plotted by feature
    countData = 0;
    
    for nComp = 1:length(selectedComp)
        
        if exist('sortCompData', 'var')
            % Add Mixing Coefficients at the top
            countData = countData + 1;
            plotData(countData).data = sortCompData(nComp).data;
            compIndex = ica_fuse_returnFileIndex(selectedComp(nComp));
            plotData(countData).axesTitle = sortCompData(nComp).title;
            plotData(countData).plotType = sortResults.plotType;
            plotData(countData).colorbarMinMaxText = sortCompData(nComp).colorbarText;
            plotData(countData).colorbarLim = sortCompData(nComp).colorbarLim;
            plotData(countData).groupNames = sortResults.selGroupNames;
            plotData(countData).selFeatureNames = sortResults.selFeatureNames;
            % Selected groups, features, component
            plotData(countData).selGroupsVal = sortResults.selGroupsVal;
            plotData(countData).selFeaturesVal = sortResults.selFeaturesVal;
            plotData(countData).compNum = selectedComp(nComp);
            plotData(countData).fusionFile = fusionFile;
            plotData(countData).threshold = threshold;
            %             if isfield(sortCompData, 'histData')
            %                 plotData(countData).histData = sortCompData(nComp).histData;
            %             end
        end
        
        % Plot For a component the features
        for nOutputFiles = 1:length(outputData)
            countData = countData + 1;
            plotData(countData).data = outputData(nOutputFiles).CompData(nComp).data;
            plotData(countData).groupCompData = outputData(nOutputFiles).CompData(nComp).groupCompData;
            plotData(countData).groupCompLegend = outputData(nOutputFiles).CompData(nComp).groupCompLegend;
            plotData(countData).meanData = outputData(nOutputFiles).meanData;
            plotData(countData).meanDataLegend = outputData(nOutputFiles).meanDataLegend;
            plotData(countData).axesTitle = outputData(nOutputFiles).CompData(nComp).axesTitle;
            plotData(countData).plotType = outputData(nOutputFiles).CompData(nComp).plotType;
            minICAIm = outputData(nOutputFiles).CompData(nComp).minICAIm;
            maxICAIm = outputData(nOutputFiles).CompData(nComp).maxICAIm;
            plotData(countData).colorbarMinMaxText = str2mat(num2str(minICAIm), num2str(maxICAIm));
            plotData(countData).colorbarLim = [outputData(nOutputFiles).CompData(nComp).minInterval, ...
                outputData(nOutputFiles).CompData(nComp).maxInterval];
            plotData(countData).groupNames = '';
            plotData(countData).textLeftRight = outputData(nOutputFiles).CompData(nComp).textLeftRight;
        end
        
    end
    clear outputData;
    
    ica_fuse_display_features('plot_data', plotData, 'color_map', cm, 'title_color', 'c', 'time_course_color', '.-c', ...
        'ica_loading_color', {'g'; 'r'; 'c'; 'm'; 'b'}, 'number_per_figure', images_per_figure, 'slice_plane', anatomical_plane);
    
    set(handles, 'userdata', displayParameters);
    set(handles, 'pointer', 'arrow');
    
catch
    
    set(handles, 'pointer', 'arrow');
    ica_fuse_errorDialog(lasterr, 'Display GUI Error', 'modal');
    %rethrow(lasterror);
    ica_fuse_displayErrorMsg;
    
end


%%%%%%% Load anatomical callback %%%%%%
function anatomicalPopupCallback(hObject, event_data, handles)

ica_fuse_setString_slices(hObject, handles);


function display_defaults_callback(hObject, event_data, handles)
% get the display defaults

displayParameters = get(handles, 'userdata');

compListH = findobj(handles, 'tag', 'selComp');

controlPara = get(compListH, 'userdata');

[controlPara, displayParameters] = ica_fuse_get_display_parameters(displayParameters, controlPara, 'on');

set(compListH, 'userdata', controlPara);

set(handles, 'userdata', displayParameters);



function doSorting(handles)
% Sorting joint components

hObject = findobj(handles, 'tag', 'sort_components');

answerString = get(hObject, 'string');
answerVal = get(hObject, 'value');

selectedString = deblank(answerString(answerVal, :));
% Add sortResults field the displayParameters structure
compListH = findobj(handles, 'tag', 'selComp');
displayParameters = get(handles, 'userdata');

if strcmpi(selectedString, 'yes')
    
    % Select all the components
    compString = get(compListH, 'string');
    set(compListH, 'value', [1:1:size(compString, 1)]);
    
    fusionFile = displayParameters.fusionFile;
    
    % Get the sorting results
    sortResults = ica_fuse_sortingGUI(fusionFile);
    
    set(compListH, 'enable', 'off');
    
    displayParameters.sortResults = sortResults;
    fprintf('\n');
    %disp('Click on display button to display the joint components ...');
else
    set(compListH, 'enable', 'on');
    if isfield(displayParameters, 'sortResults')
        % Remove the sortResults structure if no option is selected for
        % sort joint components
        displayParameters = rmfield(displayParameters, 'sortResults');
    end
    
end

set(handles, 'userdata', displayParameters);


function sortCompCallback(hObject, event_data, handles)
% Disable or enable component numbers listbox depending on whether the
% answer for drop down box is yes or no

answerString = get(hObject, 'string');
answerVal = get(hObject, 'value');

selectedString = deblank(answerString(answerVal, :));
compListH = findobj(handles, 'tag', 'selComp');

if strcmpi(selectedString, 'yes')
    % Select all the components
    compString = get(compListH, 'string');
    set(compListH, 'value', [1:1:size(compString, 1)]);
    set(compListH, 'enable', 'off');
else
    set(compListH, 'enable', 'on');
end


function createMovieCallback(hObject, event_data, handles)
% Create Movie Callback

displayParameters = get(handles, 'userdata');

fusionFile = displayParameters.fusionFile;
ica_fuse_create_movie_erp_fmri(fusionFile, displayParameters);


function resultsSummaryCallback(hObject, event_data, handles)
% Results summary callback
%

formatName = questdlg('Select results format', 'Results format', 'HTML', 'PDF', 'HTML');

if (~isempty(formatName))
    
    displayParameters = get(handles, 'userdata');
    param_file = displayParameters.fusionFile;
    
    load(param_file);
    
    if (~exist('fusionInfo', 'var'))
        error('Selected file is not a valid fusion MAT file');
    end
    
    outDir = fullfile(fileparts(param_file), [fusionInfo.setup_analysis.prefix, '_jica_results']);
    opts.outputDir = outDir;
    opts.showCode = false;
    opts.useNewFigure = false;
    opts.format = lower(formatName);
    opts.createThumbnail = true;
    if (strcmpi(opts.format, 'pdf'))
        opt.useNewFigure = false;
    end
    assignin('base', 'param_file', param_file);
    assignin('base', 'displayParameters', displayParameters);
    opts.codeToEvaluate = 'ica_fuse_jica_summary(param_file, displayParameters);';
    %publish('icatb_gica_html_report', 'outputDir', outDir, 'showCode', false, 'useNewFigure', false);
    disp('Generating reults summary. Please wait ....');
    drawnow;
    publish('ica_fuse_jica_summary', opts);
    
    close all;
    
    if (strcmpi(opts.format, 'html'))
        ica_fuse_openHTMLHelpFile(fullfile(outDir, 'ica_fuse_jica_summary.html'));
    else
        open(fullfile(outDir, 'ica_fuse_jica_summary.pdf'));
    end
    
    disp('Done');
    
end



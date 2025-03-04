function fusionFile = ica_fuse_anyway_setup_analysis(fusionFile)
% function for setting up parameters required for anyway fusion analysis
% parameters required are:
% 1. output directory
% 2. prefix for output files
% 3. data - of structure format
% 4. mask - mask for fmri data
% 5. Number of components
% All these parameters are defined in ica_fuse_anyway_fusion_params.m file.
%
% Output: fusionFile - fusionFile is a MAT file containing the necessary
% information for running the analysis

ica_fuse_defaults;

% Set inputFile variable to be empty
if ~exist('fusionFile', 'var')
    fusionFile = [];
end


% DEFINE PARAMETERS THAT ARE GOING TO BE PLOTTED.
% IN CASE OF BATCH ANALYSIS PARAMETERS WILL BE READ FROM
% THE FILE AND THE FUNCTION CALLBACKS WILL BE EXECUTED AUTOMATICALLY
inputText = ica_fuse_anyway_fusion_params;

handles_data.analysisType = 'gui';

% If the fusion file is not present select the output directory
if isempty(fusionFile)
    % select output directory
    outputDir = ica_fuse_selectEntry('typeEntity', 'directory', 'title', ...
        'Select output directory for Anyway Fusion analysis');
else
    formStr = '_anyway_fusion.mat';
    [outputDir, fileN, file_extn] = fileparts(fusionFile);
    inputString = [fileN, file_extn];
    matchIndex = ica_fuse_findstr(inputString, formStr);
    prefixToSet = inputString(1:matchIndex(1)-1);
end

handle_visibility = 'on';

% end for reading output directory for the analysis


if isempty(outputDir)
    error('Output directory for the data fusion analysis is not selected');
end

% change to the resulting output directory
cd(outputDir);


figureTag = 'Setup Fusion For Anyway Fusion Analysis';
okTag = 'Done'; cancelTag = 'Cancel';

figHandle = findobj('tag', figureTag);

if ~isempty(figHandle)
    for nFig = 1:length(figHandle)
        delete(figHandle(nFig));
    end
end

plotHelp = 1;
% Plot all the controls in a figure
[InputHandle] = ica_fuse_plot_controls_fig(inputText, figureTag, handle_visibility, okTag, cancelTag, plotHelp);


for nH = 1:length(inputText)
    helpHandle = findobj(InputHandle, 'tag', ['help_', inputText(nH).tag]);
    set(helpHandle, 'callback', {@showHelpDialog, InputHandle});
end

okHandle = findobj(InputHandle, 'tag', okTag);
cancelHandle = findobj(InputHandle, 'tag', cancelTag);

% fusionInfo is the variable that will be used for saving the parameter
% information
fusionInfo.setup_analysis.outputDir = outputDir;

% store the input text and fusionInfo in figure data
handles_data.inputText = inputText;
handles_data.fusionInfo = fusionInfo;

% set the input figure user data
set(InputHandle, 'userdata', handles_data);

%%%%%%%%%%%%%% SET UP THE OBJECT CALLBACKS IN ORDER %%%%%%%%%%%%%%%%%%%%%%%%%

% Output prefix callback
set(findobj(InputHandle, 'tag', 'prefix'), 'callback', {@prefixCallback, InputHandle});

% Set callback for data sets
set(findobj(InputHandle, 'tag', 'dataInfo'), 'callback', {@dataSelectCallback, InputHandle});

% mask callback
set(findobj(InputHandle, 'tag', 'maskFile'), 'callback', {@selectMaskCallback, InputHandle});

% number of components callback
set(findobj(InputHandle, 'tag', 'numComp'), 'callback', {@numCompCallback, InputHandle});

% done callback
set(okHandle, 'callback', {@applyCallback, InputHandle});

% cancel callback
set(cancelHandle, 'callback', {@closeCallback, InputHandle});


referenceH = findobj(InputHandle, 'tag', 'reference');
set(referenceH, 'callback', {@referenceICACallback, InputHandle});

if exist('prefixToSet', 'var')
    prefixH = findobj(InputHandle, 'tag', 'prefix');
    set(prefixH, 'string', prefixToSet);
    prefixCallback(prefixH, [], InputHandle);
end



% GUI Analysis
try
    % Wait for the input from the user in case of GUI
    set(InputHandle, 'visible', 'on');
    waitfor(InputHandle);
catch
    
end
% end for executing object callbacks in case of batch file analysis

global APPNAME;
appName = APPNAME;

clear global APPNAME;

if ~isempty(appName)
    % Get the application data
    if isappdata(0, appName)
        stackInfo = getappdata(0, appName);
        if isfield(stackInfo, 'fusionFile')
            fusionFile = stackInfo.fusionFile;
            rmappdata(0, appName);
        else
            return;
        end
        clear stackInfo;
    end
end

%%%%%%%%%%%%%%%% DEFINE FUNCTION CALLBACKS %%%%%%%%%%%%%%%%%%%%

% OUTPUT PREFIX CALLBACK
function prefixCallback(handleObj, event_data, handles)

output_prefix = get(handleObj, 'string');

tagObj = get(handleObj, 'tag');

% check the output prefix
[status, message] = ica_fuse_errorCheck(output_prefix, 'output_prefix');

if status == 0
    error(message);
end

handles_data = get(handles, 'userdata');
handles_data.fusionInfo.setup_analysis = setfield(handles_data.fusionInfo.setup_analysis, tagObj, output_prefix);

% set the figure data
set(handles, 'userdata', handles_data);

% update the selecting data files control
% check if the subject file exists or not
% and then update the controls
%if isempty(handles_data.inputFile)
% input file
inputFile = [];
oldDir = handles_data.fusionInfo.setup_analysis.outputDir;
dataFile = [handles_data.fusionInfo.setup_analysis.prefix, '_anyway_fusion.mat'];
dataFile = fullfile(oldDir, dataFile);
dataTag = findobj(handles, 'tag', 'dataInfo');
if exist(dataFile, 'file') && isempty(inputFile)
    load(dataFile);
    drawnow;
    fusionInfo.setup_analysis.outputDir = oldDir;
    if isfield(fusionInfo.setup_analysis, 'dataInfo')
        if ~isempty(fusionInfo.setup_analysis.dataInfo)
            set(dataTag, 'style', 'popup', 'string', str2mat('Yes', 'No'), 'value', 1);
            handles_data.fusionInfo = fusionInfo;
            % set the figure data
            set(handles, 'userdata', handles_data);
            dataSelectCallback(dataTag, [], handles);
        end
    end
    clear fusionInfo;
end
%end



function dataSelectCallback(handleObj, event_data, handles)
% select the data

set(handles, 'pointer', 'watch');

try
    
    % get the structure for the figure
    handles_data = get(handles, 'userdata');
    fusionInfo = handles_data.fusionInfo;
    
    % get the inputText
    inputText = handles_data.inputText;
    
    allTags = char(inputText.tag);
    
    % disable all the controls in the input figure
    ica_fuse_enable_control(allTags, handles, 'inactive');
    
    % get the directory for the analysis
    oldDir = fusionInfo.setup_analysis.outputDir;
    
    % get the prefix of the output string
    fusionInfo.setup_analysis.prefix = get(findobj(handles, 'tag', 'prefix'), 'string');
    
    % All the entered parameters are going to be stored in this file
    fusionInfo.setup_analysis.fusionFile = [fusionInfo.setup_analysis.prefix, '_anyway_fusion.mat'];
    newParamFile = fullfile(oldDir, fusionInfo.setup_analysis.fusionFile);
    
    fusionPrefix = fusionInfo.setup_analysis.prefix;
    appName = [fusionPrefix, 'setupAppData'];
    
    % get the style of the UIcontrol
    getStyle = get(handleObj, 'style');
    
    if strcmpi(getStyle, 'pushbutton')
        getSubjects = 'no';
    else
        getPopValue = get(handleObj, 'value');
        getPopString = get(handleObj, 'string');
        getSubjects = lower(deblank(getPopString(getPopValue, :)));
    end
    
    % First determine whether the subjects file exists or not
    % If it doesn''t exist then show the user the error dialog
    % else load the subject file
    if strcmp(getSubjects, 'yes')
        if ~exist(newParamFile, 'file')
            error('Fusion Info file doesn''t exist. Please use the other option');
        end
    end
    
    
    if ~strcmp(getSubjects, 'yes')
        % Open a input dialog with the following parameters
        
        % selecting data
        dataInfo = ica_fuse_dataSelection;
        
        %%%%%%%%%%%%%%%% Printing information to a file %%%%%%%%%%%%%
        file_name = fullfile(oldDir, [fusionInfo.setup_analysis.prefix, '_anyway_fusion.txt']);
        
        % Print selected data information
        ica_fuse_write_selected_data(dataInfo, file_name);
        
        fusionInfo.setup_analysis.dataInfo = dataInfo;
        
        % save the information for running the analysis
        ica_fuse_save(newParamFile, 'fusionInfo');
        
        % Remove the application data whenever the data is changed
        if isappdata(0, appName)
            rmappdata(0, appName);
        end
        
    end
    % end for selecting the data
    
    
    drawnow;
    
    % Number of groups, features and subjects.
    dataInfo = fusionInfo.setup_analysis.dataInfo;
    fusionInfo.setup_analysis.numGroups = length(dataInfo);
    fusionInfo.setup_analysis.numFeatures = length(dataInfo(1).feature);
    
    numSubjects = ica_fuse_countSubjects(dataInfo);
    
    fusionInfo.setup_analysis.numSubjects = numSubjects;
    
    % Set figure data
    cd(fusionInfo.setup_analysis.outputDir);
    
    handles_data.fusionInfo = fusionInfo;
    set(handles, 'userdata', handles_data);
    
    % show selected values for controls
    show_values_controls(handles);
    
    set(handles, 'pointer', 'arrow');
    
catch
    
    set(handles, 'pointer', 'arrow');
    
    % enable all the controls
    ica_fuse_enable_control(allTags, handles, 'on');
    
    %rethrow(lasterror);
    ica_fuse_displayErrorMsg;
    
end


function selectMaskCallback(hObject, event_data, handles)

ica_fuse_defaults;

global MRI_DATA_FILTER;

% select mask callback
handles_data = get(handles, 'userdata');
fusionInfo = handles_data.fusionInfo;
inputFile = [];

maskTag = get(hObject, 'tag');
getString = get(hObject, 'string');
getMask = get(hObject, 'value');

numFeatures = length(fusionInfo.setup_analysis.dataInfo(1).feature);
featureNames = cellstr(str2mat(fusionInfo.setup_analysis.dataInfo(1).feature.name));
modalities = cellstr(str2mat(fusionInfo.setup_analysis.dataInfo(1).feature.modality));

fusionPrefix = fusionInfo.setup_analysis.prefix;
appName = [fusionPrefix, 'setupAppData'];

if isempty(inputFile)
    % Check if the mask field exists or not
    if isfield(fusionInfo.setup_analysis, maskTag)
        if isappdata(0, appName)
            % get the answer for file
            answerFile = ica_fuse_questionDialog('title', 'Mask File', 'textbody', ...
                'Do you want to select the mask file? This will load the data again and create a mask.');
            if answerFile == 1
                rmappdata(0, appName);
            else
                return;
            end
        end
        % end for checking application data
    end
    % end for checking mask field
end

maskFile = [];

% select mask in 3D Analyze data
if getMask == 2
    if isempty(inputFile)
        % Open GUI for selecting the custom mask
        [P] = getMaskFiles(featureNames, modalities, fusionInfo.setup_analysis.dataInfo(1));
        
        if ~isempty(P)
            maskFile = P;
        else
            disp('No mask file/files selected. Setting to default mask ...');
            set(hObject, 'value', 1);
        end
    else
        inputData = ica_fuse_read_variables(inputFile, maskTag, 'cell');
        maskFile = getfield(inputData, maskTag);
        
        if length(maskFile) ~= length(modalities)
            error(['When specifying mask the length of cell array must equal the number of features']);
        end
        
        % Loop over modalities
        for nn = 1:length(modalities)
            % EEG mask
            if strcmpi(modalities{nn}, 'eeg')
                if isnumeric(maskFile{nn})
                    maskFile{nn} = num2str(maskFile{nn});
                else
                    if isempty(maskFile{nn})
                        error('Error:Mask', ['EEG indices is empty for feature %s. \nCheck the input file: %s'], ...
                            featureNames{nn}, inputFile);
                    end
                end
                
            else
                % Consider the mask as a file
                if exist(maskFile{nn}, 'file') ~= 2
                    error('Error:Mask', ['Mask file: %s doesn''t exist. \nCheck the input file: %s'], ...
                        maskFile{nn}, inputFile);
                end
            end
            % end for checking
        end
        % End loop over modalities
        clear inputData;
    end
end

% set mask field
fusionInfo.setup_analysis = setfield(fusionInfo.setup_analysis, maskTag, maskFile);

handles_data.fusionInfo = fusionInfo;
set(handles, 'userdata', handles_data);



function normalizeCallback(hObject, event_data, handles)
% Normalization callback


% get handles data
handles_data = get(handles, 'userdata');
fusionInfo = handles_data.fusionInfo;

% Get the normalization parameters
normalizeTag = get(hObject, 'tag');
normalizeVal = get(hObject, 'value');

fusionPrefix = fusionInfo.setup_analysis.prefix;
appName = [fusionPrefix, 'setupAppData'];

if isfield(fusionInfo.setup_analysis, normalizeTag)
    if isappdata(0, appName)
        
        if normalizeVal ~= (getfield(fusionInfo.setup_analysis, normalizeTag))
            % Remove the application data as the normalization scheme is
            % changed
            rmappdata(0, appName);
        end
    end
end

% set the user data
fusionInfo.setup_analysis = setfield(fusionInfo.setup_analysis, normalizeTag, normalizeVal);
handles_data.fusionInfo = fusionInfo;
set(handles, 'userdata', handles_data);


function numCompCallback(hObject, event_data, handles)

% get handles data
handles_data = get(handles, 'userdata');
fusionInfo = handles_data.fusionInfo;
compTag = get(hObject, 'tag');
compStr = get(hObject, 'string');
try
    % check if number of components is an integer
    numComp = strread(compStr, '%d');
catch
    error('Not a valid integer for number of components');
end

if numComp < 2
    error('Components must be atleast 2 to run ICA');
end

if isfield(fusionInfo.setup_analysis, 'numSubjects')
    numSubjects = fusionInfo.setup_analysis.numSubjects;
    if numComp > sum(numSubjects)
        error(['Number of components (', num2str(numComp), ...
            ') is greater than the number of data-sets (', num2str(sum(numSubjects)), ')']);
    end
end

fusionInfo.setup_analysis = setfield(fusionInfo.setup_analysis, compTag, numComp);

handles_data.fusionInfo = fusionInfo;
set(handles, 'userdata', handles_data);


function show_values_controls(handles)

% get handles data
handles_data = get(handles, 'userdata');

allTags = char(handles_data.inputText.tag);

% get fusion info
fusionInfo = handles_data.fusionInfo;

% Number of subjects
numSubjects = fusionInfo.setup_analysis.numSubjects;

% Set automatically the number of components
if min(numSubjects) < 8
    num_comp = min(numSubjects);
else
    num_comp = 8;
end

numCompH = findobj(handles, 'tag', 'numComp');

compStr = get(numCompH, 'string');

if ~isfield(fusionInfo.setup_analysis, 'numComp')
    if ~isempty(compStr)
        fusionInfo.setup_analysis.numComp = str2num(compStr);
    else
        fusionInfo.setup_analysis.numComp = num_comp;
    end
end

set(numCompH, 'enable', 'on');
set(numCompH, 'string', num2str(fusionInfo.setup_analysis.numComp(:)'));


% set mask information also
maskTag = 'maskFile';
if isfield(fusionInfo.setup_analysis, maskTag)
    maskHandle = findobj(handles, 'tag', maskTag);
    maskVal = getfield(fusionInfo.setup_analysis, maskTag);
    if ~isempty(maskVal)
        maskVal = 2;
    else
        maskVal = 1;
    end
    set(maskHandle, 'value', maskVal);
end

set(findobj(handles, 'tag', 'dataInfo'), 'style', 'popup', 'string', ...
    char('Yes', 'No'), 'value', 1);

scaleCompH = findobj(handles, 'tag', 'scale_components');
scaleComp = 'z-scores';
scaleOpts = cellstr(get(scaleCompH, 'string'));
try
    scaleComp = fusionInfo.setup_analysis.scale_components;
catch
end

scaleCompVal = strmatch(lower(scaleComp), lower(scaleOpts), 'exact');

set(scaleCompH, 'value', scaleCompVal);

% enable all the controls
ica_fuse_enable_control(allTags, handles, 'on');

handles_data.fusionInfo = fusionInfo;
set(handles, 'userdata', handles_data);

drawnow;

function dimensionEstimationCallback(hObject, event_data, handles)
% Estimate dimensionality callback

try
    
    set(handles, 'pointer', 'watch');
    ica_fuse_defaults;
    global FIG_FG_COLOR;
    global FIG_BG_COLOR;
    global AX_COLOR;
    global FWHM_VALUE;
    
    answerString = get(hObject, 'string');
    answerValue = get(hObject, 'value');
    
    handles_data = get(handles, 'userdata');
    inputFile = handles_data.inputFile;
    
    % Get the variable fusionInfo
    fusionInfo = handles_data.fusionInfo;
    
    fusionPrefix = fusionInfo.setup_analysis.prefix;
    
    % Get the mask file
    if isfield(fusionInfo.setup_analysis, 'maskFile')
        maskFile = fusionInfo.setup_analysis.maskFile;
    else
        maskFile = [];
    end
    
    fusionInfo.setup_analysis.maskFile = maskFile;
    
    % Get the normalization tag
    if ~isfield(fusionInfo.setup_analysis, 'normalize')
        fusionInfo.setup_analysis.normalize = 1;
    end
    
    normalizeVal = fusionInfo.setup_analysis.normalize;
    
    if ~isfield(fusionInfo.setup_analysis, 'dataInfo')
        error('Data is not selected for dimensionality estimation');
    end
    
    drawnow;
    
    % Estimate Components
    if strcmpi(deblank(answerString(answerValue, :)), 'yes')
        
        % Enter fwhm
        fwhm = FWHM_VALUE;
        
        % Load application data
        appName = [fusionPrefix, 'setupAppData'];
        if isappdata(0, appName)
            tmpStackInfo = getappdata(0, appName);
        else
            tmpStackInfo = ica_fuse_pre_stack_data(fusionInfo);
            setappdata(0, appName, tmpStackInfo);
        end
        
        drawnow;
        
        % Prepare data for joint ICA fusion analysis
        stackInfo = ica_fuse_prepare_data('dataInfo', fusionInfo.setup_analysis.dataInfo, 'mask_ind', tmpStackInfo.mask_ind, 'normalize_scheme', ...
            normalizeVal, 'voxels', tmpStackInfo.voxels, 'estim_data', 1);
        
        clear tmpStackInfo;
        
        drawnow;
        
        estimationInfo = stackInfo.estimationInfo;
        
        fusionInfo.setup_analysis.estimationInfo = estimationInfo;
        
        figHandle = figure('name', 'Dimensionality Estimation', 'color', FIG_BG_COLOR, 'resize', 'off', 'windowstyle', 'modal');
        
        for nPlot = 1:length(estimationInfo)
            sh = subplot(length(estimationInfo), 1, nPlot);
            plot(1:length(estimationInfo(nPlot).mdl), estimationInfo(nPlot).mdl, 'c', 'parent', sh);
            axis(sh, 'tight');
            title([estimationInfo(nPlot).feature_name , ' Estimated Comps = ', num2str(estimationInfo(nPlot).comp)], 'color', FIG_FG_COLOR, 'parent', sh);
            ylabel('MDL Value', 'parent', sh);
            set(sh, 'YColor', FIG_FG_COLOR);
            set(sh, 'XColor', FIG_FG_COLOR);
            set(sh, 'color', AX_COLOR);
        end
        
        xlabel('No. Of Subjects', 'parent', sh);
        
        drawnow;
        
        
        if (isfield(stackInfo, 'ccaEstimationInfo'))
            ccaEstimationInfo = stackInfo.ccaEstimationInfo;
            fusionInfo.setup_analysis.ccaEstimationInfo = ccaEstimationInfo;
            dimMsgString = [char(ccaEstimationInfo.name), repmat(' : ', length(ccaEstimationInfo), 1), num2str([ccaEstimationInfo.comp]')];
            dimMsgString = cellstr(dimMsgString);
            textBody = {'PCA-CCA estimation is done on pair of features and dimensionality estimation is shown below:';'';''};
            textBody = [textBody;dimMsgString];
            ica_fuse_dialogBox('title', 'PCA-CCA Estimation', 'textType', 'large', 'textbody', textBody);
            disp(char(textBody));
        end
        
        
        drawnow;
        
        % Estimate the data
        %[est_comp, mdlVal] = ica_fuse_estimate_dimension((stackInfo.data)', fwhm);
        
        clear stackInfo;
        
        %         plotX.y = mdlVal;
        %         plotX.x = (1:length(mdlVal));
        %         plotX.title = 'Plot of MDL where minimum is the estimated dimensionality.';
        %         msgString = ['The estimated independent components is found to be ', num2str(est_comp), ...
        %             ' using the MDL criteria.'];
        %         disp(msgString);
        %         fprintf('\n');
        %         if isempty(inputFile)
        %             ica_fuse_dialogBox('title', 'Estimated Components', 'textBody', msgString, 'textType', 'large', 'plotbutton', plotX);
        %         end
    end
    
    % set the user data
    handles_data.fusionInfo = fusionInfo;
    set(handles, 'userdata', handles_data, 'pointer', 'arrow');
    
catch
    
    set(handles, 'pointer', 'arrow');
    %rethrow(lasterror);
    ica_fuse_displayErrorMsg;
    
end


function zScoresCallback(hObject, event_data, handles)
% Z-scores callback

handles_data = get(handles, 'userdata');
inputFile = handles_data.inputFile;
fusionInfo = handles_data.fusionInfo;

% Get the value
value = get(hObject, 'value');
allStrings = get(hObject, 'string');

if strcmpi(deblank(allStrings(value, :)), 'z-scores')
    value = 1;
else
    value = 0;
end

tagToAttach = get(hObject, 'tag');
% set the field to structure fusionInfo
fusionInfo.setup_analysis = setfield(fusionInfo.setup_analysis, tagToAttach, value);

% Set fusionInfo field to handles or figure
handles_data.fusionInfo = fusionInfo;
set(handles, 'userdata', handles_data);

function applyCallback(hObject, event_data, handles)
% apply callback

ica_fuse_defaults;
global FUSION_INFO_MAT_FILE;

try
    
    set(handles, 'pointer', 'watch');
    
    handles_data = get(handles, 'userdata');
    
    inputFile = [];
    fusionInfo = handles_data.fusionInfo;
    
    % Check the data field first
    if ~isfield(fusionInfo.setup_analysis, 'dataInfo')
        error(['Data is not selected for the analysis']);
    end
    
    % Check the mask field
    if ~isfield(fusionInfo.setup_analysis, 'maskFile')
        fusionInfo.setup_analysis.maskFile = [];
    end
    
    % Get the normalization parameter
    if ~isfield(fusionInfo.setup_analysis, 'normalize')
        fusionInfo.setup_analysis.normalize = 1;
    end
    
    % Check the z-scores field
    if ~isfield(fusionInfo.setup_analysis, 'scale_components')
        fusionInfo.setup_analysis.scale_components = 'z-scores';
    end
    
    % Type of PCA
    typePCAH = findobj(handles, 'tag', 'type_pca');
    PCAString = get(typePCAH, 'string');
    PCAVal = get(typePCAH, 'value');
    typePCA = lower(deblank(PCAString(PCAVal, :)));
    
    if (~strcmpi(typePCA, 'reference') && ~strcmpi(typePCA, 'mccar'))
        fusionInfo.setup_analysis.reference = [];
    end
    
    % Type of PCA
    fusionInfo.setup_analysis.type_pca = typePCA;
    
    % ICA algorithm
    algoH = findobj(handles, 'tag', 'algorithm');
    
    algorithm = get(algoH, 'value');
    
    fusionInfo.setup_analysis.algorithm = algorithm;
    
    % number of groups
    numGroups = length(fusionInfo.setup_analysis.dataInfo);
    fusionInfo.setup_analysis.numGroups = numGroups;
    
    
    % number of tasks
    numFeatures = length(fusionInfo.setup_analysis.dataInfo(1).feature);
    fusionInfo.setup_analysis.numFeatures = numFeatures;
    
    
    fusionPrefix = fusionInfo.setup_analysis.prefix;
    numSubjects = fusionInfo.setup_analysis.numSubjects;
    
    
    % Get the number of components
    compHandle = findobj(handles, 'tag', 'numComp');
    
    try
        compStr = get(compHandle, 'string');
        numComp = strread(compStr, '%d');
    catch
        error('Not a valid integer for number of components');
    end
    
    numComp = numComp(:);
    
    if (length(numComp) == 1)
        numComp = repmat(numComp, 1, numFeatures);
    end
    
    if (length(numComp) ~= numFeatures)
        error(['Number of components (', num2str(length(numComp)), ') doesn''t match the number of features (', num2str(numFeatures), ')']);
    end
    
    % Type of ICA
    typeICAH = findobj(handles, 'tag', 'type_ica');
    ICAString = get(typeICAH, 'string');
    ICAVal = get(typeICAH, 'value');
    typeICA = lower(deblank(ICAString(ICAVal, :)));
    
    % Type of ICA
    fusionInfo.setup_analysis.type_ica = typeICA;
    
    % Number of runs
    numRunsH = findobj(handles, 'tag', 'num_ica_runs');
    
    % number of IC
    fusionInfo.setup_analysis.numComp = numComp;
    
    appName = [fusionPrefix, 'setupAppData'];
    
    if numComp > sum(numSubjects)
        error(['Number of components (', num2str(numComp), ...
            ') is greater than the number of data-sets (', num2str(sum(numSubjects)), ')']);
    end
    
    scaleCompH = findobj(handles, 'tag', 'scale_components');
    scaleCompVal = get(scaleCompH, 'value');
    scaleCompStr = char(get(scaleCompH, 'string'));
    fusionInfo.setup_analysis.scale_components = deblank(scaleCompStr(scaleCompVal, :));
    
    drawnow;
    
    %% Get required information from application name
    if (~isappdata(0, appName))
        stackInfo = ica_fuse_pre_stack_data(fusionInfo);
    else
        stackInfo = getappdata(0, appName);
        rmappdata(0, appName);
    end
    
    set(handles, 'pointer', 'arrow');
    
    handle_visibility = 'off';
    
    if isempty(inputFile)
        handle_visibility = 'on';
    end
    
    % Store dimensions
    newDims = stackInfo.dims;
    
    %% Check EEG modality
    newDims = ica_fuse_compute_new_dims(newDims, cellstr(str2mat(fusionInfo.setup_analysis.dataInfo(1).feature.modality)), stackInfo.mask_ind);
    
    % Sum of sample sizes of features
    spatialDim = sum(newDims);
    
    % Store some fields
    fusionInfo.setup_analysis.mask_ind = stackInfo.mask_ind;
    fusionInfo.setup_analysis.dims = stackInfo.dims;
    fusionInfo.setup_analysis.voxels = stackInfo.voxels;
    fusionInfo.setup_analysis.newDims = newDims;
    
    clear stackInfo;
    
    
    ICA_Options = getICAOptions;
    
    if (~isempty(ICA_Options))
        fusionInfo.setup_analysis.algorithm_opts = ICA_Options;
    end
    
    
    fusionInfo.run_analysis.isInitialized = 0;
    
    % save information to a file
    fusionInfo.setup_analysis.fusionFile = [fusionInfo.setup_analysis.prefix, '_anyway_fusion.mat'];
    fusionFile = fullfile(fusionInfo.setup_analysis.outputDir, [fusionInfo.setup_analysis.fusionFile]);
    
    % save the information for running the analysis
    ica_fuse_save(fusionFile, 'fusionInfo');
    
    disp(['Fusion information is saved in file: ', fusionFile]);
    
    stackInfo.fusionFile = fusionFile;
    
    % set application data
    setappdata(0, appName, stackInfo);
    
    clear stackInfo;
    
    fprintf('\n');
    disp('Please run the analysis using the same fusion information file');
    
    delete(handles);
    
    global APPNAME;
    APPNAME = appName;
    
catch
    
    if exist('handles', 'var')
        set(handles, 'pointer', 'arrow');
    end
    %rethrow(lasterror);
    ica_fuse_displayErrorMsg;
    
end


function maskFiles = getMaskFiles(featureNames, modalities, dataInfo)
% Get mask for each feature if possible

[maskFiles] = ica_fuse_selectMask(featureNames, modalities, dataInfo);


function typePCACallback(hObject, event_data, handles)
% Type of PCA callback

getStr = get(hObject, 'string');
getVal = get(hObject, 'value');

handles_data = get(handles, 'userdata');
inputFile = handles_data.inputFile;
if ~isfield(handles_data.fusionInfo.setup_analysis, 'dataInfo')
    error('Please select the data before selecting PCA type');
end

dataInfo = handles_data.fusionInfo.setup_analysis.dataInfo;
groupNames = str2mat(dataInfo.name);
numSubjects = handles_data.fusionInfo.setup_analysis.numSubjects;
numFeatures = handles_data.fusionInfo.setup_analysis.numFeatures;

selectedStr = deblank(getStr(getVal, :));

reference = [];
if (strcmpi(selectedStr, 'reference') || strcmpi(selectedStr, 'mccar'))
    oldDir = pwd;
    if isempty(inputFile)
        % Open GUI
        if  strcmpi(selectedStr, 'mccar')
            refFile = ica_fuse_selectEntry('title', 'Select reference file for mCCAR', 'filter', '*.asc;*.dat', 'typeEntity', 'file');
            reference = ica_fuse_load_ascii_or_mat(refFile);
            if (length(reference) ~= sum(numSubjects))
                error('Error:Reference', 'Length of reference vector (%d) is not equal to the total number of subjects (%d)', length(reference), sum(numSubjects));
            end
        else
            reference = ica_fuse_get_reference_vector(groupNames, numSubjects);
        end
    else
        % Read from file
        keywd = 'reference';
        try
            inputData = ica_fuse_read_variables(inputFile, keywd, {'double'});
            reference = getfield(inputData, keywd);
        catch
        end
        
        if ~isempty(reference)
            if (numel(reference) == length(reference))
                reference = reference(:);
                if (length(reference) ~= sum(numSubjects))
                    error('Error:Reference', 'Length of reference vector (%d) is not equal to the total number of subjects (%d)', length(reference), sum(numSubjects));
                end
            end
        end
        
    end
    cd(oldDir);
    handles_data.fusionInfo.setup_analysis.reference = reference;
end

%% CCA callback
if (strcmpi(selectedStr, 'cca') || strcmpi(selectedStr, 'mcca'))
    
    if (strcmpi(selectedStr, 'cca') && (numFeatures ~= 2))
        error('CCA can be used only with two features');
    end
    
    if (strcmpi(selectedStr, 'mcca') && (numFeatures == 1))
        error('Please select more than one feature if you want to use MCCA');
    end
    
    if (isempty(inputFile))
        
        %% Open dialog box
        if (isfield(handles_data.fusionInfo.setup_analysis, 'cca_opts'))
            cca_opts = handles_data.fusionInfo.setup_analysis.cca_opts;
        else
            %cca_opts.numPC1 = min(numSubjects);
            %cca_opts.numPC2 = min(numSubjects);
            cca_opts.numPC = repmat(min(numSubjects), 1, handles_data.fusionInfo.setup_analysis.numFeatures);
        end
        
        if (~isfield(cca_opts, 'numPC'))
            cca_opts.numPC = [cca_opts.numPC1, cca_opts.numPC2];
        end
        
        % dialog Title
        dlg_title = 'CCA Options';
        
        numParameters = 1;
        
        % numPC1
        inputText(numParameters).promptString = 'Enter principal components in a vector for features';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(cca_opts.numPC);
        inputText(numParameters).answerType = 'numeric';
        inputText(numParameters).tag = 'numPC';
        inputText(numParameters).enable = 'on';
        inputText(numParameters).value = 0;
        
        %         numParameters = numParameters + 1;
        %
        %         % numPC2
        %         inputText(numParameters).promptString = 'Enter No. Of Principal Components For Feature 2';
        %         inputText(numParameters).uiType = 'edit';
        %         inputText(numParameters).answerString = num2str(cca_opts.numPC2);
        %         inputText(numParameters).answerType = 'numeric';
        %         inputText(numParameters).tag = 'numPC2';
        %         inputText(numParameters).enable = 'on';
        %         inputText(numParameters).value = 0;
        
        %% Get answers
        answers = ica_fuse_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');
        clear inputText;
        
        if (~isempty(answers))
            if (numel(find(cellfun('isempty', answers) == 1)) > 0)
                error('One or more answers for CCA options is empty or not a valid integer');
            end
            cca_opts.numPC = answers{1};
        end
        
    else
        % Read from file
        keywd = 'cca_opts';
        try
            inputData = ica_fuse_read_variables(inputFile, keywd, {'struct'});
            cca_opts = getfield(inputData, keywd);
        catch
            cca_opts.numPC = repmat(min(numSubjects), 1, numFeatures);
        end
    end
    
    
    if (~isfield(cca_opts, 'numPC'))
        cca_opts.numPC = [cca_opts.numPC1, cca_opts.numPC2];
    end
    
    if any(cca_opts.numPC > sum(numSubjects))
        error(['Number of components selected for one of the features is greater than no. of data-sets (', num2str(sum(numSubjects)), ')']);
    end
    
    handles_data.fusionInfo.setup_analysis.cca_opts = cca_opts;
    
end

set(handles, 'userdata', handles_data);


function typeICACallback(handleObj, event_data, handles)
% Get ICA analysis type

handles_data = get(handles, 'userdata');
inputFile = handles_data.inputFile;

if (~isempty(inputFile))
    type_ica = 'average';
    try
        keywd = 'type_ica';
        inputData = ica_fuse_read_variables(inputFile, keywd, {'character'});
        type_ica = getfield(inputData, keywd);
    catch
    end
    inds = strmatch(type_ica, char('average', 'icasso'), 'exact');
    if (~isempty(inds))
        set(handleObj, 'value', inds);
    end
end

function numRunsCallback(handleObj, event_data, handles)
% Get ICA analysis type

handles_data = get(handles, 'userdata');
inputFile = handles_data.inputFile;

if (~isempty(inputFile))
    try
        keywd = 'num_ica_runs';
        inputData = ica_fuse_read_variables(inputFile, keywd, {'integer'});
        num_ica_runs = getfield(inputData, keywd);
        set(handleObj, 'string', num2str(num_ica_runs));
    catch
    end
end


function showHelpDialog(hObject, event_data, handles)
% Show help dialog box

inputTag = get(hObject, 'tag');

inputTag = strrep(inputTag, 'help_', '');

switch (inputTag)
    
    case 'prefix'
        
        titleStr = 'Prefix';
        D(1).string = 'Enter a valid variable name as application data and files will be stored with this name.';
        
    case 'dataInfo'
        
        titleStr = 'Data';
        D(1).string = 'Select the data for the joint ICA analysis. After the data selection, a drop down box will appear with options as ''Yes'' and ''No''. No means data can be selected again.';
        
    case 'maskFile'
        
        titleStr = 'Mask';
        D(1).string = 'Mask used for the analysis. There are two options like ''Default Mask'' and ''Select Mask''. Explanation of each option is given below:';
        D(length(D) + 1).string = '';
        D(length(D) + 1).string = '1. Default Mask - Default mask uses non-zero and not Nan voxels for MRI data whereas for EEG, indices specified in EEG_DATA_INDICES (See ica_fuse_defaults.m) are used.';
        D(length(D) + 1).string = '2. Select Mask - Mask must be specified for each modality.';
        
    case 'normalize'
        
        titleStr = 'Feature normalization';
        D(1).string = 'Each feature is normalized separately before stacking the data. Options are ''Default'', ''Norm2'', ''Std'' and ''None''. Each option is explained below:';
        D(length(D) + 1).string = '';
        D(length(D) + 1).string = '1. Default - Square root of mean of squared data is used.';
        D(length(D) + 1).string = '2. Norm2 - Norm2 is used.';
        D(length(D) + 1).string = '3. Std - Standard deviation is used.';
        D(length(D) + 1).string = '4. None - No normalization.';
        
    case 'estimate_components'
        
        titleStr = 'Dimensionality Estimation';
        D(1).string = ['For MRI data-sets, estimation code is based on i.i.d sampling and Minimum Description Length (MDL) is applied after estimating the samples.', ...
            'For other modalities, standard MDL estimation is used. Also, CCA based dimensionality estimation is used to determine components based on a pair of features.'];
        
    case 'scale_components'
        
        titleStr = 'Scaling Components';
        D(1).string = 'You have the option to scale components to Z-scores.';
        D(length(D) + 1).string = '1. None - No scaling is done.';
        D(length(D) + 1).string = '2. Z-scores - Components are scaled to Z-scores.';
        
    case 'numComp'
        
        titleStr = 'No. of independent components';
        D(1).string = 'Enter the no. of indepenent components you want to extract from the data.';
        
    case 'type_pca'
        
        titleStr = 'Type of PCA';
        D(1).string = 'There are four options like ''Reference'', ''Standard'', ''CCA'' and ''MCCA''.';
        D(length(D) + 1).string = '';
        D(length(D) + 1).string = '1. Standard - Eigen vectors are obtained using eigen decomposition on the data.';
        D(length(D) + 1).string = '2. Reference - Information from groups is used to extract eigen vectors from the data.';
        D(length(D) + 1).string = '3. CCA - Canonical correlation analysis is used to do data reduction. This is same as CCA + joint ICA.';
        D(length(D) + 1).string = '4. MCCA - Multimodality canonical correlation analysis is used to do data reduction. This is same as MCCA + joint ICA.';
        D(length(D) + 1).string = '3. MCCAR - Reference information is used in Multimodality canonical correlation analysis followed by ICA on the features.';
        
        
    case 'type_ica'
        
        titleStr = 'Type of ICA analysis';
        D(1).string = 'There are two options like ''Average'' and ''ICASSO''.';
        D(length(D) + 1).string = '';
        D(length(D) + 1).string = '1. Average - ICA is run multiple times. Components are averaged across runs.';
        D(length(D) + 1).string = '2. ICASSO - ICA is run multiple times using random initialization. Stable ICA run estimates are used in further calculations. ICASSO information is saved in file *icasso*results.mat.';
        
    case 'num_ica_runs'
        
        titleStr = 'Num of ICA runs';
        D(1).string = 'Number of times you want ICA to be run on the data.';
        
    case 'algorithm'
        
        titleStr = 'ICA Algorithm';
        D(1).string = 'Currently, 12 ICA algorithms are available in the toolbox like Infomax, Fast ICA, Erica, Simbec, Evd, Jade Opac, Amuse, SDD ICA, CCICA, Combi, EBM and ERBM.';
        
    case 'reference'
        titleStr = 'Reference file';
        D(1).string = 'If selected yes, reference anyway fusion analysis will be done. Select a reference file containing a vector whose length equals the number of subjects in the analysis';
    otherwise
        
        titleStr = 'Unknown';
        D(1).string = 'Unknown option specified';
        
end
% End for checking

msgStr = char(D.string);
ica_fuse_dialogBox('title', titleStr, 'textBody', msgStr, 'textType', 'large');

function closeCallback(handleObj, event_data, handles)
% closes the figure window

% Close the current figure
delete(handles);


function referenceICACallback(hObject, event_data, handles)
% Reference ICA callback

options_list = cellstr(get(hObject, 'string'));
options_val = get(hObject, 'value');
sel_option = deblank(options_list{options_val});

handles_data = get(handles, 'userdata');
fusionInfo = handles_data.fusionInfo;

numSubjects =  sum(fusionInfo.setup_analysis.numSubjects);

reference  = [];
if strcmpi (sel_option, 'yes')
    refFile = ica_fuse_selectEntry('title', 'Select reference file for anyway fusion analysis', 'filter', '*.asc;*.dat;*.txt', 'typeEntity', 'file');
    reference = ica_fuse_load_ascii_or_mat(refFile);
    reference = reference(:);
    if (length(reference) ~= sum(numSubjects))
        error('Error:Reference', 'Length of reference vector (%d) is not equal to the total number of subjects (%d)', length(reference), sum(numSubjects));
    end
end
fusionInfo.setup_analysis.reference_ica = reference;
handles_data.fusionInfo = fusionInfo;
set(handles, 'userdata', handles_data);

function ICA_Options = getICAOptions


numParameters = 1;

% dialog Title
dlg_title = 'Select the Options for the Anyway Fusion algorithm';

% define all the input parameters in a structure
inputText(numParameters).promptString = 'Select max IVA steps for pre-alignment';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(50);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'ivamaxstep';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 0;


numParameters = numParameters + 1;

% define all the input parameters in a structure
inputText(numParameters).promptString = 'Select lambda regularizer';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(0.2);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'iva_lambda';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 0;

% Input dialog box
answer = ica_fuse_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');

% ICA options with flags and the values corresponding to it
ICA_Options = cell(1, 2*length(answer));

for i = 1:length(answer)
    ICA_Options{2*i - 1} = inputText(i).tag;
    ICA_Options{2*i} = answer{i};
end

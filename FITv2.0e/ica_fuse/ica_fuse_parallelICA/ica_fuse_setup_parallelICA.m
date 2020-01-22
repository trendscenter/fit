function fusionFile = ica_fuse_setup_parallelICA(inputFile, parallel_ica_fusion_file)
% function for setting up parameters required for parallel ICA fusion
% parameters required are:
% 1. output directory
% 2. prefix for output files
% 3. data - of structure format
% 4. mask - mask for fmri data
% 6. Number of components
% 7. Num of times ICA will run
%
% Output: parallel_ica_fusion_file - parallel_ica_fusion_file is a MAT file containing the necessary
% information for running parallel ICA fusion

ica_fuse_defaults;
global PARALLEL_ICA_INFO_MAT_FILE;

% Set inputFile variable to be empty
if ~exist('inputFile', 'var')
    inputFile = [];
end

% Set inputFile variable to be empty
if ~exist('parallel_ica_fusion_file', 'var')
    parallel_ica_fusion_file = [];
end

% DEFINE PARAMETERS THAT ARE GOING TO BE PLOTTED.
% IN CASE OF BATCH ANALYSIS PARAMETERS WILL BE READ FROM
% THE FILE AND THE FUNCTION CALLBACKS WILL BE EXECUTED AUTOMATICALLY
inputText = ica_fuse_define_parameters_paraICA(inputFile);


% Store the inputFile information in figure user data.
handles_data.inputFile = inputFile;

% Read Output directory for the analysis
if ~isempty(inputFile)
    
    [inputPath, inputFile, extn] = fileparts(inputFile);
    
    if isempty(inputPath)
        inputPath = pwd;
    end
    
    inputFile = fullfile(inputPath, [inputFile, extn]);
    
    % check if the input file exists or not
    if ~exist(inputFile, 'file')
        error(['Input file: ', inputFile, ' doesn''t exist']);
    end
    
    % add analysis field 'batch' or 'gui'
    handles_data.analysisType = 'batch';
    handle_visibility = 'off';
    
    % read output directory for the analysis
    keywd = 'outputDir';
    inputData = ica_fuse_read_variables(inputFile, keywd, 'directory');
    % get the output directory field
    outputDir = getfield(inputData, keywd);
    clear inputData;
    % set the input text for controls
    
else
    
    handles_data.analysisType = 'gui';
    
    % If the fusion file is not present select the output directory
    if isempty(parallel_ica_fusion_file)
        % select output directory
        outputDir = ica_fuse_selectEntry('typeEntity', 'directory', 'title', ...
            'Select output directory for parallel ICA fusion');
    else
        formStr = [PARALLEL_ICA_INFO_MAT_FILE, '.mat'];
        [outputDir, fileN, file_extn] = fileparts(parallel_ica_fusion_file);
        inputString = [fileN, file_extn];
        matchIndex = ica_fuse_findstr(inputString, formStr);
        prefixToSet = inputString(1:matchIndex(1)-1);
    end
    
    handle_visibility = 'on';
    
end
% end for reading output directory for the analysis


if isempty(outputDir)
    error('Output directory for the parallel ICA fusion analysis is not selected');
end

% change to the resulting output directory
cd(outputDir);


figureTag = 'Setup Parallel ICA Fusion Analysis';
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

% help on setup ICA GUI
fitHelpTitle = uimenu('parent', InputHandle, 'label', 'FIT-Help');
fitHelpMenu = uimenu(fitHelpTitle, 'label', 'Setup Analysis For Parallel ICA', 'callback', ...
    'ica_fuse_openHTMLHelpFile(''fit_paraICA_setup_analysis.htm'');');

for nH = 1:length(inputText)
    helpHandle = findobj(InputHandle, 'tag', ['help_', inputText(nH).tag]);
    set(helpHandle, 'callback', {@showHelpDialog, InputHandle});
end

% for nH = 1:length(inputText)
%     helpHandle = findobj(InputHandle, 'tag', ['help_', inputText(nH).tag]);
%     set(helpHandle, 'callback', ['feval(@showHelpDialog, ''', inputText(nH).tag, ''')']);
% end

okHandle = findobj(InputHandle, 'tag', okTag);
cancelHandle = findobj(InputHandle, 'tag', cancelTag);

% paraICAInfo is the variable that will be used for saving the parameter
% information
paraICAInfo.setup_analysis.outputDir = outputDir;

% store the input text and paraICAInfo in figure data
handles_data.inputText = inputText;
handles_data.paraICAInfo = paraICAInfo;

% set the input figure user data
set(InputHandle, 'userdata', handles_data);

%%%%%%%%%%%%%% SET UP THE OBJECT CALLBACKS IN ORDER %%%%%%%%%%%%%%%%%%%%%%%%%

% Output prefix callback
set(findobj(InputHandle, 'tag', 'prefix'), 'callback', {@prefixCallback, InputHandle});

% Set callback for data sets
set(findobj(InputHandle, 'tag', 'dataInfo'), 'callback', {@dataSelectCallback, InputHandle});

% mask callback
set(findobj(InputHandle, 'tag', 'maskFile'), 'callback', {@selectMaskCallback, InputHandle});

% Dimensionality estimation callback
set(findobj(InputHandle, 'tag', 'estimate_components'), 'callback', {@dimEstCallback, InputHandle});

% number of components callback
set(findobj(InputHandle, 'tag', 'modality1_numComp'), 'callback', {@numCompCallback, InputHandle});

set(findobj(InputHandle, 'tag', 'SNP_numComp'), 'callback', {@numCompCallback, InputHandle});

set(findobj(InputHandle, 'tag', 'type_parallel_ica'), 'callback', {@typeParallelICACallback, InputHandle});

set(findobj(InputHandle, 'tag', 'type_pca'), 'callback', {@typePCACallback, InputHandle});

% Set type of ICA callback
set(findobj(InputHandle, 'tag', 'type_ica'), 'callback', {@typeICACallback, InputHandle});

% Set number of runs callback
set(findobj(InputHandle, 'tag', 'num_ica_runs'), 'callback', {@numRunsCallback, InputHandle});


% done callback
set(okHandle, 'callback', {@applyCallback, InputHandle});

% cancel callback
set(cancelHandle, 'callback', {@closeCallback, InputHandle});


if exist('prefixToSet', 'var')
    prefixH = findobj(InputHandle, 'tag', 'prefix');
    set(prefixH, 'string', prefixToSet);
    prefixCallback(prefixH, [], InputHandle);
end

% Execute type of parallel ICA callback
%typeParallelICACallback(findobj(InputHandle, 'tag', 'type_parallel_ica'), [], InputHandle);

% Execute the object callbacks automatically for batch file analysis
if ~isempty(inputFile)
    
    %%%% EXECUTE OBJECT CALLBACKS %%%%%%
    objectsToExecute = find([inputText.executeCallback] == 1);
    for nInput = 1:length(objectsToExecute)
        currentH = findobj(InputHandle, 'tag', inputText(objectsToExecute(nInput)).tag);
        ica_fuse_executeCallback(currentH);
    end
    %%%%% END FOR EXECUTING OBJECT CALLBACKS %%%%%
    
    try
        %%%%% EXECUTE DONE CALLBACK %%%%%%%%%
        applyCallback(okHandle, [], InputHandle);
    catch
        if exist('InputHandle', 'var')
            if ishandle(InputHandle)
                delete(InputHandle);
            end
        end
        %rethrow(lasterror);
        ica_fuse_displayErrorMsg;
    end
    
else
    
    % GUI Analysis
    try
        % Wait for the input from the user in case of GUI
        set(InputHandle, 'visible', 'on');
        waitfor(InputHandle);
    catch
        
    end
    
end
% end for executing object callbacks in case of batch file analysis


fusionFile = [];
if isappdata(0, 'paraICA_APPData')
    fusionFile = getappdata(0, 'paraICA_APPData');
    rmappdata(0, 'paraICA_APPData');
end



%%%%%%%%%%%%%%%% DEFINE FUNCTION CALLBACKS %%%%%%%%%%%%%%%%%%%%

% OUTPUT PREFIX CALLBACK
function prefixCallback(handleObj, event_data, handles)

ica_fuse_defaults;
global PARALLEL_ICA_INFO_MAT_FILE;

output_prefix = get(handleObj, 'string');

tagObj = get(handleObj, 'tag');

% check the output prefix
[status, message] = ica_fuse_errorCheck(output_prefix, 'output_prefix');

if status == 0
    error(message);
end

handles_data = get(handles, 'userdata');
handles_data.paraICAInfo.setup_analysis = setfield(handles_data.paraICAInfo.setup_analysis, tagObj, output_prefix);

% set the figure data
set(handles, 'userdata', handles_data);

% update the selecting data files control
% check if the subject file exists or not
% and then update the controls
if isempty(handles_data.inputFile)
    oldDir = handles_data.paraICAInfo.setup_analysis.outputDir;
    dataFile = [handles_data.paraICAInfo.setup_analysis.prefix, PARALLEL_ICA_INFO_MAT_FILE, '.mat'];
    dataFile = fullfile(oldDir, dataFile);
    dataTag = findobj(handles, 'tag', 'dataInfo');
    if exist(dataFile, 'file')
        load(dataFile);
        paraICAInfo.setup_analysis.outputDir = oldDir;
        if isfield(paraICAInfo.setup_analysis, 'dataInfo')
            set(dataTag, 'style', 'popup', 'string', char('Yes', 'No'), 'value', 1);
            handles_data.paraICAInfo = paraICAInfo;
            % set the figure data
            set(handles, 'userdata', handles_data);
            dataSelectCallback(dataTag, [], handles);
        end
    end
end


function dataSelectCallback(handleObj, event_data, handles)
% select the data

set(handles, 'pointer', 'watch');
try
    % Load defaults
    ica_fuse_defaults;
    
    global PARALLEL_ICA_INFO_MAT_FILE;
    
    global PARALLEL_ICA_SEL_DATA_TXT_FILE;
    
    % get the structure for the figure
    handles_data = get(handles, 'userdata');
    paraICAInfo = handles_data.paraICAInfo;
    
    % input file
    inputFile = handles_data.inputFile;
    
    % analysis type
    analysisType = handles_data.analysisType;
    
    % get the inputText
    inputText = handles_data.inputText;
    
    allTags = char(inputText.tag);
    
    % disable all the controls in the input figure
    ica_fuse_enable_control(allTags, handles, 'inactive');
    
    % get the directory for the analysis
    oldDir = paraICAInfo.setup_analysis.outputDir;
    
    % get the prefix of the output string
    paraICAInfo.setup_analysis.prefix = get(findobj(handles, 'tag', 'prefix'), 'string');
    
    % All the entered parameters are going to be stored in this file
    paraICAInfo.setup_analysis.fusionFile = [paraICAInfo.setup_analysis.prefix, PARALLEL_ICA_INFO_MAT_FILE, '.mat'];
    newParamFile = fullfile(oldDir, paraICAInfo.setup_analysis.fusionFile);
    
    fusionPrefix = paraICAInfo.setup_analysis.prefix;
    
    % get the style of the UIcontrol
    getStyle = get(handleObj, 'style');
    
    if strcmp(lower(getStyle), 'pushbutton')
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
    
    set(handles, 'pointer', 'watch');
    
    if ~strcmp(getSubjects, 'yes')
        % Open a input dialog with the following parameters
        
        % selecting data
        [dataInfo, numSubjects] = ica_fuse_select_parallelICA_data(inputFile);
        
        % Modalities files
        for nn = 1:length(dataInfo)
            modality1_files(nn).name = char(dataInfo(nn).feature(1).files.name);
            modality2_files(nn).name = char(dataInfo(nn).feature(2).files.name);
        end
        
        
        paraICAInfo.setup_analysis.dataInfo = dataInfo;
        paraICAInfo.setup_analysis.numSubjects = numSubjects;
        
        %%%%%%%%%%%%%%%% Printing information to a file %%%%%%%%%%%%%
        file_name = fullfile(oldDir, [paraICAInfo.setup_analysis.prefix, PARALLEL_ICA_SEL_DATA_TXT_FILE, '.txt']);
        
        disp(['Saving data information in file ', file_name]);
        fid = fopen(file_name, 'w+');
        try
            ica_fuse_printString(fid, ['Selected ', dataInfo(1).feature(1).modality, ' data files are as follows:'], ...
                char(modality1_files.name));
            fprintf(fid, '\n');
            ica_fuse_printString(fid, ['Selected ', dataInfo(1).feature(2).modality, ' data files are as follows:'], ...
                char(modality2_files.name));
            fclose(fid);
        catch
            fclose(fid);
        end
        
        disp(['Done saving data information']);
        
        % save the information for running the analysis
        ica_fuse_save(newParamFile, 'paraICAInfo');
        
    end
    % end for selecting the data
    
    dataInfo = paraICAInfo.setup_analysis.dataInfo;
    numGroups = length(dataInfo);
    numFeatures = length(dataInfo(1).feature);
    if isfield(dataInfo, 'numSubjects')
        numSubjects = dataInfo.numSubjects;
    else
        numSubjects = zeros(1, length(dataInfo));
        for nn = 1:numGroups
            tempFiles = char(dataInfo(nn).feature(1).files.name);
            timePoints = ica_fuse_get_countTimePoints(tempFiles);
            numSubjects(nn) = timePoints;
        end
    end
    
    featureNames = cellstr(char(dataInfo(1).feature.name));
    
    % Prompts for PC's
    comp1PromptH = findobj(handles, 'tag', ['prompt', 'modality1_numComp']);
    comp2PromptH = findobj(handles, 'tag', ['prompt', 'modality2_numComp']);
    
    % Set strings
    if ~isempty(comp1PromptH)
        set(comp1PromptH, 'string', ['Number of PC for feature ', featureNames{1}]);
    end
    
    % Set strings
    if ~isempty(comp2PromptH)
        set(comp2PromptH, 'string', ['Number of PC for feature ', featureNames{2}]);
    end
    
    paraICAInfo.setup_analysis.numGroups = numGroups;
    paraICAInfo.setup_analysis.numFeatures = numFeatures;
    paraICAInfo.setup_analysis.numSubjects = numSubjects;
    
    % Set figure data
    cd(paraICAInfo.setup_analysis.outputDir);
    
    handles_data.paraICAInfo = paraICAInfo;
    set(handles, 'userdata', handles_data);
    
    if isempty(inputFile)
        % show selected values for controls
        show_values_controls(handles);
    end
    
    set(handles, 'pointer', 'arrow');
catch
    
    set(handles, 'pointer', 'arrow');
    % enable all the controls
    ica_fuse_enable_control(allTags, handles, 'on');
    
    %rethrow(lasterror);
    ica_fuse_displayErrorMsg;
    
end


function show_values_controls(handles)

% get handles data
handles_data = get(handles, 'userdata');

allTags = char(handles_data.inputText.tag);

% get fusion info
paraICAInfo = handles_data.paraICAInfo;

% Number of subjects
numSubjects = paraICAInfo.setup_analysis.numSubjects;

% Set automatically the number of components
if min(numSubjects) < 8
    num_comp = min(numSubjects);
else
    num_comp = 8;
end

numCompH = findobj(handles, 'tag', 'numComp');

compStr = get(numCompH, 'string');

if ~isfield(paraICAInfo.setup_analysis, 'numComp')
    numComp = [];
    try
        numComp = [numComp, paraICAInfo.setup_analysis.modality1_numComp];
    catch
    end
    
    try
        numComp = [numComp, paraICAInfo.setup_analysis.modality2_numComp];
    catch
    end
    
else
    numComp = paraICAInfo.setup_analysis.numComp;
end

if (isempty(numComp))
    if ~isempty(compStr)
        numComp = str2num(compStr);
    else
        numComp = ones(1, paraICAInfo.setup_analysis.numFeatures)*num_comp;
    end
end

if (length(numComp) < paraICAInfo.setup_analysis.numFeatures)
    tmpNumComp = ones(1, paraICAInfo.setup_analysis.numFeatures)*num_comp;
    tmpNumComp(1:length(numComp)) = numComp;
    numComp = tmpNumComp;
else
    numComp = numComp(1:paraICAInfo.setup_analysis.numFeatures);
end

paraICAInfo.setup_analysis.numComp = numComp;

set(numCompH, 'string', num2str(paraICAInfo.setup_analysis.numComp));



% modality1CompH = findobj(handles, 'tag', 'modality1_numComp');
%
% compStr = get(modality1CompH, 'string');
%
% if ~isfield(paraICAInfo.setup_analysis, 'modality1_numComp')
%     if ~isempty(compStr)
%         paraICAInfo.setup_analysis.modality1_numComp = str2num(compStr);
%     else
%         paraICAInfo.setup_analysis.modality1_numComp = num_comp;
%     end
% end
%
% set(modality1CompH, 'enable', 'on');
%
% set(modality1CompH, 'string', num2str(paraICAInfo.setup_analysis.modality1_numComp));
%
%
% modality2CompH = findobj(handles, 'tag', 'modality2_numComp');
%
% compStr = get(modality2CompH, 'string');
%
% if ~isfield(paraICAInfo.setup_analysis, 'modality2_numComp')
%     if ~isempty(compStr)
%         paraICAInfo.setup_analysis.modality2_numComp = str2num(compStr);
%     else
%         paraICAInfo.setup_analysis.modality2_numComp = num_comp;
%     end
% end
%
% set(modality2CompH, 'enable', 'on');
% set(modality2CompH, 'string', num2str(paraICAInfo.setup_analysis.modality2_numComp));

% Number of ICA runs
numICARunsH = findobj(handles, 'tag', 'num_ica_runs');
if isfield(paraICAInfo.setup_analysis, 'numICARuns')
    numICARuns = paraICAInfo.setup_analysis.numICARuns;
    set(numICARunsH, 'enable', 'on');
    set(numICARunsH, 'string', num2str(numICARuns));
end

% set mask information also
maskTag = 'maskFile';
if isfield(paraICAInfo.setup_analysis, maskTag)
    maskHandle = findobj(handles, 'tag', maskTag);
    maskVal = getfield(paraICAInfo.setup_analysis, maskTag);
    if ~isempty(maskVal)
        maskVal = 2;
    else
        maskVal = 1;
    end
    set(maskHandle, 'value', maskVal);
end

set(findobj(handles, 'tag', 'dataInfo'), 'style', 'popup', 'string', ...
    char('Yes', 'No'), 'value', 1);


% Select type of parallel ICA
parallelICAH = findobj(handles, 'tag', 'type_parallel_ica');
parallelICAStrings = lower(get(parallelICAH, 'string'));
if isfield(paraICAInfo.setup_analysis, 'type_parallel_ica')
    parallelICAStr = paraICAInfo.setup_analysis.type_parallel_ica;
    parallelICAVal = strmatch(lower(parallelICAStr), parallelICAStrings, 'exact');
    set(parallelICAH, 'value', parallelICAVal);
else
    parallelICAVal = get(parallelICAH, 'value');
    parallelICAStr = deblank(parallelICAStrings(parallelICAVal, :));
    paraICAInfo.setup_analysis.type_parallel_ica = parallelICAStr;
end


if isfield(paraICAInfo.setup_analysis, 'type_pca')
    type_pca = paraICAInfo.setup_analysis.type_pca;
else
    type_pca = 'reference';
    paraICAInfo.setup_analysis.type_pca = type_pca;
end

% Type of PCA
typePCAH = findobj(handles, 'tag', 'type_pca');
PCAString = get(typePCAH, 'string');

matchIndex = strmatch(type_pca, lower(PCAString), 'exact');

if strcmpi(parallelICAStr, 'aa') ||  strcmpi(parallelICAStr, 'aa-ref')
    set(typePCAH, 'value', matchIndex);
end


if isfield(paraICAInfo.setup_analysis, 'type_ica')
    type_ica = paraICAInfo.setup_analysis.type_ica;
else
    type_ica = 'average';
    paraICAInfo.setup_analysis.type_ica = type_ica;
end

% Type of ICA
typeICAH = findobj(handles, 'tag', 'type_ica');
ICAString = get(typeICAH, 'string');

matchIndex = strmatch(type_ica, lower(ICAString), 'exact');

set(typeICAH, 'value', matchIndex);

% enable all the controls
ica_fuse_enable_control(allTags, handles, 'on');

typeParallelICACallback(parallelICAH, [], handles);

handles_data.paraICAInfo = paraICAInfo;
set(handles, 'userdata', handles_data);

function selectMaskCallback(hObject, event_data, handles)
% Select mask callback

% Load defaults
ica_fuse_defaults;

global MRI_DATA_FILTER;

% select mask callback
handles_data = get(handles, 'userdata');
paraICAInfo = handles_data.paraICAInfo;
inputFile = handles_data.inputFile;

maskTag = get(hObject, 'tag');
getString = get(hObject, 'string');
getMask = get(hObject, 'value');

numFeatures = length(paraICAInfo.setup_analysis.dataInfo(1).feature);
featureNames = cellstr(char(paraICAInfo.setup_analysis.dataInfo(1).feature.name));
modalities = cellstr(char(paraICAInfo.setup_analysis.dataInfo(1).feature.modality));

maskFile = [];

% select mask in 3D Analyze data
if getMask == 2
    if isempty(inputFile)
        % Open GUI for selecting the custom mask
        [P] = getMaskFiles(featureNames, modalities, paraICAInfo.setup_analysis.dataInfo(1));
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
            if strcmpi(modalities{nn}, 'eeg') | strcmpi(modalities{nn}, 'gene')
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
paraICAInfo.setup_analysis = setfield(paraICAInfo.setup_analysis, maskTag, maskFile);

handles_data.paraICAInfo = paraICAInfo;
set(handles, 'userdata', handles_data);

function numCompCallback(hObject, event_data, handles)

% get handles data
handles_data = get(handles, 'userdata');
paraICAInfo = handles_data.paraICAInfo;
compTag = get(hObject, 'tag');
compStr = get(hObject, 'string');
try
    % check if number of components is an integer
    numComp = strread(compStr, '%d');
catch
    error('Not a valid integer for number of components');
end

if numComp < 2
    error('Number of components must be atleast 2 to run ICA');
end

if isfield(paraICAInfo.setup_analysis, 'numSubjects')
    numSubjects = paraICAInfo.setup_analysis.numSubjects;
    if numComp > numSubjects
        error(['Number of components (', num2str(numComp), ...
            ') is greater than the number of data-sets (', num2str(numSubjects), ')']);
    end
end

paraICAInfo.setup_analysis = setfield(paraICAInfo.setup_analysis, compTag, numComp);

handles_data.paraICAInfo = paraICAInfo;
set(handles, 'userdata', handles_data);

function applyCallback(hObject, event_data, handles)
% apply callback

ica_fuse_defaults;
global PARALLEL_ICA_INFO_MAT_FILE;

try
    
    set(handles, 'pointer', 'watch');
    
    handles_data = get(handles, 'userdata');
    
    inputFile = handles_data.inputFile;
    paraICAInfo = handles_data.paraICAInfo;
    
    handle_visibility = 'on';
    if ~isempty(inputFile)
        handle_visibility = 'off';
    end
    
    % Check the data field first
    if ~isfield(paraICAInfo.setup_analysis, 'dataInfo')
        error(['Data is not selected for the analysis']);
    end
    
    % Check the mask field
    if ~isfield(paraICAInfo.setup_analysis, 'maskFile')
        paraICAInfo.setup_analysis.maskFile = [];
    end
    
    
    % Data preprocessing
    paraICAInfo.setup_analysis.preproc_type = 'none';
    try
        preprocH = findobj(handles, 'tag', 'preproc_type');
        preprocOpts = cellstr (get(preprocH, 'string'));
        preprocVal = get(preprocH, 'value');
        paraICAInfo.setup_analysis.preproc_type = preprocOpts{preprocVal};
    catch
    end
    
    % Get the number of components
    compHandle = findobj(handles, 'tag', 'numComp');
    
    try
        compStr = get(compHandle, 'string');
        numComp = str2num(compStr);
        numComp = numComp(:)';
    catch
        error('Not a valid integer for number of components');
    end
    
    paraICAInfo.setup_analysis.numComp = numComp;
    
    %     % Get the number of components
    %     compHandle = findobj(handles, 'tag', 'modality1_numComp');
    %
    %     try
    %         compStr = get(compHandle, 'string');
    %         modality1_numComp = strread(compStr, '%d');
    %     catch
    %         error('Not a valid integer for number of components');
    %     end
    %
    %     % number of modality1 IC
    %     paraICAInfo.setup_analysis.modality1_numComp = modality1_numComp;
    %
    %     % Get the number of components
    %     compHandle = findobj(handles, 'tag', 'modality2_numComp');
    %
    %     try
    %         compStr = get(compHandle, 'string');
    %         modality2_numComp = strread(compStr, '%d');
    %     catch
    %         error('Not a valid integer for number of components');
    %     end
    %
    %     % number of modality1 IC
    %     paraICAInfo.setup_analysis.modality2_numComp = modality2_numComp;
    
    modalities = cellstr(char(paraICAInfo.setup_analysis.dataInfo(1).feature.modality));
    
    % Select type of parallel ICA
    parallelICAH = findobj(handles, 'tag', 'type_parallel_ica');
    getVal = get(parallelICAH, 'value');
    getString = get(parallelICAH, 'string');
    type_parallel_ica = deblank(getString(getVal, :));
    paraICAInfo.setup_analysis.type_parallel_ica = type_parallel_ica;
    
    if strcmpi(type_parallel_ica, 'aa-ref')
        if ~strcmpi(modalities{2}, 'gene');
            error('In order to use AA-ref algorithm, second modality must be gene');
        end
        
        if (~isempty(inputFile))
            try
                keywd = 'snp_ref_file';
                inputData = ica_fuse_read_variables(inputFile, keywd, {'file'});
                snp_ref_file = getfield(inputData, keywd);
            catch
            end
        end
        
        if (~exist('snp_ref_file', 'var'))
            snp_ref_file = ica_fuse_selectEntry('title', 'Select SNP reference file/files in order to use AA-ref algorithm', 'typeEntity', 'file', ...
                'typeSelection', 'multiple', 'filter', '*.asc;*.dat');
            drawnow;
            if (isempty(snp_ref_file))
                error('SNP reference file/files is/are not selected');
            end
        end
        
        snp_ref = ica_fuse_loadData(snp_ref_file);
        
        snp_ref = reshape(snp_ref(:, 2, :), size(snp_ref, 1), size(snp_ref, 3));
        snp_ref = snp_ref';
        
        
        tmpA = load(paraICAInfo.setup_analysis.dataInfo(1).feature(2).files(1).name);
        if (length(tmpA) ~= size(snp_ref, 2))
            error(['SNP reference vector length (', num2str(size(snp_ref, 2)), ') must match the data vector length (', num2str(length(tmpA)), ')']);
        end
        clear tmpA;
        
        
        maskInds = [];
        try
            maskInds = str2num(paraICAInfo.setup_analysis.maskFile{2});
        catch
        end
        
        if (~isempty(maskInds))
            snp_ref = snp_ref(:, maskInds);
        end
        
    end
    
    % Type of PCA
    typePCAH = findobj(handles, 'tag', 'type_pca');
    PCAString = get(typePCAH, 'string');
    PCAVal = get(typePCAH, 'value');
    typePCA = lower(deblank(PCAString(PCAVal, :)));
    
    % Type of PCA
    paraICAInfo.setup_analysis.type_pca = typePCA;
    
    
    
    % Type of ICA
    typeICAH = findobj(handles, 'tag', 'type_ica');
    ICAString = get(typeICAH, 'string');
    ICAVal = get(typeICAH, 'value');
    typeICA = lower(deblank(ICAString(ICAVal, :)));
    
    % Type of ICA
    paraICAInfo.setup_analysis.type_ica = typeICA;
    
    % Number of runs
    numRunsH = findobj(handles, 'tag', 'num_ica_runs');
    numRuns = str2num(get(numRunsH, 'string'));
    if numRuns < 1
        numRuns = 1;
    end
    
    if strcmpi(typeICA, 'icasso') && numRuns < 2
        error('You need to select atleast 2 runs inorder to run ICASSO');
    end
    
    paraICAInfo.setup_analysis.numICARuns = numRuns;
    
    fusionPrefix = paraICAInfo.setup_analysis.prefix;
    
    groupNames = char(paraICAInfo.setup_analysis.dataInfo.name);
    numSubjects = paraICAInfo.setup_analysis.numSubjects;
    
    for n = 1:length(numComp)
        if numComp(n) > numSubjects
            error(['Number of components (', num2str(numComp(n)), ...
                ') is greater than the number of data-sets (', num2str(numSubjects), ')']);
        end
    end
    
    %     if modality1_numComp > numSubjects
    %         error(['Number of components (', num2str(modality1_numComp), ...
    %             ') is greater than the number of data-sets (', num2str(numSubjects), ')']);
    %     end
    %
    %     if modality2_numComp > numSubjects
    %         error(['Number of components (', num2str(modality2_numComp), ...
    %             ') is greater than the number of data-sets (', num2str(numSubjects), ')']);
    %     end
    
    % [ICA_Options] = ica_fuse_paraICAOptions([min([modality1_numComp, modality2_numComp]), sum(numSubjects)], handle_visibility);
    
    ICA_Options = ica_fuse_paraICAOptions([min(numComp), sum(numSubjects)], handle_visibility, numComp);
    
    if strcmpi(type_parallel_ica, 'aa-ref')
        ICA_Options(end+1:end+4) = {'ref_snp', snp_ref, 'lrate', [1e-3, 1e-3]};
    end
    
    paraICAInfo.setup_analysis.ICA_Options = ICA_Options;
    
    % save information to a file
    paraICAInfo.setup_analysis.fusionFile = [paraICAInfo.setup_analysis.prefix, PARALLEL_ICA_INFO_MAT_FILE, '.mat'];
    fusionFile = fullfile(paraICAInfo.setup_analysis.outputDir, [paraICAInfo.setup_analysis.fusionFile]);
    
    % save the information for running the analysis
    ica_fuse_save(fusionFile, 'paraICAInfo');
    
    disp(['Parallel ICA fusion information is saved in file: ', fusionFile]);
    
    setappdata(0, 'paraICA_APPData', fusionFile);
    
    delete(handles);
    
catch
    
    if exist('handles', 'var')
        if ishandle(handles)
            set(handles, 'pointer', 'arrow');
        end
    end
    %rethrow(lasterror);
    ica_fuse_displayErrorMsg;
    
end


function maskFiles = getMaskFiles(featureNames, modalities, dataInfo)
% Get mask for each feature if possible

ica_fuse_defaults;
global MRI_DATA_FILTER;

[maskFiles] = ica_fuse_selectMask(featureNames, modalities, dataInfo);



function typeParallelICACallback(hObject, event_data, handles)
% Type of Parallel ICA callback


getStr = get(hObject, 'string');
getVal = get(hObject, 'value');

typePCAH = findobj(handles, 'tag', 'type_pca');

if strcmpi(deblank(getStr(getVal, :)), 'aa') || strcmpi(deblank(getStr(getVal, :)), 'aa-ref')
    set(typePCAH, 'enable', 'on');
else
    set(typePCAH, 'enable', 'off');
end

function typePCACallback(hObject, event_data, handles)
% Type of PCA callback

getStr = get(hObject, 'string');
getVal = get(hObject, 'value');

handles_data = get(handles, 'userdata');
inputFile = handles_data.inputFile;
dataInfo = handles_data.paraICAInfo.setup_analysis.dataInfo;
groupNames = char(dataInfo.name);
numSubjects = handles_data.paraICAInfo.setup_analysis.numSubjects;

selectedStr = deblank(getStr(getVal, :));

reference = [];
if strcmpi(selectedStr, 'reference')
    oldDir = pwd;
    if isempty(inputFile)
        % Open GUI
        reference = ica_fuse_get_reference_vector(groupNames, numSubjects);
    else
        % Read from file
        keywd = 'reference';
        try
            inputData = ica_fuse_read_variables(inputFile, keywd, {'double'});
            reference = getfield(inputData, keywd);
        catch
        end
        
        if ~isempty(reference)
            reference = reference(:);
            if (length(reference) ~= sum(numSubjects))
                error('Error:Reference', 'Length of reference vector (%d) is not equal to the total number of subjects (%d)', length(reference), sum(numSubjects));
            end
        end
        
    end
    cd(oldDir);
    handles_data.paraICAInfo.setup_analysis.reference = reference;
end

set(handles, 'userdata', handles_data);


function closeCallback(handleObj, event_data, handles)
% closes the figure window

% Close the current figure
delete(handles);


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
        D(1).string = 'Select the data for the parallel ICA analysis. After the data selection, a drop down box will appear with options as ''Yes'' and ''No''. No means data can be selected again.';
        
    case 'maskFile'
        
        titleStr = 'Mask';
        D(1).string = 'Mask used for the analysis. There are two options like ''Default Mask'' and ''Select Mask''. Explanation of each option is given below:';
        D(length(D) + 1).string = '';
        D(length(D) + 1).string = '1. Default Mask - Default mask uses non-zero and not Nan voxels for MRI data whereas for SNP''s all indices are included.';
        D(length(D) + 1).string = '2. Select Mask - Mask must be specified for each modality. For SNP''s include only the indices you want for the analysis.';
        
    case 'numComp'
        
        titleStr = 'Num of PC for features';
        D(1).string = 'Enter number of principal components to be extracted from the data for all the features in a row vector like 5, 7';
        
    case 'type_parallel_ica'
        
        titleStr = 'Type of parallel ICA';
        D(1).string = 'There are three options like ''AA'', ''AS'' and ''AA-ref''.';
        D(length(D) + 1).string = '';
        D(length(D) + 1).string = '1. AA option uses correlation measure between mixing coefficient of modality 1 with mixing coefficient of modality 2.';
        D(length(D) + 1).string = '2. AS option uses correlation measure between mixing coefficient of modality 1 with source of modality 2.';
        D(length(D) + 1).string = '3. AA-ref option uses AA algorithm with SNP constraint.';
        
    case 'type_pca'
        
        titleStr = 'Type of PCA';
        D(1).string = 'There are two options like ''Reference'' and ''Standard''. Reference option uses information from groups to project eigen vectors to that dimension.';
        
        
    case 'type_ica'
        
        titleStr = 'Type of ICA analysis';
        D(1).string = 'There are two options like ''Average'' and ''ICASSO''.';
        D(length(D) + 1).string = '';
        D(length(D) + 1).string = '1. Average - ICA is run multiple times. Components are averaged across runs.';
        D(length(D) + 1).string = '2. ICASSO - ICA is run multiple times using random initialization. Stable ICA run estimates are used in further calculations. ICASSO information is saved in file *icasso*results.mat.';
        
    case 'num_ica_runs'
        
        titleStr = 'Num of ICA runs';
        D(1).string = 'Number of times you want ICA to be run on the data.';
        
    case 'estimate_components'
        
        titleStr = 'Dimensionality Estimation';
        D(1).string = ['For MRI data-sets, estimation code is based on i.i.d sampling and Minimum Description Length (MDL) is applied after estimating the samples.', ...
            'For other modalities, standard MDL estimation is used. Also, CCA based dimensionality estimation is used to determine components based on a pair of features.'];
        
    case 'preproc_type'
        
        titleStr = 'Preprocessing Type';
        D(1).string = 'Options are none and Z-scores. If you selected z-scores, data for each modality is scaled to z-scores.';
        
        
    otherwise
        
        titleStr = 'Unknown';
        D(1).string = 'Unknown option specified';
        
end
% End for checking

msgStr = char(D.string);
ica_fuse_dialogBox('title', titleStr, 'textBody', msgStr, 'textType', 'large');


function dimEstCallback(hObject, event_data, handles)
%% Dimensionality estimation callback
%

ica_fuse_defaults;
global FIG_FG_COLOR;
global FIG_BG_COLOR;
global AX_COLOR;

% select mask callback
handles_data = get(handles, 'userdata');

inputFile = handles_data.inputFile;

if (isempty(inputFile))
    
    drawnow;
    
    paraICAInfo = handles_data.paraICAInfo;
    
    if ~isfield(paraICAInfo.setup_analysis, 'dataInfo')
        error('Please select data in order to use dimensionality estimation');
    end
    
    maskFile = [];
    if isfield(paraICAInfo.setup_analysis, 'maskFile')
        maskFile = paraICAInfo.setup_analysis.maskFile;
    end
    
    % Create mask
    mask_ind = ica_fuse_createMask(paraICAInfo.setup_analysis.dataInfo, maskFile);
    
    % Get featureInfo
    featureInfo = ica_fuse_get_feature_info(paraICAInfo.setup_analysis.dataInfo);
    
    numFeatures = length(featureInfo);
    
    estimationInfo = repmat(struct('comp', [], 'mdl', [], 'aic', [], 'feature_name', ''), 1, length(featureInfo));
    for nF = 1:length(featureInfo)
        
        disp(['Doing dimensionality estimation for feature ', featureInfo(nF).feature_name, ' ...']);
        
        % Apply mask
        featureData = ica_fuse_applyMask(featureInfo(nF), mask_ind(nF));
        
        % Do iid sampling for fmri, smri data
        doSampling = (strcmpi(featureInfo(nF).modality, 'fmri') || strcmpi(featureInfo(nF).modality, 'smri'));
        
        estimationInfo(nF).feature_name = featureInfo(nF).feature_name;
        
        [estimationInfo(nF).comp, estimationInfo(nF).mdl, estimationInfo(nF).aic] = ica_fuse_estim_dim(featureData.data, 'maskvec', mask_ind(nF).ind, 'doSampling', doSampling);
        
        fprintf('\n');
        
        disp(['Components estimated to be ', num2str(estimationInfo(nF).comp)]);
        
        fprintf('\n');
        
        clear featureData;
        
    end
    
    drawnow;
    
    handles_data.paraICAInfo.setup_analysis.estimationInfo = estimationInfo;
    
    set(handles, 'userdata', handles_data);
    
    fprintf('Done with dimensionality estimation \n');
    
    
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
    
    if (length(featureInfo) > 1)
        
        disp('Using PCA-CCA based estimation ..');
        allCombinations = nchoosek(1:numFeatures, 2);
        allFeatureNames = cellstr(char(featureInfo.feature_name));
        ccaEstimationInfo = repmat(struct('name', [], 'comp', []), 1, size(allCombinations, 1));
        for nF = 1:size(allCombinations, 1)
            featureA = allCombinations(nF, 1);
            featureB = allCombinations(nF, 2);
            % feature 1 and 2
            disp(['Loading features ', ccaEstimationInfo(nF).name, '...']);
            xData = ica_fuse_applyMask(featureInfo(featureA), mask_ind(featureA));
            yData = ica_fuse_applyMask(featureInfo(featureB), mask_ind(featureB));
            ccaEstimationInfo(nF).name = [allFeatureNames{featureA}, ' & ', allFeatureNames{featureB}];
            ccaEstimationInfo(nF).comp = ica_fuse_model_est_cca(xData.data', yData.data');
            
        end
        
        handles_data.paraICAInfo.setup_analysis.ccaEstimationInfo = ccaEstimationInfo;
        dimMsgString = [char(ccaEstimationInfo.name), repmat(' : ', length(ccaEstimationInfo), 1), num2str([ccaEstimationInfo.comp]')];
        dimMsgString = cellstr(dimMsgString);
        textBody = {'PCA-CCA estimation is done on pair of features and dimensionality estimation is shown below:';'';''};
        textBody = [textBody;dimMsgString];
        ica_fuse_dialogBox('title', 'PCA-CCA Estimation', 'textType', 'large', 'textbody', textBody);
        disp(char(textBody));
    end
    
    set(handles, 'userdata', handles_data);
    drawnow;
    
    
end



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
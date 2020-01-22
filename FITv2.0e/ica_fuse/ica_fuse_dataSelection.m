function [groupsData] = ica_fuse_dataSelection(inputFile, isCCA)
% Include in data selection appropriate selection criteria

% load defaults;
ica_fuse_defaults;
global AVAILABLE_MODALITIES;

if ~exist('isCCA', 'var')
    isCCA = 0;
end

% make input file empty if it does not exist
if ~exist('inputFile', 'var')
    inputFile = [];
end

% GET THE NUMBER OF GROUPS AND FEATURES
if isempty(inputFile)
    
    % open a dialog box for getting the number of groups and features
    numParameters = 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Number of groups';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '2';
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'numGroups';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    if (~isCCA)
        
        numParameters = numParameters + 1;
        
        % define all the input parameters in a structure
        inputText(numParameters).promptString = 'Number of features';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '2';
        inputText(numParameters).answerType = 'numeric';
        inputText(numParameters).tag = 'numFeatures';
        inputText(numParameters).enable = 'on';
        inputText(numParameters).value = 1;
        
    end
    
    % Input dialog box
    answer_groups_features = ica_fuse_inputDialog('inputtext', inputText, 'Title', ...
        'Enter no. of groups and features');
    
    clear inputText;
    
    if isempty(answer_groups_features)
        error('Figure window to get the number of groups and features was quit');
    end
    
else
    
    % read number of groups and features
    keywds = {'groupNames', 'featureNames'};
    % read output directory for the analysis
    inputData = ica_fuse_read_variables(inputFile, keywds, {'cell', 'cell'});
    groupNames = getfield(inputData, 'groupNames'); % group names
    featureNames = getfield(inputData, 'featureNames'); % feature names
    
    %%% Force feature names not to contain & %%%%%%
    for nFea = 1:length(featureNames)
        matchFeatureIndex = ica_fuse_findstr(featureNames{nFea}, '&');
        if ~isempty(matchFeatureIndex)
            error(['Feature name for feature ', num2str(nFea), ' should not contain & character.']);
        end
    end
    %%%% end for forcing the feature names not to contain & %%%%
    
    answer_groups_features{1} = length(groupNames); answer_groups_features{2} = length(featureNames);
    clear inputData;
    
end
% END FOR GETTING THE NUMBER OF GROUPS AND FEATURES

% get the number of groups and features
numGroups = answer_groups_features{1};

if (~isCCA)
    numFeatures = answer_groups_features{2};
else
    numFeatures = 2;
end


if isempty(numGroups)
    error('Number of groups is not a valid integer');
end

if isempty(numFeatures)
    error('Number of features is not a valid integer');
end

if numGroups*numFeatures < 2
    error('Need atleast two features or groups to do data fusion');
end

if (isCCA && (numFeatures ~= 2))
    error('CCA + ICA works with only two features');
end

% READ MODALITY, GROUP NAMES, FEATURE NAMES
if ~isempty(inputFile)
    
    % get the modalities
    keywd = 'modality';
    inputData = ica_fuse_read_variables(inputFile, keywd, {'cell'});
    answersInput = getfield(inputData, keywd);
    
    if length(answersInput) ~= numFeatures
        error('Number of modalities must equal the number of features');
    end
    
    if isempty(answersInput)
        error(['Check the modality variable in the file: ', inputFile]);
    end
    
    for nInput = 1:length(answersInput)
        if ~(strcmpi(answersInput{nInput}, 'fmri') || strcmpi(answersInput{nInput}, 'eeg') || strcmpi(answersInput{nInput}, 'smri'))
            error('Error:Modality', ['Unrecognized modality (%s) in file:%s.\nCurrently available  modalities are fmri, eeg and smri.'], ...
                answersInput{nInput}, inputFile);
        end
    end
    
    % get the group names
    for nGroup = 1:numGroups
        inputLength = length(answersInput) + 1;
        answersInput{inputLength} = groupNames{nGroup};
    end
    
    % get the feature names
    for nFeature = 1:numFeatures
        inputLength = length(answersInput) + 1;
        answersInput{inputLength} = featureNames{nFeature};
    end
    
else
    
    inputTextLength = 2*numFeatures + numGroups;
    
    inputText = repmat(struct('promptString', [], 'uiType', 'edit', 'answerString', [], 'answerType', 'string', 'tag', [], 'enable', [], ...
        'value', 1), 1, inputTextLength);
    
    textCount = 0;
    % get the type of features
    for nFeature = 1:numFeatures
        textCount = textCount + 1;
        % define all the input parameters in a structure
        inputText(nFeature).promptString = ['Select modality for feature ', num2str(nFeature), ':'];
        inputText(nFeature).uiType = 'popup';
        inputText(nFeature).answerString = AVAILABLE_MODALITIES;
        inputText(nFeature).answerType = 'string';
        inputText(nFeature).tag = ['modality', num2str(nFeature)];
        inputText(nFeature).enable = 'on';
        inputText(nFeature).value = 1;
        
    end
    % end loop over features
    
    %textCount = length(inputText);
    
    % form input text structure for getting input from the user
    for nGroup = 1:numGroups
        textCount = textCount + 1;
        % define all the input parameters in a structure
        inputText(textCount).promptString = ['Name ', 'Group', ' ', num2str(nGroup), ':'];
        inputText(textCount).uiType = 'edit';
        inputText(textCount).answerString = ['Group', num2str(nGroup)];
        inputText(textCount).answerType = 'string';
        inputText(textCount).tag = ['Group', num2str(nGroup)];
        inputText(textCount).enable = 'on';
        inputText(textCount).value = 1;
    end
    
    %textCount = length(inputText);
    
    % form input text structure for getting input from the user
    for nFeature = 1:numFeatures
        textCount = textCount + 1;
        % define all the input parameters in a structure
        inputText(textCount).promptString = ['Name ', 'Feature', ' ', num2str(nFeature), ':'];
        inputText(textCount).uiType = 'edit';
        inputText(textCount).answerString = ['Feature', num2str(nFeature)];
        inputText(textCount).answerType = 'string';
        inputText(textCount).tag = ['Feature', num2str(nFeature)];
        inputText(textCount).enable = 'on';
        inputText(textCount).value = 1;
    end
    
    
    % Input dialog box
    answersInput = ica_fuse_inputDialog('inputtext', inputText, 'Title', ...
        'Select modality and give namings for groups and features ...');
    
    if isempty(answersInput)
        error('Figure window to get the type of modality, naming for features and groups was quit.');
    end
    
    
    countFeatures = 0;
    %%% Force feature names not to contain & %%%%%%
    for nFea = length(answersInput) - numFeatures + 1:length(answersInput)
        countFeatures = countFeatures + 1;
        matchFeatureIndex = ica_fuse_findstr(answersInput{nFea}, '&');
        if ~isempty(matchFeatureIndex)
            error(['Feature name for feature ', num2str(countFeatures), ' should not contain & character.']);
        end
    end
    %%%% end for forcing the feature names not to contain & %%%%
    
end
% END FOR READING MODALITY, GROUP NAMES, FEATURE NAMES

clear inputText;

allModalities = answersInput(1:numFeatures);
allModalities = repmat(allModalities, 1, numGroups);

% loop over groups
for nGroup = 1:numGroups
    groupsData(nGroup).name = answersInput{nGroup + numFeatures};
    % loop over features
    for nFeature = 1:numFeatures
        groupsData(nGroup).feature(nFeature).modality = answersInput{nFeature};
        groupsData(nGroup).feature(nFeature).name = answersInput{nFeature + numGroups + numFeatures};
    end
    % end loop over features
end
% end loop over groups

clear answersInput;

%% File selection
if (isempty(inputFile))
    dataSelectionOption = ica_fuse_questionDialog('title', 'Data Selection', 'textbody', 'Do you want to use file selection window for selecting files for each group and feature?');
    dataSelectionOption = dataSelectionOption + 1;
    drawnow;
else
    
    dataSelectionOption = 1;
    try
        keywd = 'dataSelectionOption';
        inputData = ica_fuse_read_variables(inputFile, keywd, {'integer'});
        dataSelectionOption = getfield(inputData, keywd);
        clear inputData;
    catch
    end
    
end


%% Use userinput of file names rather than selecting based on file patterns
if (dataSelectionOption == 2)
    
    if isempty(inputFile)
        %% GUI
        
        listBoxStrings = repmat(struct('name', ''), 1, numGroups*numFeatures);
        
        countF = 0;
        for nGroup = 1:numGroups
            for nFeature = 1:numFeatures
                countF = countF + 1;
                listBoxStrings(countF).name = [groupsData(nGroup).name, ' ', groupsData(nGroup).feature(nFeature).name];
            end
        end
        
        fileTypes = repmat({'any'}, 1, numGroups*numFeatures);
        
        isMRI = strcmpi(allModalities, 'fmri') | strcmpi(allModalities, 'smri');
        
        fileTypes(isMRI) = {'image'};
        
        tmp_files = ica_fuse_ui_select_data('title', 'Select data for fusion', 'num_subjects', ...
            numGroups, 'num_sessions', numFeatures, 'filter_string', '*.img;*.nii;*.asc', 'type_file_selection', 'multiple', 'fileType', fileTypes, ...
            'subject_string', listBoxStrings, 'files_specification', 'equal_sessions');
        
        if (isempty(tmp_files))
            error('Files are not selected for fusion');
        end
        
        countF = 0;
        for nGroup = 1:numGroups
            for nFeature = 1:numFeatures
                countF  = countF + 1;
                tmp = tmp_files(countF).name;
                for ii = 1:size(tmp, 1)
                    groupsData(nGroup).feature(nFeature).files(ii).name = deblank(tmp(ii, :));
                end
            end
        end
        
        clear tmp_files;
        
    else
        %% File input
        
        keywd = 'files';
        inputData = ica_fuse_read_variables(inputFile, keywd, {'cell'});
        grp_files = getfield(inputData, keywd);
        clear inputData;
        
        if (isempty(grp_files))
            error('Files are not selected for fusion');
        end
        
        if (size(grp_files, 1) ~= numGroups)
            error(['No. of rows in the files variable must equal the no. of groups(', num2str(numGroups), ')']);
        end
        
        if (size(grp_files, 2) ~= numFeatures)
            error(['No. of cols in the files variable must equal the no. of features(', num2str(numFeatures), ')']);
        end
        
        countF = 0;
        for nGroup = 1:numGroups
            for nFeature = 1:numFeatures
                
                countF = countF + 1;
                
                tmp = deblank(str2mat(grp_files{nGroup, nFeature}));
                
                if isempty(tmp)
                    error(['No files passed for group ', groupsData(nGroup).name, ' feature ', groupsData(nGroup).feature(nFeature).name]);
                end                                
                
                if (size(tmp, 1) == 1)
                    [tmp_path, tmpF, extn] = fileparts(tmp);
                    tmp = ica_fuse_listFiles_inDir(tmp_path, [tmpF, extn]);
                    if isempty(tmp)
                        error(['Files not found for (', fullfile(tmp, [tmpF, extn]), ')']);
                    end
                    tmp = ica_fuse_fullFile('files', tmp, 'directory', tmp_path);
                end
                
                tmp = ica_fuse_rename_4d_file(tmp);
                
                if (nFeature == 1)
                    oldFileCount = size(tmp, 1);
                else
                    if (oldFileCount ~= size(tmp, 1))
                        error(['No. of files is not the same between the features for group ', groupsData(nGroup).name]);
                    end
                end
                
                for ii = 1:size(tmp, 1)
                    groupsData(nGroup).feature(nFeature).files(ii).name = deblank(tmp(ii, :));
                end
                
            end
        end
    end
    
    
else
    
    % GET ANSWER FOR FILE PATTERN
    if numGroups > 1
        if ~isempty(inputFile)
            keywd = 'answerFilePattern';
            inputData = ica_fuse_read_variables(inputFile, keywd, {'integer'});
            answerFilePattern = getfield(inputData, keywd);
            clear inputData;
            % check if it is not 1 or 0
            if answerFilePattern ~= 0 && answerFilePattern ~= 1
                disp(['Answer for file pattern must be 1 or 0. By default setting the answer to 1 ...']);
                answerFilePattern = 1;
            end
        else
            % get the answer for file pattern
            answerFilePattern = ica_fuse_questionDialog('title', 'File pattern', 'textbody', ...
                'Is the file pattern same between groups?');
        end
    else
        answerFilePattern = 1;
    end
    % END FOR GETTING FILE PATTERN
    
    
    % add file pattern to the structure groups data
    [groupsData] = getFilePattern(numGroups, numFeatures, answerFilePattern, groupsData, inputFile);
    
    
    answerDir = 0;
    % Ask question about directory
    if numGroups > 1
        
        if isempty(inputFile)
            % get the answer for file pattern
            answerDir = ica_fuse_questionDialog('title', 'Directories ...', 'textbody', 'Is the data organized in one group folder?');
        else
            % read the variable in inputFile
            keywd = 'answerDir';
            inputData = ica_fuse_read_variables(inputFile, keywd, {'integer'});
            answerDir = getfield(inputData, keywd);
            if answerDir > 1
                answerDir = 1;
            end
            clear inputData;
        end
    end
    % end for asking the question about directory
    
    % If the answerDir is 1 then select the first group folder only
    if answerDir == 1
        num_groups = 1;
    else
        % select the respective group folder
        num_groups = numGroups;
    end
    
    startPath = pwd;
    % select the directories for groups
    for nGroup = 1:num_groups
        
        if isempty(inputFile)
            getDirGroups{nGroup} = ica_fuse_selectEntry('typeEntity', 'directory', 'startPath', startPath, 'title', ...
                ['Select directory for ', groupsData(nGroup).name]);
        else
            keywd = ['group', num2str(nGroup), '_dir'];
            inputData = ica_fuse_read_variables(inputFile, keywd, {'directory'});
            temp = getfield(inputData, keywd);
            getDirGroups{nGroup} = deblank(temp);
            clear temp; clear inputData; clear keywd;
        end
        startPath = getDirGroups{nGroup};
        if isempty(getDirGroups{nGroup})
            error(['Directory for group ', num2str(nGroup), ' is not selected.']);
        end
        
        if ~exist(startPath, 'dir')
            error(['Directory for group ', groupsData(nGroup).name, ' doesn''t exist']);
        end
        
    end
    % END FOR READING DIRECTORY INFORMATION
    
    % replicate the directories over groups if there is one directory for all
    % groups and features
    if length(getDirGroups) == 1
        firstDir = getDirGroups{1}; clear getDirGroups;
        for nGroup = 1:numGroups
            getDirGroups{nGroup} = firstDir;
        end
    end
    % end for replicating over groups
    
    drawnow;
    
    % PERFORM THE FINAL STEP FOR READING FILES
    [groupsData] = getFiles(groupsData, getDirGroups);
    
    
    % Image dimensions for a feature between the groups must be the same
    %
    % Loop over features
    for nF = 1:numFeatures
        
        clear allFiles;
        clear currentFiles;
        
        % Loop over groups
        for nG = 1:numGroups
            currentFiles(nG).name = str2mat(groupsData(nG).feature(nF).files.name);
        end
        % End loop over groups
        
        allFiles = str2mat(currentFiles.name);
        
        try
            % check the image dimensions & no. of images
            [countTimePoints, extns, dims] = ica_fuse_get_countTimePoints(allFiles);
        catch
            errorMsg = lasterr;
            if ica_fuse_findstr(errorMsg, 'Data dimensions should be the same for all files')
                strIndices = ica_fuse_findstr(errorMsg, 'Please check');
                newStr = errorMsg(strIndices(1):end);
                error('Error:Data_dimensions', ['Image dimensions must be the same within a feature (%s).\n%s'], ...
                    groupsData(1).feature(nF).name, newStr);
            else
                ica_fuse_displayErrorMsg;
            end
        end
        
        % Check no. of images between features
        if nF == 1
            nTimePoints = countTimePoints;
        else
            if countTimePoints ~= nTimePoints
                error('Error:Num_images', ['No. of images must be the same between features.\n Check no. of images for feature %s'], ...
                    groupsData(1).feature(nF).name);
            end
        end
        % End for checking no. of images between features
        
    end
    % End loop over features
    
end


%%%%%%%%%%%%%%% DEFINE SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%

function [groupStr] = getNaming(numGroups, dlgTitle)
% open input dialog box to name

% form input text structure for getting input from the user
for ii = 1:numGroups
    % define all the input parameters in a structure
    inputText(ii).promptString = ['Name ', dlgTitle, ' ', num2str(ii), ':'];
    inputText(ii).uiType = 'edit';
    inputText(ii).answerString = [dlgTitle, num2str(ii)];
    inputText(ii).answerType = 'string';
    inputText(ii).tag = [dlgTitle, ' ', num2str(ii)];
    inputText(ii).enable = 'on';
    inputText(ii).value = 1;
end

if ~exist('dlgTitle', 'var')
    dlgTitle = 'Group';
end

% Input dialog box
groupStr = ica_fuse_inputDialog('inputtext', inputText, 'Title', ['Naming ', dlgTitle]);

if isempty(groupStr)
    error(['figure window with title: Naming ', dlgTitle ' is closed.']);
end


function [groupsData] = getFilePattern(numGroups, numFeatures, checkSameFileP, groupsData, inputFile)
% add file pattern to the structure groupsData

ica_fuse_defaults;
global MRI_DATA_FILTER;
global EEG_DATA_FILTER;
global BEH_DATA_FILTER;


if ~isempty(inputFile)
    handle_visibility = 'off';
else
    handle_visibility = 'on';
end

if checkSameFileP == 1
    num_groups = 1;
else
    num_groups = numGroups;
end

% loop over groups
for nGroup = 1:num_groups
    % loop over features
    if  strcmpi(handle_visibility, 'off')
        keywd = ['group', num2str(nGroup), '_file_pattern'];
        % read input
        inputData = ica_fuse_read_variables(inputFile, keywd, {'cell'});
        answer_file_pattern = getfield(inputData, keywd);
        clear keywd; clear inputData;
        if length(answer_file_pattern) ~= numFeatures
            error('Number of file patterns must equal the number of features');
        end
    end
    
    
    for nFeature = 1:numFeatures
        
        if strcmpi( handle_visibility, 'off')
            fileP = answer_file_pattern{nFeature};
        else
            
            if strcmpi(groupsData(nGroup).feature(nFeature).modality, 'fmri') || ...
                    strcmpi(groupsData(nGroup).feature(nFeature).modality, 'smri')
                fileP = MRI_DATA_FILTER;
            elseif strcmpi(groupsData(nGroup).feature(nFeature).modality, 'behavioral')
                fileP = BEH_DATA_FILTER;
            else
                fileP = EEG_DATA_FILTER;
            end
        end
        % define all the input parameters in a structure
        inputText(nFeature).promptString = ['Select file pattern for ', groupsData(nGroup).name, ' ', ...
            groupsData(nGroup).feature(nFeature).name, ' feature'];
        inputText(nFeature).uiType = 'edit';
        inputText(nFeature).answerString = fileP;
        inputText(nFeature).answerType = 'string';
        inputText(nFeature).tag = ['feature', num2str(nFeature)];
        inputText(nFeature).enable = 'on';
        inputText(nFeature).value = 1;
    end
    % end loop over features
    
    clear answer_file_pattern;
    
    % Input dialog box
    answer_file_pattern = ica_fuse_inputDialog('inputtext', inputText, 'Title', ...
        ['Select file pattern for group ', num2str(nGroup)], 'handle_visibility', handle_visibility);
    
    if isempty(answer_file_pattern)
        if strcmpi( handle_visibility, 'off')
            error(['error in reading file pattern for group ', num2str(nGroup)]);
        else
            error(['file pattern for group ', num2str(nGroup), ' input dialog box was quit']);
        end
    end
    
    % loop over number of features
    for nFeature = 1:numFeatures
        currentP = answer_file_pattern{nFeature};
        if strcmpi(groupsData(nGroup).feature(nFeature).modality, 'fmri') || ...
                strcmpi(groupsData(nGroup).feature(nFeature).modality, 'smri')
            
            if ~(strcmpi(currentP(end-3:end), '.img') || strcmpi(currentP(end-3:end), '.nii') )
                error(['Image format for group ', groupsData(nGroup).name, ' feature ',  ...
                    groupsData(nGroup).feature(nFeature).name, ' must be .img or .nii']);
            end
        end
    end
    % end for loop over number of features
    
    answerFileP(nGroup).name = answer_file_pattern;
    
    clear inputText;
    clear answer_file_pattern;
    
end
% end loop over groups

% replicate the file pattern over groups
if num_groups == 1
    answerFileP = repmat(answerFileP, numGroups, 1);
end

% replicate the file pattern over groups
for nGroup = 1:numGroups
    
    % loop over features
    for nFeature = 1:numFeatures
        groupsData(nGroup).feature(nFeature).filePattern = answerFileP(nGroup).name{nFeature};
    end
    % end loop over features
    
end
% end loop over groups


function [groupsData] = getFiles(groupsData, getDirGroups)
% search files in the directory

% Loop over groups
for nGroup = 1:length(groupsData)
    
    % get the number of features
    numFeatures = length(groupsData(nGroup).feature);
    currentDir = getDirGroups{nGroup};
    
    % Loop over features
    for nFeature = 1:numFeatures
        
        fileP = groupsData(nGroup).feature(nFeature).filePattern;
        [selFiles] = ica_fuse_listFiles_inDir(currentDir, fileP);
        
        % get the selected files
        if ~isempty(selFiles)
            % loop over the selected files
            for ii = 1:size(selFiles, 1)
                selectedFiles(ii).name = fullfile(currentDir, deblank(selFiles(ii, :)));
            end
            
        else
            
            % First read the data in sub-folders and then sub-sub folders
            try
                disp(['Reading data from sub folders in directory ', currentDir]);
                selectedFiles = ica_fuse_get_sub_data(currentDir, fileP, 'data_in_subject_folder', 0);
            catch
                % try for sub-sub folders
                try
                    disp(['Reading data from sub sub-folders in directory ', currentDir]);
                    selectedFiles = ica_fuse_get_sub_data(currentDir, fileP, 'data_in_subject_subfolder', 0);
                catch
                    error('Error:Data', ['Data doesn''t exist with the file pattern (%s) for group %s and feature %s'], ...
                        fileP, groupsData(nGroup).name, groupsData(nGroup).feature(nFeature).name);
                end
                % end for try for sub-sub folders
                
            end
            % end for try for sub folders
            
        end
        % end for getting the selected files
        
        % Loop over the selected files
        for ii = 1:length(selectedFiles)
            currentFile = deblank(selectedFiles(ii).name(1, :));
            
            % TEST IF THE SELECTED EEG FILES ARE ASCII OR NOT
            if strcmpi(groupsData(nGroup).feature(nFeature).modality, 'eeg')
                ica_fuse_check_ascii_file(currentFile);
            else
                % check the header file for analyze images
                ica_fuse_check_hdr_file(currentFile);
            end
            
            % set the file name to the groups data structure
            groupsData(nGroup).feature(nFeature).files(ii).name = deblank(selectedFiles(ii).name(1, :));
            
        end
        % end loop over the selected files
        
        clear selectedFiles;
        
    end
    % end loop over features
    
end
% end loop over groups


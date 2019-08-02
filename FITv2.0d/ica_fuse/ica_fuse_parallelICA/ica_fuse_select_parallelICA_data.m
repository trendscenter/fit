function [dataInfo, numSubjects] = ica_fuse_select_parallelICA_data(inputFile)
% Function for selecting data using the information of file
% pattern


ica_fuse_defaults;
global MRI_DATA_FILTER;
global GENE_DATA_FILTER;
global PARALLEL_ICA_MODALITIES;

if ~exist('inputFile', 'var')
    inputFile = [];
end

numFeatures = 2;

if isempty(inputFile)
    
    numParameters = 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Number of groups';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '1';
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'numGroups';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Number of features';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(numFeatures);
    inputText(numParameters).answerType = 'numeric';
    inputText(numParameters).tag = 'numFeatures';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    
    % Input dialog box
    answer_groups = ica_fuse_inputDialog('inputtext', inputText, 'Title', 'Enter no. of groups');
    
    if isempty(answer_groups)
        error('Figure window to get the number of groups was quit');
    end
    
    % get the number of groups
    numGroups = answer_groups{1};
    numFeatures = answer_groups{2};
    
    if (numFeatures > 3 || numFeatures < 2)
        error('Number of features you can enter is either 2 or 3');
    end
    
    
    clear inputText;
    
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
        inputText(nFeature).answerString = PARALLEL_ICA_MODALITIES;
        inputText(nFeature).answerType = 'string';
        inputText(nFeature).tag = ['modality', num2str(nFeature)];
        inputText(nFeature).enable = 'on';
        inputText(nFeature).value = 1;
        
    end
    % end loop over features
    
    
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
    
    % Modalities
    modalities = answersInput(1:numFeatures);
    groupNamings = answersInput(numFeatures + 1:numFeatures + numGroups);
    featureNamings = answersInput(numFeatures + numGroups + 1:end);
    
else
    
    % read number of groups and features
    keywds = {'groupNames', 'featureNames', 'modality'};
    % read output directory for the analysis
    inputData = ica_fuse_read_variables(inputFile, keywds, {'cell', 'cell', 'cell'});
    groupNamings = getfield(inputData, 'groupNames'); % group names
    featureNamings = getfield(inputData, 'featureNames'); % feature names
    modalities = getfield(inputData, 'modality'); % modalities        
    
    clear inputData;
    
    numFeatures = length(featureNamings);
    
    numGroups = length(groupNamings);
    
    if ((length(featureNamings) > 3) || (length(featureNamings) < 2))
        error('Number of features you can enter is either 2 or 3');
    end
    
    if ((length(modalities) > 3) || (length(modalities) < 2))
        error('Number of modalities you can enter is either 2 or 3');
    end
    
end


% loop over groups
for nGroup = 1:numGroups
    dataInfo(nGroup).name = groupNamings{nGroup};
    % loop over features
    for nfeature = 1:numFeatures
        dataInfo(nGroup).feature(nfeature).modality = modalities{nfeature};
        dataInfo(nGroup).feature(nfeature).name = featureNamings{nfeature};
    end
    % end loop over features
end
% end loop over groups


%% File selection
if (isempty(inputFile))
    dataSelectionOption = ica_fuse_questionDialog('title', 'Data Selection', 'textbody', 'Do you want to use file selection window for selecting files?');
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


allModalities = repmat(modalities, 1, numGroups);

%% Use userinput of file names rather than selecting based on file patterns
if (dataSelectionOption == 2)
    
    if isempty(inputFile)
        %% GUI
        
        listBoxStrings = repmat(struct('name', ''), 1, numGroups*numFeatures);
        
        countF = 0;
        for nGroup = 1:numGroups
            for nFeature = 1:numFeatures
                countF = countF + 1;
                listBoxStrings(countF).name = [dataInfo(nGroup).name, ' ', dataInfo(nGroup).feature(nFeature).name];
            end
        end
        
        fileTypes = repmat({'any'}, 1, numGroups*numFeatures);
        
        isMRI = strcmpi(allModalities, 'fmri') | strcmpi(allModalities, 'smri');
        
        fileTypes(isMRI) = {'image'};
        
        tmp_files = ica_fuse_ui_select_data('title', 'Select data for fusion', 'num_subjects', ...
            numGroups, 'num_sessions', numFeatures, 'filter_string', '*.img;*.asc;*.nii', 'type_file_selection', 'multiple', ...
            'subject_string', listBoxStrings, 'files_specification', 'equal_sessions', 'fileType', fileTypes);
        
        if (isempty(tmp_files))
            error('Files are not selected for fusion');
        end
        
        countF = 0;
        for nGroup = 1:numGroups
            for nFeature = 1:numFeatures
                countF  = countF + 1;
                tmp = tmp_files(countF).name;
                for ii = 1:size(tmp, 1)
                    dataInfo(nGroup).feature(nFeature).files(ii).name = deblank(tmp(ii, :));
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
                    error(['No files passed for group ', dataInfo(nGroup).name, ' feature ', dataInfo(nGroup).feature(nFeature).name]);
                end                                
                
                if (size(tmp, 1) == 1)
                    [tmp_path, tmpF, extn] = fileparts(tmp);
                    tmp = ica_fuse_listFiles_inDir(tmp_path, [tmpF, extn]);
                    if isempty(tmp)
                        error(['Files not found for (', fullfile(tmp, [tmpF, extn]), ')']);
                    end
                end
                
                tmp = ica_fuse_rename_4d_file(tmp);
                
                if (nFeature == 1)
                    oldFileCount = size(tmp, 1);
                else
                    if (oldFileCount ~= size(tmp, 1))
                        error(['No. of files is not the same between the features for group ', dataInfo(nGroup).name]);
                    end
                end
                
                for ii = 1:size(tmp, 1)
                    dataInfo(nGroup).feature(nFeature).files(ii).name = deblank(tmp(ii, :));
                end
                
            end
        end
    end
    
    
else
    
    
    % check the question regarding the file pattern
    if numGroups > 1
        if isempty(inputFile)
            answerFilePattern = ica_fuse_questionDialog('title', 'File pattern', 'textbody', ...
                'Is the file pattern same between groups?');
        else
            keywd = 'answerFilePattern';
            inputData = ica_fuse_read_variables(inputFile, keywd, {'integer'});
            answerFilePattern = getfield(inputData, keywd);
            clear inputData;
        end
    else
        answerFilePattern = 1;
    end
    
    % add file pattern to the structure groups data
    [dataInfo] = getFilePattern(numGroups, numFeatures, answerFilePattern, dataInfo, inputFile);
    
    
    answerDir = 1;
    num_groups = 1;
    if (numGroups > 1)
        if isempty(inputFile)
            % get the answer for group dir
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
    
    if (answerDir == 0)
        num_groups = numGroups;
    end
    
    startPath = pwd;
    
    getDirGroups = cell(num_groups, 1);
    % select the directories for groups
    for nGroup = 1:num_groups
        if isempty(inputFile)
            getDirGroups{nGroup} = ica_fuse_selectEntry('typeEntity', 'directory', 'startPath', startPath, 'title', ...
                ['Select directory for ', dataInfo(nGroup).name]);
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
    end
    
    if length(getDirGroups) == 1
        getDirGroups = repmat(getDirGroups, numGroups, 1);
    end
    
    drawnow;
    
    % get the files per group
    
    [dataInfo] = getFiles(dataInfo, getDirGroups);
    
end

numSubjects = zeros(1, length(dataInfo));

% % Modality 1 and Modality 2 files
for nn = 1:numGroups
    tempFiles = str2mat(dataInfo(nn).feature(1).files.name);
    numSubjects(nn) = size(tempFiles, 1);
    modality1_files(nn).name = tempFiles;
    modality2_files(nn).name = str2mat(dataInfo(nn).feature(2).files.name);
end

timePoints = sum(numSubjects);

size_data = size(str2mat(modality2_files.name), 1);

if size_data ~= timePoints
    error('Error:Data', ['Number of subjects of ', dataInfo(1).feature(1).modality,  ' data (%d) doesn''t match with the ', dataInfo(1).feature(2).modality, ...
        ' (%d) data'], size_data, timePoints);
end


% Select locus names
temp = ica_fuse_loadData(deblank(modality2_files(1).name(1, :)), 1);
numTimePoints = size(temp, 1);


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


function [dataInfo] = getFilePattern(numGroups, numFeatures, checkSameFileP, dataInfo, inputFile)
% add file pattern to the structure dataInfo

ica_fuse_defaults;
global GENE_DATA_FILTER;
global EEG_DATA_FILTER;

if checkSameFileP == 1
    num_groups = 1;
else
    num_groups = numGroups;
end

handle_visibility = 'on';
if ~isempty(inputFile)
    handle_visibility = 'off';
end

fileP = '*.img';

% loop over groups
for nGroup = 1:num_groups
    
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
    
    % loop over features
    for nfeature = 1:numFeatures
        if strcmpi(handle_visibility, 'off')
            fileP = answer_file_pattern{nfeature};
        else
            if strcmpi(dataInfo(nGroup).feature(nfeature).modality, 'fmri') || strcmpi(dataInfo(nGroup).feature(nfeature).modality, 'smri')
                fileP = '*.img';
            elseif strcmpi(dataInfo(nGroup).feature(nfeature).modality, 'gene')
                fileP = GENE_DATA_FILTER;
            elseif strcmpi(dataInfo(nGroup).feature(nfeature).modality, 'eeg')
                fileP = EEG_DATA_FILTER;
            end
        end
        % define all the input parameters in a structure
        inputText(nfeature).promptString = ['Select file pattern for ', dataInfo(nGroup).name, ' ', ...
            dataInfo(nGroup).feature(nfeature).name, ' feature'];
        inputText(nfeature).uiType = 'edit';
        inputText(nfeature).answerString = fileP;
        inputText(nfeature).answerType = 'string';
        inputText(nfeature).tag = ['feature', num2str(nfeature)];
        inputText(nfeature).enable = 'on';
        inputText(nfeature).value = 1;
    end
    % end loop over features
    
    % Input dialog box
    answer_file_pattern = ica_fuse_inputDialog('inputtext', inputText, 'Title', ['Select file pattern for group ', ...
        num2str(nGroup)], 'handle_visibility', handle_visibility);
    
    
    if isempty(answer_file_pattern)
        error(['file pattern for group ', num2str(nGroup), ' input dialog box was quit']);
    end
    
    answerFileP(nGroup).name = answer_file_pattern;
    
    clear inputText;
    
end
% end loop over groups

% replicate the file pattern over groups
if num_groups == 1
    answerFileP = repmat(answerFileP, numGroups, 1);
end

% replicate the file pattern over groups
for nGroup = 1:numGroups
    
    % loop over features
    for nfeature = 1:numFeatures
        dataInfo(nGroup).feature(nfeature).filePattern = answerFileP(nGroup).name{nfeature};
    end
    % end loop over features
    
end
% end loop over groups

function [dataInfo] = getFiles(dataInfo, getDirGroups)
% search files in the directory

for nGroup = 1:length(dataInfo)
    % get the number of features
    numFeatures = length(dataInfo(nGroup).feature);
    currentDir = getDirGroups{nGroup};
    % list files for the current feature for the group
    for nfeature = 1:numFeatures
        fileP = dataInfo(nGroup).feature(nfeature).filePattern;
        [selectedFiles] = ica_fuse_listFiles_inDir(currentDir, fileP);
        if ~isempty(selectedFiles)
            for ii = 1:size(selectedFiles, 1)
                dataInfo(nGroup).feature(nfeature).files(ii).name = fullfile(currentDir, deblank(selectedFiles(ii, :)));
            end
            
        else
            try
                disp(['Reading data from sub folders in directory ', currentDir]);
                selectedFiles = ica_fuse_get_sub_data(currentDir, fileP, 'data_in_subject_folder', 0);
            catch
                disp(['Reading data from sub sub-folders in directory ', currentDir]);
                selectedFiles = ica_fuse_get_sub_data(currentDir, fileP, 'data_in_subject_subfolder', 0);
            end
            for ii = 1:length(selectedFiles)
                dataInfo(nGroup).feature(nfeature).files(ii).name = deblank(selectedFiles(ii).name(1, :));
            end
        end
        clear selectedFiles;
        
    end
end
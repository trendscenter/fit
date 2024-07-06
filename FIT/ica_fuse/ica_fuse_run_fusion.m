function ica_fuse_run_fusion(inputFile)
% Run joint ICA fusion analysis using the input batch file. This function
% gives you the flexibility that data files with the matching file pattern can be nested
% anywhere in a directory.
%
% NOTE: Run this function only on UNIX as this uses find command to get
% the files.
%
% Inputs:
%
% inputFile - Input batch file.
%

if ~isunix
    error('ica_fuse_run_fusion function can only be run on UNIX operating system');
end

% Load defaults
ica_fuse_defaults;
global FUSION_INFO_MAT_FILE;

[pathstr, fName, extn] = fileparts(inputFile);

if isempty(pathstr)
    pathstr = fileparts(which(inputFile));
end

addpath(pathstr);

inputFile = fullfile(pathstr, [fName, extn]);

eval(fName);

% Fusion file
fusionFile = fullfile(outputDir, [prefix, FUSION_INFO_MAT_FILE, '.mat']);

% Do error check
if (length(modality) ~= length(featureNames))
    error('Error:NumModalities', 'Number of modalities (%d) does not match number of features (%d)', length(modality), ...
        length(featureNames));
end


% Initialise files
files = cell(length(groupNames), length(featureNames));

% Initialise data Info
dataInfo = repmat(struct('name', ''), 1, length(groupNames));

% Loop over groups
for nG = 1:length(groupNames)

    numFiles = zeros(1, length(featureNames));

    % Check the file pattern
    if answerFilePattern
        currentFilePattern = group1_file_pattern;
    else
        currentFilePattern = eval(['group', num2str(nG), '_file_pattern']);
    end

    % Group directory
    if answerDir
        currentDir = group1_dir;
    else
        currentDir = eval(['group', num2str(nG), '_dir']);
    end

    dataInfo(nG).name = groupNames{nG};

    % Loop over features
    for nF = 1:length(featureNames)
        disp(['Listing files of group ', num2str(nG), ' feature ', num2str(nF)]);
        dataInfo(nG).feature(nF).name = featureNames{nF};
        dataInfo(nG).feature(nF).filePattern = currentFilePattern{nF};
        dataInfo(nG).feature(nF).modality = modality{nF};
        txtFile = fullfile(outputDir, ['group', num2str(nG), '_feature', num2str(nF), '_files.txt']);
        sysCommand = ['find ', currentDir, ' -name ''', currentFilePattern{nF}, ''' > ', txtFile];
        [status] = system(sysCommand);
        if (status == 1)
            error(['Files not found with file pattern ', currentFilePattern{nF}, ' for group ', currentDir]);
        end
        current_files = textread(txtFile, '%s', 'delimiter', '\n');
        files{nG, nF} = current_files;
        clear current_files
        disp('Done listing');
        fprintf('\n');
    end

    check = find(numFiles == numFiles(1));
    if length(check) ~= length(numFiles)
        error(['Number of files for group ', currentDir, ' between features are unequal']);
    end
    % End loop over features

end
% End loop over groups


% Loop over groups
for nG = 1:length(groupNames)
    % Loop over features
    for nF = 1:length(featureNames)
        current_files = files{nG, nF};
        % Loop over subjects
        for nSub = 1:length(current_files)
            dataInfo(nG).feature(nF).files(nSub).name = current_files{nSub};
        end
        % End loop over subjects
        clear current_files;
    end
    % End loop over features
end
% End loop over groups

clear files;

% Add dataInfo field
fusionInfo.setup_analysis.dataInfo = dataInfo;

clear dataInfo;

% Save fusion file
ica_fuse_save(fusionFile, 'fusionInfo');

% Setup analysis
[fusionFile] = ica_fuse_setup_analysis(inputFile, fusionFile);

% Run analysis
ica_fuse_run_analysis(fusionFile);
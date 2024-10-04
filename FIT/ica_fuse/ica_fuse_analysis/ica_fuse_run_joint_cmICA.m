function ica_fuse_run_joint_cmICA(fileName)
%% Run Joint cmICA
%

ica_fuse_defaults;
global GREY_MATTER;
global WHITE_MATTER;

if (isempty(GREY_MATTER) && (WHITE_MATTER))
    error('Grey and white matter masks are not specified in ica_fuse_defaults');
end

if (~exist('fileName', 'var'))
    fileName = ica_fuse_selectEntry('title', 'Select joint cmICA file', 'filter', '*cmica_info.mat', 'typeEntity', 'file');
end

drawnow;

if (isempty(fileName))
    error('cmICAInfo file is not selected for analysis');
end

if (ischar(fileName))
    load(fileName);
    if (~exist('cmICAInfo', 'var'))
        error(['Selected file ', fileName, ' is not a valid joint cmICA file']);
    end
    [outputDir, fName, extn] = fileparts(fileName);
    
else
    cmICAInfo = fileName;
    outputDir = cmICAInfo.outputDir;
end

if (isempty(outputDir))
    outputDir = pwd;
end

if (length(cmICAInfo.dataInfo) < 2)
    error('Need to select at least two features in order to do joint cmICA');
end

%% Get params
outputPrefix = cmICAInfo.output_prefix;
numPC1 = cmICAInfo.numPC1;
numComp = cmICAInfo.numComp;
algorithm = cmICAInfo.algorithm;
dataInfo = cmICAInfo.dataInfo;
cmICAInfo.outputDir = outputDir;

%% Record analysis information in a log file
diaryFile = fullfile(outputDir, [outputPrefix, '_cmica_results.log']);
diary(diaryFile);

tic;

fprintf('\n');

disp('Starting Joint cmICA Analysis ...');

fprintf('\n');

%% Data Reduction
ica_fuse_cmICA_dataReduction(cmICAInfo);

%% Calculate ICA
ica_fuse_cmICA_calculateICA(cmICAInfo);

%% Calculate ICA
ica_fuse_cmICA_backReconstruct(cmICAInfo);

%% Scale components
ica_fuse_cmICA_scaleComponents(cmICAInfo);

fusionFile = fullfile(outputDir, [outputPrefix, '_joint_cmica_info.mat']);

fprintf('\n');
disp(['Fusion information is saved in file: ', fusionFile]);

t_end = toc;

disp(['Time taken to run the analysis is ', num2str(t_end), ' seconds']);
fprintf('\n');
disp(['All the analysis information is stored in log file: ', diaryFile]);
fprintf('\n');

diary('off');


diary('off');

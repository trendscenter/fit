function fusionInfo = ica_fuse_backReconstruct(fusionInfo, data)
%% Back reconstruction step. All the information will be stored in MAT
% files. The MAT file name is dependent on the global variable
% BACK_RECONSTRUCT_FILE (see ica_fuse_defaults.m) and combination number.
%
% Inputs:
% 1. fusionInfo - Fusion Info data structure
% 2. data - Original data
%
% Outputs:
% fusionInfo - Fusion Info data structure with additional field backReconstructFiles.
%

ica_fuse_defaults;
global BACK_RECONSTRUCT_FILE;

% Output directory
outputDir = fusionInfo.run_analysis.outputDir;
combinationName = fusionInfo.run_analysis.currentCombName;
comb_number = fusionInfo.run_analysis.currentComb;

ica_algo = ica_fuse_icaAlgorithm;
ica_algo = cellstr(ica_algo);
algorithmName = ica_algo{fusionInfo.run_analysis.algorithm};

%% Load PCA file
pcaFile = fullfile(outputDir, fusionInfo.run_analysis.pcaFiles(comb_number).name);
load(pcaFile, 'dewhiteM');

%% Load ICA file
icaFile = fullfile(outputDir, fusionInfo.run_analysis.icaFiles(comb_number).name);
load(icaFile, 'A', 'icasig');

disp('-----------------------------------------------------------------------------------------------');
disp(['Doing Back Reconstruction for ',  combinationName]);
disp('-----------------------------------------------------------------------------------------------');

fprintf('\n');

%% Multiply with dewhitening matrix
if ~iscell(dewhiteM)
    A = dewhiteM*A;
else
    
    if (~iscell(A))
        A = repmat({A}, length(dewhiteM), 1);
    end
    
    newA = cell(1, length(dewhiteM));
    for n = 1:length(dewhiteM)
        newA{n} = dewhiteM{n}*A{n};
    end
    
    A = newA;
    clear newA;
    
end

if (strcmpi(algorithmName, 'iva-g') || strcmpi(algorithmName, 'iva-ggd'))
    
    icasigN = icasig;
    AN = A;
    clear A icasig;
    icasig = cat(1, AN{:});
    icasig = icasig';
    A = cell(size(icasigN));
    for nM = 1:length(icasigN)
        A{nM} = icasigN{nM}';
    end
    
end

%%%%%%%% Saving ICA Information %%%%%%%%%%%%%%%
brFile = [fusionInfo.run_analysis.prefix, BACK_RECONSTRUCT_FILE, 'comb_', num2str(comb_number), '.mat'];
fusionInfo.run_analysis.backReconstructFiles(comb_number).name = brFile;
fusionInfo.run_analysis.backReconstructFiles(comb_number).combinationName = combinationName;
brFile = fullfile(outputDir, brFile);


%% Number of groups
numGroups = fusionInfo.run_analysis.numGroups;

%% Number of subjects
numSubjects = fusionInfo.run_analysis.numSubjects;

%% Group back reconstruction
if ~iscell(A)
    groups_icasig = ica_fuse_backReconstruct_groups(A, data, numGroups, numSubjects);
else
    groups_icasig = cell(1, length(A));
    endN = 0;
    for n = 1:length(A)
        startN = endN + 1;
        endN = endN + fusionInfo.run_analysis.newDims(n);
        groups_icasig{n} = ica_fuse_backReconstruct_groups(A{n}, data(:, startN:endN), numGroups, numSubjects);
    end
    groups_icasig = cat(1, groups_icasig{:});
end

ica_fuse_save(brFile, 'A', 'icasig', 'groups_icasig', 'combinationName');

disp(['Back reconstruction information for ', combinationName, ' is saved in ', brFile]);

fprintf('\n');

%% Save fusion file
fusionFile = fusionInfo.run_analysis.fusionFile;
ica_fuse_save(fusionFile, 'fusionInfo');


disp('-----------------------------------------------------------------------------------------------');
disp(['Done Back Reconstruction for ',  combinationName]);
disp('-----------------------------------------------------------------------------------------------');

fprintf('\n');
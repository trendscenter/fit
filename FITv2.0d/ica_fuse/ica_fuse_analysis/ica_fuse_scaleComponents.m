function [icasig, A, groups_icasig, fusionInfo] = ica_fuse_scaleComponents(fusionInfo, meanData, meanDataGroups, featureDataLength)
%% Scale components based on regression and flip the sign of the component
% if the beta weight is negative. The information will be stored in MAT
% files depending on the global variable SCALE_COMP_FILE (See
% ica_fuse_defaults.m).
%
% Inputs:
% 1. fusionInfo - Fusion information file
% 2. meanData - Mean data for each features
% 3. meanDataGroups - Mean data for each feature and group
%
% Outputs:
% 1. icasig - Scaled icasig
% 2. FusionInfo variable
%

% Output directory
outputDir = fusionInfo.run_analysis.outputDir;
combinationName = fusionInfo.run_analysis.currentCombName;
comb_number = fusionInfo.run_analysis.currentComb;

disp('-----------------------------------------------------------------------------------------------');
disp(['Doing scaling components for ',  combinationName]);
disp('-----------------------------------------------------------------------------------------------');

fprintf('\n');


ica_fuse_defaults;
global SCALE_COMP_FILE;
global DETREND_NUMBER;
global FLIP_SIGN_COMPONENTS;

%% Load back reconstruction file
brFile = fullfile(outputDir, fusionInfo.run_analysis.backReconstructFiles(comb_number).name);
load(brFile, 'A', 'icasig', 'groups_icasig');

%featureNormPara = fusionInfo.run_analysis.featureNormPara;

%stdParameters = fusionInfo.run_analysis.stdParameters;

%numComp = fusionInfo.run_analysis.numComp;

numSubjects = fusionInfo.run_analysis.numSubjects;

numGroups = fusionInfo.run_analysis.numGroups;

ConvertToZ = fusionInfo.run_analysis.z_scores;

groupNames = cellstr(str2mat(fusionInfo.run_analysis.dataInfo.name));

% feature names (Force the character not have & character)
%featureNames = strread(combinationName, '%s', 'delimiter', '&');
featureNames = cellstr(str2mat(fusionInfo.run_analysis.dataInfo(1).feature.name));

% all Modalities
modalities = cellstr(str2mat(fusionInfo.run_analysis.dataInfo(1).feature.modality));

if isempty(FLIP_SIGN_COMPONENTS)
    flip_components = 1;
else
    flip_components = FLIP_SIGN_COMPONENTS;
end

if flip_components
    flipFlag = 'flip';
else
    flipFlag = 'no-flip';
end

if ~isempty(groups_icasig)
    % Loop over groups
    for nGroup = 1:numGroups
        % Make sure that joint component and back-reconstructed have
        % the same sign
        for ii = 1:size(icasig, 1)
            data1 = detrend(icasig(ii, :), 0);
            data2 = squeeze(detrend(groups_icasig(:, ii, nGroup), 0));
            correlationValue = data1*data2/(norm(data1)*norm(data2));
            if sign(correlationValue) == -1
                str = ['Difference in Sign Between Reference Image and Group ', groupNames{nGroup}, ...
                    ' Component ', num2str(ii)];
                disp([str, ' -> Changing Sign of Group ', groupNames{nGroup}, ' Component ', num2str(ii)]);
                groups_icasig(:, ii, nGroup) = sign(correlationValue)*...
                    groups_icasig(:, ii, nGroup);
            end
        end
    end
    % End loop over groups
end

disp('Scaling joint components to percent signal change ...');
% [icasig, betaWeights, betaWeightStr, featureCompStr] = ica_fuse_scaleCompData(data, icasig, featureDataLength, ...
%     fusionInfo, 'mean_data', [], flipFlag);

[icasig, betaWeights, betaWeightStr, featureCompStr] = ica_fuse_scaleCompData(meanData, icasig, featureDataLength, featureNames, flipFlag, ConvertToZ);

% Loop over components
if strcmpi(flipFlag, 'flip')
    
    if (~iscell(A))
        
        A = A .* repmat(sign(betaWeights(1).value(:)'), size(A, 1), 1);
        
    else
        
        for i = 1:length(A)
            A{i} = A{i} .* repmat(sign(betaWeights(i).value(:)'), size(A{i}, 1), 1);
        end
        
    end
    
end
% End loop over components

disp('Done scaling joint components');

% % Number of components & voxels
% numComp = size(icasig, 1);
% voxels = size(icasig, 2);

if ~isempty(groups_icasig)
    fprintf('\n');
    
    disp('Scaling group back-reconstructed components to percent signal change ...');
    
    % Loop over groups
    for nGroups = 1:numGroups
        fprintf('\n');
        disp(['Doing scaling on group ', groupNames{nGroups}, ' ...']);
        groupInd = ica_fuse_get_groupInd(nGroups, numSubjects);
        currentData = groups_icasig(:, :, nGroups);
        currentData = currentData';
        currentData = ica_fuse_scaleCompData(meanDataGroups{nGroups}, currentData, featureDataLength, featureNames, flipFlag, ConvertToZ);
        
        groups_icasig(:, :, nGroups) = currentData';
        disp(['Done scaling on group ', groupNames{nGroups}]);
        fprintf('\n');
        
    end
    
    fprintf('\n');
    
    fprintf('\n');
    
    % End loop over groups
    disp('Done scaling group back-reconstructed components');
    fprintf('\n');
end

%%%%%%%% Saving ICA Information %%%%%%%%%%%%%%%
scaleFile = [fusionInfo.run_analysis.prefix, SCALE_COMP_FILE, 'comb_', num2str(comb_number), '.mat'];
fusionInfo.run_analysis.scaleCompFiles(comb_number).name = scaleFile;
fusionInfo.run_analysis.scaleCompFiles(comb_number).combinationName = combinationName;
scaleFile = fullfile(outputDir, scaleFile);

% save the data
ica_fuse_save(scaleFile, 'icasig', 'A', 'groups_icasig', 'combinationName', 'betaWeights');

fprintf('\n');
disp(['Scaling components information for ', combinationName, ' is saved in ', scaleFile]);

fprintf('\n');

% Print all the beta weights information to a file
% Store the data information
numPara = 1;
varStruct(numPara).tag = 'Component';
varStruct(numPara).value = featureCompStr;

numPara = numPara + 1;
varStruct(numPara).tag = 'Beta Weight (Component fit to mean data of that feature)';
varStruct(numPara).value = betaWeightStr;

clear betaWeightStr; clear featureCompStr;

% Form strings to print to a file
titlePrint = ['Beta weight information of combination ', combinationName];

betaFile = [fusionInfo.run_analysis.prefix, '_beta_weight_comb_', num2str(comb_number), '.txt'];
betaFile = fullfile(outputDir, betaFile);

ica_fuse_printToFile(betaFile, varStruct, titlePrint, 'column_wise', 'append');

disp(['Beta weight information for ', combinationName, ' is saved in ', betaFile])

fprintf('\n');

%% Save fusion file
fusionFile = fusionInfo.run_analysis.fusionFile;
ica_fuse_save(fusionFile, 'fusionInfo');

disp('-----------------------------------------------------------------------------------------------');
disp(['Done scaling components for ',  combinationName]);
disp('-----------------------------------------------------------------------------------------------');

fprintf('\n');
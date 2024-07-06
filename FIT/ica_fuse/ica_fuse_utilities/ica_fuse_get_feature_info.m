function featureInfo = ica_fuse_get_feature_info(dataInfo, selGroups, selFeatures)
%% Get the files of selected groups and features
%
% Inputs:
% 1. dataInfo - Data information structure
% 2. selGroups - Selected Groups
% 3. selFeatures - Selected features
% 
% Outputs:
% featureInfo - featureInfo data structure
%

%% Number of groups
numGroups = length(dataInfo);

%% Number of features
numFeatures = length(dataInfo(1).feature);

if (~exist('selGroups', 'var'))
    selGroups = (1:numGroups);
end

if (~exist('selFeatures', 'var'))
    selFeatures = (1:numFeatures);
end

%% Initialise featureInfo data structure
featureInfo = repmat(struct('files', [], 'feature_name', [], 'modality', []), 1, length(selFeatures));

%% Loop over no. of selected features
for nFeature = 1:length(selFeatures)
    tempF = repmat(struct('name', []), 1, length(selGroups));
    %% Loop over groups
    for nGroup = 1:length(selGroups)
        tempF(nGroup).name = str2mat(dataInfo(selGroups(nGroup)).feature(selFeatures(nFeature)).files.name);
    end
    %% End loop over groups
    totalFiles = str2mat(tempF.name);
    clear tempF;
    featureInfo(nFeature).files = totalFiles;    
    featureInfo(nFeature).feature_name = dataInfo(1).feature(selFeatures(nFeature)).name;
    featureInfo(nFeature).modality = dataInfo(1).feature(selFeatures(nFeature)).modality;
    clear totalFiles;
end
%% End loop over no. of selected features




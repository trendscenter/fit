function varargout = ica_fuse_calculate_histogram(histParameters, histogramCriteria, threshold, combNumber, selFeatures, selComp)
% Calculate histogram of features or components.
%
% Inputs:
% 1. histParameters - Structure containing necessary fields to generate
% histograms.
% 2. histogramCriteria - Feature or component.
% 3. threshold - Z-threshold
%
% Feature Histograms: Feature histogram is calculated by Z-thresholding the best
% component which is obtained by sorting (ttest2 on mixing coefficients or spatial divergence).
% The resulting voxels are sorted in descending order and this is used as a mask to generate
% histograms of individual subjects.
%
% Component Histograms: Component voxels are sorted in descending order and
% histograms are calculated.
%
% Output:
% varargout - varargout contains histData where histData is a structure
% containing fields like histData(1).group(1).data. Length of histData depends on number of
% feature combinations. Length of histData(1).group depends on number of
% groups.
%
% Note: For feature histograms mean (meanHist variable) of each group is also returned.
% meanHist is a structure containing data as field. Length of meanHist
% depends on the number of combinations.
%

% if isempty(which('ttest2.m'))
%     error('Need statistics toolbox to plot histograms for subjects');
% end

%% Load defaults
ica_fuse_defaults;

global Z_THRESHOLD_HISTOGRAM;

if ~exist('histParameters', 'var')
    error('Parameters structure required for generating histograms is not provided');
end

%% Histogram parameters
if ~isfield(histParameters, 'outputDir')
    outputDir = pwd;
else
    outputDir = histParameters.outputDir;
end

% Group names
groupNames = histParameters.groupNames;

% dataInfo
dataInfo = histParameters.dataInfo;

% Back reconstruct files
back_reconstruct_files = histParameters.backReconstructFiles;

% Voxel indices
maskIndices = histParameters.mask_ind;

maskIndices = ica_fuse_form_maskInd(maskIndices, dataInfo);

% Number of subjects
numSubjects = histParameters.numSubjects;

if length(numSubjects) ~= length(dataInfo)
    numSubjects = repmat(numSubjects, 1, length(dataInfo));
end

% Feature data length
featureDataLength = histParameters.featureDataLength;

% Normalization type
normalizeVal = histParameters.normalize;
max_voxels = histParameters.voxels;
all_comb = histParameters.all_comb;

% Feature normalization parameters
featureNormPara = histParameters.featureNormPara;

selGroups = histParameters.selGroupsVal;

if ~exist('threshold', 'var')
    threshold = Z_THRESHOLD_HISTOGRAM;
end

if ~exist('histogramCriteria', 'var')
    histogramCriteria = 'feature';
end

drawnow;

%% Check if the components need to be sorted
sort_comp = 1;
if exist('selComp', 'var')
    sort_comp = 0;
end

%% Check if data already exists
if isfield(histParameters, 'dataN')
    dat = histParameters.dataN;
end

histParameters.dataN = [];

selSubjects = numSubjects(selGroups);

%% MAT file combination number
if exist('combNumber', 'var')
    combInd = combNumber;
else
    combInd = 1:length(back_reconstruct_files);
end

%% Initialise histogram data
histData = repmat(struct('group', [], 'bestComp', [], 'combinationName', []), 1, length(combInd));
meanHist = repmat(struct('data', []), 1, length(combInd));

countBrFile = 0;

%% Loop over number of combinations
for nn = combInd

    countBrFile = countBrFile + 1;

    %% Get combination numbers
    currentCombNum = all_comb{nn};

    if isempty(currentCombNum)
        continue;
    end

    %% Current back reconstruction file
    currentFile = fullfile(outputDir, deblank(back_reconstruct_files(nn).name));

    dataLength = featureDataLength(nn).Length;

    load(currentFile, 'icasig', 'combinationName');

    % Relative indices of features
    if (exist('selFeatures', 'var'))
        featureNum = selFeatures;
    else
        featureNum = 1:length(dataLength);
    end

    %% Truncate component voxels
    startN = 1;
    endN = 0;
    allFeaInd = zeros(1, sum(dataLength(featureNum)));
    for nFea = 1:length(featureNum)
        [featureStartInd, featureEndInd] = ica_fuse_get_featureInd(featureNum(nFea), dataLength);
        tempInd = (featureStartInd:featureEndInd);
        endN =  endN + length(tempInd);
        allFeaInd(startN:endN) = tempInd;
        clear tempInd;
        startN = endN + 1;
    end

    clear startN endN featureStartInd featureEndInd;

    icasig = icasig(:, allFeaInd);


    % Truncate data length
    dataLength = dataLength(featureNum);


    if ica_fuse_findstr(histogramCriteria, 'feature')

        %% Load data if it doesn't exist
        if ~exist('dat', 'var')
            % Stack information
            stackInfo = ica_fuse_prepare_data('dataInfo', dataInfo, 'mask_ind', maskIndices, 'normalize_scheme', ...
                normalizeVal, 'voxels', max_voxels, 'sel_features', currentCombNum(featureNum), ...
                'sel_groups', selGroups);

            dat = stackInfo.data;
            clear stackInfo;
        end
    end

    selectedFeatures = str2mat(strread(combinationName, '%s', 'delimiter', '&'));
    histParameters.selectedFeatureVal = featureNum;
    histParameters.selectedFeature = selectedFeatures(featureNum, :);

    % Sort components
    if sort_comp
        [sortResults] = ica_fuse_sort_components(histParameters, nn, 0);
        selComp = sortResults.sorted_comp(1);
        clear sortResults;
    end


    if ica_fuse_findstr(histogramCriteria, 'feature')
        % Use only selected component
        icasig = icasig(selComp, :);
        [histData(countBrFile).group, meanHist(countBrFile).data, histData(countBrFile).hIndex] = ica_fuse_compute_hist(dat, dataLength, selSubjects, 'feature', 'icasig', icasig, ...
            'sel_comp', selComp, 'feature_names', histParameters.selectedFeature, 'group_names', groupNames(selGroups, :), 'threshold', threshold);
    else
        % Load back reconstructed components
        load(currentFile, 'groups_icasig');
        % Select only the component of interest
        groups_icasig = squeeze(groups_icasig(allFeaInd, selComp, selGroups))';
        histData(countBrFile).group = ica_fuse_compute_hist(groups_icasig, dataLength, selSubjects, 'component', 'feature_names', histParameters.selectedFeature, 'group_names', ...
            groupNames(selGroups, :));
    end

    clear dat;

    histData(countBrFile).bestComp = selComp;
    histData(countBrFile).combinationName = combinationName;
    clear groups_icasig A;

    fprintf('\n');

end
%% End loop over number of combinations


%% Return Output
varargout{1} = histData;

%% Return mean of histograms when feature is selected
if ica_fuse_findstr(histogramCriteria, 'feature')
    varargout{2} = meanHist;
end
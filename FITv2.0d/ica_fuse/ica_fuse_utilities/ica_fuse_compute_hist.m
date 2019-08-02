function varargout = ica_fuse_compute_hist(data, dataLength, subjectsGroups, histogramType, varargin)
%% Function to compute histograms. There are two types of histograms like
% feature and component. 
% 
% Inputs:
% 1. Data - Data is dependent on the type of histogram.
%       a. Feature - For feature histograms, original data is required. The
%       size of the data is total subjects by total voxels. Total subjects is
%       the sum of subjectsGroups vector. Total voxels is equal to sum of
%       dataLength vector
%       b. Component - For component histograms, group components is
%       required. The size of group components is equal to number of
%       groups by totalVoxels.
%
% 2. dataLength - Vector containing voxels information for each feature. The length of vector dataLength 
% is equal to the number of features.
%
% 3. subjectsGroups - Number of subjects for each group in a row vector.
% The length of number of subjects vector is equal to the number of groups.
% Sum of the subjectsGroups vector is equal to the total number of
% subjects.
%
% 4. histogramType - There are two types of histograms like feature and
% component. Options are 'feature' and 'component'.
%
% 5. varargin - Arguments must be passed in pairs. For feature histograms,
% you need to provide the component data you want to use in a row vector. 
%

%% Load defaults
ica_fuse_defaults;

global Z_THRESHOLD_HISTOGRAM;
global NUM_BINS;

%% Optional parameters number
optionalVarsNum = length(varargin);

%% Do error check
if (mod(optionalVarsNum, 2) ~= 0)
    error('Error:HistErr', 'Variables after histogramType must be passed in pairs');
end

selComp = '';

%% Loop over variables
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'icasig')
        icasig = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'feature_names')
        featureNames = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'sel_comp')
        selComp = num2str(varargin{i + 1});
    elseif strcmpi(varargin{i}, 'group_names')
        groupNames = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'threshold')
        threshold = varargin{i + 1};
    end
end

%% Check data length
totalVoxels = sum(dataLength);
totalSubjects = sum(subjectsGroups);
numGroups = length(subjectsGroups);
numFeatures = length(dataLength);

%% Check optional parameters
if strcmpi(histogramType, 'feature')
    if ~exist('icasig', 'var')
        error('Error:HistErr', 'You need to provide the component data when using feature histograms');
    end               
    %% Reshape component signal to row vector
    if exist('icasig', 'var')
        if (size(icasig, 2) ~= totalVoxels)
            icasig = icasig(:)';
        end
        if (length(icasig) ~= totalVoxels)
            error('Error:HistErr', 'Total number of voxels (%d) is not equal to the component data voxels (%d)\n', totalVoxels, size(icasig, 2));
        end
    end
    
    %% Do error check on data
    if (size(data, 1) ~= totalSubjects)
        error('Error:HistErr', 'No. of rows of data (%d) must equal the total number of subjects (%d)\n', size(data, 1), totalSubjects);
    end
    
else
    
    %% Do error check on data
    if (size(data, 1) ~= numGroups)
        error('Error:HistErr', 'No. of rows of data (%d) must equal the number of groups (%d)\n', size(data, 1), numGroups);
    end
    
end


if (size(data, 2) ~= totalVoxels)
    error('Error:HistErr', 'No. of cols of data (%d) must equal the total no. of voxels (%d)\n', size(data, 2), totalVoxels);
end

%% Error check for no. of features
if ~exist('featureNames', 'var')
    % Default feature names
    featureNames = [repmat('Feature ', numFeatures, 1), num2str((1:numFeatures)')];
end

if ~exist('groupNames', 'var')
    % Default group names
    groupNames = [repmat('Group ', numGroups, 1), num2str((1:numGroups)')];
end

groupNames = cellstr(groupNames);
featureNames = cellstr(featureNames);

if (length(featureNames) ~= numFeatures)
    error('Error:HistErr', 'No. of feature names (%d) must match the number of features (%d)\n', length(featureNames), numFeatures);
end

if ~exist('threshold', 'var')
    threshold = Z_THRESHOLD_HISTOGRAM;
end


% Initialise histogram data for groups
groupHist = repmat(struct('data', []), 1, numGroups);

% Feature histograms
if strcmpi(histogramType, 'feature')
    
    fprintf('\n');
    
    disp(['Using component ', selComp, ' indices to calculate histograms of subjects ']);
    disp('and the Z-threshold will be applied ...');
    
    indices = cell(1, numFeatures);
    
    %% Loop over number of features
    for countFeature = 1:numFeatures
        [featureStartInd, featureEndInd] = ica_fuse_get_featureInd(countFeature, dataLength);
        group_task = icasig(featureStartInd:featureEndInd);
        % Convert component voxel values to z-scores
        group_task = group_task ./ std(group_task(:));
        % Sort component voxel values max to min
        [g1, xIndex] = sort(abs(group_task(:)));
        xIndex = xIndex(end:-1:1);
        ind = abs(group_task(xIndex)) >= threshold;
        xIndex = xIndex(ind);
        if isempty(xIndex)
            error('Error:ZThreshold', ['Cannot find the voxels in the component for feature %s with a z-threshold of %d.', ...
                    '\nPlease set a threshold lower than %d.'], featureNames{countFeature}, threshold, ...
                max(abs(round(group_task(:)))));
        end
        clear ind;
        if countFeature == 1
            numVoxels = length(xIndex);
        else
            numVoxels = min([numVoxels, length(xIndex)]);
        end
        indices{countFeature} = xIndex;
        clear group_task;
    end
    %% End loop over features
    
    %% Use top component voxels based on information of all features
    hIndex = cell(1, numFeatures);
    for countFeature = 1:numFeatures
        [featureStartInd, featureEndInd] = ica_fuse_get_featureInd(countFeature, dataLength);
        ind = indices{countFeature};
        ind = ind(1:numVoxels);
        indices{countFeature} = ind;
        temp = data(:, featureStartInd:featureEndInd);
        temp = temp(:, ind);
        % Compute histogram indices based on number of bins
        xind = linspace(min(temp(:)), max(temp(:)), NUM_BINS);
        hIndex{countFeature} = xind;
        clear ind xind temp;
    end
    %%% End for getting top component voxels %%%%%
    
    %% Loop over selected groups
    for nGroup = 1:numGroups
        groupInd = ica_fuse_get_groupInd(nGroup, subjectsGroups);
        countN = 0;
        % Loop over subjects
        for nSub = groupInd
            countN = countN + 1;
            hData = cell(1, numFeatures);
            % Loop over selected features
            for countFeature = 1:numFeatures
                [featureStartInd, featureEndInd] = ica_fuse_get_featureInd(countFeature, dataLength);
                temp = data(nSub, featureStartInd:featureEndInd);
                hData{countFeature} = temp(indices{countFeature});
            end
            % End loop over selected features
            
            disp(['Calculating cross task histograms for group ', groupNames{nGroup}, ' subject ', num2str(countN)]);
            
            % Generate ND histogram
            out = ica_fuse_histnd(hData, hIndex);
            groupHist(nGroup).data(:, countN) = out(:);
            size_out = size(out);
            clear hData out;
        end
        % End loop over selected subjects
        groupHist(nGroup).data = reshape(groupHist(nGroup).data, [size_out, length(groupInd)]);
    end
    %% End loop over selected groups
    
    clear data;
    
    %% Mean histogram
    meanHist = cell(numGroups, 1);
    
    % Calculate mean histogram
    for nGroup = 1:numGroups
        meanHist{nGroup} = mean(groupHist(nGroup).data, length(size(groupHist(nGroup).data)));
    end
    
    % Store histogram for each subject and group
    varargout{1} = groupHist;
    % Store mean histogram for each group
    varargout{2} = meanHist;
    % Store voxel values
    varargout{3} = hIndex;
    
else
    
    % Format String
    formatStr = sprintf(repmat('%s ', 1, numFeatures), featureNames{:});
    
    % Component histograms
    disp(['Calculating cross task histograms for features ', formatStr, '...']);
    
    % Loop over groups
    for nGroup = 1:numGroups
        hData = cell(1, numFeatures);
        % Loop over number of features
        for countFeature = 1:numFeatures
            [featureStartInd, featureEndInd] = ica_fuse_get_featureInd(countFeature, dataLength);
            hData{countFeature} = data(nGroup, featureStartInd:featureEndInd);
        end
        % End loop over features
        clear countFeature;
        % Form ND histogram where N refers to the number of features
        out = ica_fuse_histnd(hData);
        groupHist(nGroup).data = out;
        clear hData out;
    end
    % End loop over selected groups
    
    % Store histogram for each subject and group
    varargout{1} = groupHist;
    
end
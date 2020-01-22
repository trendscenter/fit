function [dataN] = ica_fuse_applyMask(featureInfo, mask_ind)
% Apply Mask

ica_fuse_defaults;
global EEG_DATA_INDICES;

dataN = repmat(struct('data', [], 'xAxis', [], 'files', [], 'feature_name', []), 1, length(featureInfo));

for nFeature = 1:length(featureInfo)
    ind = mask_ind(nFeature).ind; % Mask indices for the corresponding feature
    length_ind = length(find(ind ~= 0));
    if strcmpi(featureInfo(nFeature).modality, 'fmri') || strcmpi(featureInfo(nFeature).modality, 'smri')
        % get data for each feature
        tempV = ica_fuse_spm_vol(featureInfo(nFeature).files);
        currentData = zeros(length_ind, length(tempV));
        % Loop over files
        for nFiles = 1:length(tempV)
            temp = ica_fuse_read_vols(tempV(nFiles));
            temp(isnan(temp)) = 0;
            temp = temp(ind);
            currentData(:, nFiles) = temp(:);
        end
        % End loop over files
        xAxis = 0;
    else
        %currentData = featureInfo(nFeature).data;
        currentData = getCurrentData(featureInfo(nFeature));
        % get the y axis from ascii file
        yAxis = squeeze(currentData(:, 2, :));
        xAxis = squeeze(currentData(:, 1, 1));
        if length(ind) > length(yAxis)
            error(['Check the mask for ', featureInfo(nFeature).modality, ' data']);
        end

        xAxis = xAxis(ind, :);
        yAxis = yAxis(ind, :);
        clear currentData;
        currentData = yAxis;
        clear yAxis;
    end

    currentData = currentData';

    dataN(nFeature).data = currentData;
    dataN(nFeature).xAxis = xAxis;
    dataN(nFeature).files = featureInfo(nFeature).files;
    dataN(nFeature).feature_name = featureInfo(nFeature).feature_name;
    clear xAxis;
    clear currentData;
end

function [currentData, size_data] = getCurrentData(featureInfo)
%% Load feature data

% Load data
currentData = ica_fuse_loadData(featureInfo.files);

size_data = size(currentData);

if length(size_data) == 3
    size_data(4) = 1;
end
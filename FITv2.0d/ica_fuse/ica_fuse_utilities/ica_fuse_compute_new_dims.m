function new_dims = ica_fuse_compute_new_dims(dims, modalities, mask_ind)
%% Recompute new dimensions
%
% Inputs:
% 1. dims - Dimensions of each feature
% 2. modalities - Modalities
% 3. mask_ind - Mask indices for each feature in a data structure
%
% Outputs:
% 1. new_dims - New dimensions
%

%% Initialise new dimensions
new_dims = dims;

%% Voxels is assumed to be the maximum of all
voxels = max(dims);

%% Convert to cell array
if (~iscell(modalities))
    modalities = cellstr(modalities);
end

modalities = lower(modalities);

%% Check EEG modality
checkEEG = strmatch('eeg', modalities, 'exact');
if (~isempty(checkEEG))
    checkEEG = checkEEG(:)';
    %% Interpolate the EEG Data
    for nn = checkEEG
        % EEG Data length
        eegDataLength = dims(nn);
        if (eegDataLength < voxels)
            %% Resample data
            new_dims(nn) = length(ica_fuse_resample(mask_ind(nn).ind, voxels, eegDataLength));
        end
    end
    %% End for interpolating the EEG data
end
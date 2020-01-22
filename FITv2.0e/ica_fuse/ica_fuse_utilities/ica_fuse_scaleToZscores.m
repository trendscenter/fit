function [data] = ica_fuse_scaleToZscores(modality, data)
% Converts data to z-scores depending upon the modality
% If modality is fMRI, sMRI then only the data points not zero are taken
% into account.
%
% Inputs:
% 1. modality - Options are 'fmri', 'eeg' or 'smri'
% 2. data - Data is a 2D matrix such that rows are observations.
%
% Output:
% scaledData - data converted to z-scores

if ~exist('modality', 'var')
    error('Modality variable must be present. Options are fmri, smri or eeg.');
end

if ~exist('data', 'var')
    error('Data variable must be present. Data must be a 2D matrix where rows are observations.');
end

[data] = ica_fuse_convertToZScores(data);

% if strcmpi(modality, 'eeg')
%     % convert eeg data to z-scores
%     [data] = ica_fuse_convertToZScores(data);
% else
%     % treat data as images
%     [data] = ica_fuse_convertImageToZScores(data);
% end


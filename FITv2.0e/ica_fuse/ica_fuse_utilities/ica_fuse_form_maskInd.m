function [mask_ind] = ica_fuse_form_maskInd(mask_ind, dataInfo)
% Return mask_ind as structure

ica_fuse_defaults;
global EEG_DATA_INDICES;

if ~isstruct(mask_ind)
    % Number of features
    numFeatures = length(dataInfo(1).feature);

    temp = mask_ind;
    clear mask_ind;
    for nn = 1:numFeatures
        if strcmpi(dataInfo(1).feature(nn).modality, 'eeg')
            mask_ind(nn).ind = EEG_DATA_INDICES;
        else
            [diffTimePoints, extns, dims] = ica_fuse_get_countTimePoints(dataInfo(1).feature(nn).files(1).name);
            dat = zeros(dims);
            dat(temp) = 1;
            mask_ind(nn).ind = (dat ~= 0);
            clear dat;
        end
    end

end
function [currentData, normParameter] = ica_fuse_normalize_data(currentData, normalizeVal)
% Normalize data

% Add the normalization schemes
if normalizeVal == 1
    % Default normalization scheme
    meanData = mean(currentData(:).^2);
    normParameter = sqrt(meanData);
elseif normalizeVal == 2
    normParameter = norm(currentData(:), 2);
elseif normalizeVal == 3
    normParameter = std(currentData(:));
else
    normParameter = 1;
    % Add other schemes

end


currentData = (currentData) / normParameter;
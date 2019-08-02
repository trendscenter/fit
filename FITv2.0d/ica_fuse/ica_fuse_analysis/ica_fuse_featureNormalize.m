function [dataN, normalizationPara] = ica_fuse_featureNormalize(dataN, normalizeVal)
%% Normalize data
%
% Inputs:
% 1. dataN - Data
% 2. normalizeVal - Type of normalization
%
% Output
% 1. dataN - Normalized data
% 2. normalizationPara - Normalization parameters

disp('Doing feature normalization ...');
fprintf('\n');
if normalizeVal == 1
    disp('Using square root of mean of squared data to normalize features ...');
elseif normalizeVal == 2
    disp('Using norm 2 of data to normalize features ...');
elseif normalizeVal == 3
    disp('Using standard deviation of data to normalize features ...');
else
    disp('Features are not normalized ...');
end
fprintf('\n');

normalizationPara = zeros(1, length(dataN));

%% Loop over features
for nF = 1:length(dataN)
    currentData = dataN(nF).data;
    [currentData, normalizationPara(nF)] = ica_fuse_normalize_data(currentData, normalizeVal);
    dataN(nF).data = currentData;
    clear currentData;
end
%% End loop over features

disp('Done feature normalization');
fprintf('\n');
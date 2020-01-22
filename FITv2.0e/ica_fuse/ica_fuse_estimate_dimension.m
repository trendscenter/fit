function [b1, mdlVal, aicVal] = ica_fuse_estimate_dimension(data, fwhm)
% Purpose: esimates data dimensionality using MDL criteria
%
% Inputs:
% 1. files - character array or the masked data
% 2. fwhm - smoothness factor (vector like [8 8 11])
%
% Output:
% 1. b1 - Estimated components
% 2. mdlVal - MDL vector of dimension (time points - 1)
% 3. aicVal - AIC vector of dimension (time points - 1)



if ~exist('data', 'var')
    data = ica_fuse_selectEntry('typeSelection', 'multiple', 'filter', '*.img; *.nii', 'title', 'Select data files');
end

if ischar(data)
    
    files = data; clear data;
    
    disp('Loading data ...');
    % load data
    [data] = ica_fuse_loadData(files);   
    
    size_data = size(data);
    
    % make fwhm empty if 
    if ~exist('fwhm', 'var')
        fwhm = [];
    else
        if length(fwhm) == 1
            fwhm = repmat(fwhm, 1, 3);
        end   
    end
    
    data = reshape(data, [prod(size_data(1:3)), size_data(4)]);
    
    % mean for data
    meanData = mean(data, 2);
    
    % convert data to 3D array
    meanData = reshape(meanData, size_data(1:3));
    
    % mean data mask_ind
    mask_ind = find(meanData(:) >= mean(meanData(:)));
    
    data = data(mask_ind, :);
    
end

if length(fwhm) == 1
    fwhm = repmat(fwhm, 1, 3);
elseif length(fwhm) == 2
    fwhm(3) = 1;
elseif length(fwhm) > 3
    error(['Check fwhm vector']);
end

Xdim = size(data, 1) / fwhm(1) / fwhm(2) / fwhm(3);


tdim = size(data,2);

% do PCA
[V, Lambda, tmp, lam] = ica_fuse_v_pca(data, 1, tdim, 0, 'untranspose', 'no');

clear data; clear V; clear Lambda; clear tmp;

aicVal = zeros(1, tdim - 1);
mdlVal = aicVal;

for numc = 1:(tdim-1)
    tempLam = lam(numc + 1:end);
    % calculating nthroot of the eigen values
    storeLam = tempLam.^(1/(tdim - numc));
    L(numc) = log((tdim-numc)*prod(storeLam)) - log(sum(tempLam));
    %L(numc) = log(tdim - numc) + (sum(log(tempLam)) / (tdim - numc)) -  log(sum(tempLam));

    geomMean(numc) = (tdim - numc)*prod(storeLam);
    arithmeticMean(numc) = sum(tempLam);
    clear storeLam;

    aicVal(numc) = -2*Xdim*(tdim-numc)*L(numc) + 2*numc*(2*tdim-numc);
    mdlVal(numc) = -Xdim*(tdim-numc)*L(numc) + 0.5*numc*(2*tdim-numc)*log(Xdim);
    kicVal(numc) = -2*Xdim*(tdim - numc)*L(numc) + 3*numc*(2*tdim - numc);
end

% return the number of estimated components using different criterias
[a1, b1] = min(mdlVal);

[a2, b2] = min(aicVal);
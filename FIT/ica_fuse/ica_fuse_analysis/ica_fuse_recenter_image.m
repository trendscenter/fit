function [v, offset] = ica_fuse_recenter_image(v)
%% Center image distribution to zero (Based on code from Elena Allen)
%
% Inputs:
% v - Input image
%
% Outputs:
% v - Centered image
%

SMOOTHINGVALUE = 1.1;

mask_ind = (v ~= 0);
v2 = v(mask_ind);

if (isempty(mask_ind))
    error('No non-zero voxels found');
end

v2 = v2(:);
nPoints = ceil(length(v2)/10);
hist_bins = linspace(min(v2), max(v2), nPoints);
%hist_bins = min(v2):.1:max(v2);
[dist,bins] = hist(v2, hist_bins);
dist = ica_fuse_gauss_smooth1D(dist, SMOOTHINGVALUE);
[mv, mind] = max(dist);
offset = bins(mind);
v(mask_ind) = v2-offset;


function data = ica_fuse_gauss_smooth1D(data, FWHM)
% function data = gauss_smooth(data,FWHM)
% smooth in xy plane with a Gaussian Kernel
% replaces input image with smoothed image
%
% Input:
% 1. data - data can be double or object of class complex_data
% and is of size nrows by ncols where each column indicates the
% observation.
% 2. FWHM - Smoothing factor
%
% Output:
% data - smoothed data

if (numel(data) == length(data))
    data = data(:);
end

% check the data
if isreal(data)
    % check if data is of real type
    [data] = smoothData(data, FWHM);
elseif isa(data, 'complex_data')
    % detect if data is of complex data class
    % get the field names
    fieldNames = fieldnames(data);
    
    % set the fields to the object
    for ii = 1:length(fieldNames)
        % get the field data
        currentData = smoothData(getfield(data, fieldNames{ii}), FWHM);
        % set the data to the object
        data = setfield(data, fieldNames{ii}, currentData);
    end
else
    error('Unknown data type');
end
% end for checking

function [dataIn] = smoothData(dataIn, FWHM)
% sub function to smooth the data

% loop over cols
for ii = 1:size(dataIn, 2)
    %data(:, ii) = gauss_smooth1D(data(:, ii), FWHM);
    data = dataIn(:, ii);
    % data = gauss_smooth(data,2);
    tdim = length(data);
    %create Gaussian Kernel
    s  = [FWHM];
    
    s  = s/sqrt(8*log(2));				% FWHM -> Gaussian parameter
    
    x  = round(6*s(1));
    kdimx = 2*x+1;
    x = [-x:x];
    x  = exp(-(x).^2/(2*(s(1)).^2));
    x  = x/sum(x);
    
    gauss = zeros(kdimx,1);
    for i = 1:kdimx,
        gauss(i) = sqrt(x(i));
    end;
    gauss = gauss/max(gauss);
    gauss = gauss/sum(gauss);
    [val ind] = max(gauss);
    %ind = round(ind/2);
    data = mean(data)+conv(data-mean(data),gauss);
    %data = data(ind-1:tdim+ind-1);
    data = data(ind:tdim+ind-1);
    dataIn(:, ii) = data;
    clear data;
end
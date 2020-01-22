function [y] = ica_fuse_resample(x, p, q)
% Resample data using resample function. 
%
% Input:
% 1. x - data to be sampled 
% 2. p - factor to interpolate
% 3. q - factor to decimate
%
% Output:
% y - resampled data

interpFactor = ceil(p/q);

if interpFactor > 1
    y = ica_fuse_interp(x, interpFactor);
    %y = resample(x, p, q);
else
    y = x;
end

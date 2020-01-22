function [scaledData, minVal, maxVal] = ica_fuse_scaleImage(tmin, tmax, imageData, imType, imageValues, imRange)
% Scales image depending upon functional or anatomical
%
% Input:
% 1. tmin - colorbar minimum
% 2. tmax - colorbar maximum
% 3. imageData - image data
% 4. imType - 'functional' or 'anatomical'
%
% Output:
% scaledData - scaled image


if (exist('imRange', 'var') && (length(imRange) > 1))
    minVal = min(imRange);
    maxVal = max(imRange);
else
    % Get the actual min and max of image data
    maxVal = max(imageData(:));
    minVal = min(imageData(:));
end

if isempty(maxVal)
    maxVal = 0;
end

if isempty(minVal)
    minVal = 0;
end

% For functional images set the minimum and maximum within the colorbar
% values
if strcmpi(imType, 'functional')
    
    if imageValues == 1
        if maxVal < abs(minVal)
            maxVal = abs(minVal);
        else
            minVal = -maxVal;
        end
    elseif imageValues == 2
        minVal = 0;
    elseif imageValues == 3
        minVal = 0;
    else
        maxVal = 0;
    end
    
end

rangeVal = maxVal - minVal;
trange = tmax - tmin;

if rangeVal ~= 0
    scaledData = (((imageData-minVal)./rangeVal)./(1/trange))+tmin;
else
    scaledData = imageData;
end

% This will be required to plot the text of colorbar
minVal = round(10*minVal)/10;
maxVal = round(10*maxVal)/10;
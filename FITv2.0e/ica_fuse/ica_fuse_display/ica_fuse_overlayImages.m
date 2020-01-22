function [im, minICAIm, maxICAIm, minInterval, maxInterval] = ica_fuse_overlayImages(func_im, anat, structDIM, imageValues, imRange)
% Overlay images
% Scale structural images to be higher than the functional images
%


% Initialise min interval and max interval
minInterval = 0;
maxInterval = 100;

% Number of images
numImages = size(func_im, 1);

% minimum and maximum for structural images
newMin = numImages*maxInterval + 1;
newMax = (numImages + 1)*maxInterval;

% scale structural higher
anat = ica_fuse_scaleImage(newMin, newMax, anat(:), 'anatomical');

% Initialise im with anatomical image
im = anat;

%%%%%%%%%%% Loop over functional images %%%%%%%%%%%%%%%%%
for nImages = 1:numImages
    % get the current image
    currentIm = func_im(nImages, :, :, :);
    % Make the current im column vector
    currentIm = currentIm(:);

    unitColor = ica_fuse_range(currentIm)/(min([256 64*(numImages + 1)])/(numImages + 1));

    % Update new minimum and maximum
    newMin = (nImages - 1)*maxInterval + minInterval;
    newMax = nImages*maxInterval - unitColor;

    % Get the indices that are not zero
    indices(nImages).ind = find(currentIm ~= 0);

    % scale functional image
    [currentIm, minICAIm, maxICAIm] = ica_fuse_scaleImage(newMin, newMax, currentIm(:), 'functional', imageValues, imRange);

    currentIm = currentIm(indices(nImages).ind);

    % Overlay images
    if ~isempty(indices(nImages).ind)
        im(indices(nImages).ind) = currentIm;
    end

    clear currentIm;

end
%%%%%%%%%%% End loop over functional Images %%%%%%%%%%%%%%%%%

% reshape the anatomical image to 3D
im = reshape(im, structDIM);


function [minVal, maxVal] = getMinMax(imageData, imageValues)

% Get the actual min and max of image data
maxVal = max(imageData(:));
minVal = min(imageData(:));

if isempty(maxVal)
    maxVal = 0;
end

if isempty(minVal)
    minVal = 0;
end

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
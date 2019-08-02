function [funcImg, minICAIm, maxICAIm, minInterval, maxInterval] = ica_fuse_create_composite(anat, funcData, ...
    imageValues, anatomicalView, structDIM, imRange)
% creates a composite image from structural and functional image
%
% Input:
% 1. anat - anatomical image
% 2. funcData 
% 3. imageValues - 1 means positive and negative, 2 means positive, 3 means
% Absoulte, 4 means Negative
% 4. anatomicalView - Anatomical view
%
% Output:
% 1. funcImg - functional image converted to montage
% 2. minICAIm - Minimum value of functional image
% 3. maxICAIm - Max value of functional image
% 4. minInterval - Minimum interval
% 5. maxInterval - Max interval

% Initialise all vars
maxICAIm = zeros(1, size(funcData, 1));
minICAIm = zeros(1, size(funcData, 1));

%--get image in correct plane
if ica_fuse_findstr(lower(anatomicalView), 'sagittal')
    % sagittal plane
    permuteOrder = [2 3 1];
elseif ica_fuse_findstr(lower(anatomicalView), 'coronal')
    % coronal plane
    permuteOrder = [1 3 2];
else
    permuteOrder = [1 2 3];
end

% Loop over number of functional images
for nImages = 1:size(funcData, 1)
    status = 0;
    % Get the current image
    func_im = funcData(nImages, :, :, :);  
    
    % Overlay images
    [im, minICAIm(nImages), maxICAIm(nImages), minInterval, maxInterval] = ica_fuse_overlayImages(func_im, anat, structDIM, imageValues, imRange);                   
    
    im = reshape(im, structDIM);
 
    % Permute the images
    im = permute(im, permuteOrder);     
    DIM = structDIM(permuteOrder);
    % Get the montage
    [montage_im] = ica_fuse_return_montage(im, DIM);
    clear im;
    % store in the output cell array
    funcImg(nImages, :, :) = montage_im;
    clear montage_im;
end
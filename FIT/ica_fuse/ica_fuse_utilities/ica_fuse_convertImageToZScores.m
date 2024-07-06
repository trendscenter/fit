function [im] = ica_fuse_convertImageToZScores(icasig)
% converts images to z scores

% Number of images
numberOfImages = size(icasig,1);

icasig(isnan(icasig)) = 0;

% loop over images
for i=1:numberOfImages

    % get the current image
    v = icasig(i, :);    
    mask_ind = (v ~= 0);
    v(mask_ind) = detrend(v(mask_ind), 0);
    v2 = v(mask_ind);
    
    vstd = std(v2);
    if vstd ~= 0
        v=v./(eps + vstd);
    else
        disp('Not converting to z-scores as division by zero warning may occur.');
    end

    icasig(i,:) = v;
end

im = icasig;



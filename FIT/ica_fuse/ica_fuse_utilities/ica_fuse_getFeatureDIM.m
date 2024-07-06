function [dims, voxels] = ica_fuse_getFeatureDIM(mask_ind)
% Return the dimensions of each feature

dims = zeros(1, length(mask_ind));
for nFeatures = 1:length(mask_ind)
    ind = find(mask_ind(nFeatures).ind ~= 0);
    dims(nFeatures) = length(ind);
end

% Number of voxels
voxels = max(dims);
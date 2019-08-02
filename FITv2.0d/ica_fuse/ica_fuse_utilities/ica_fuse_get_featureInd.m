function [featureStartInd, featureEndInd] = ica_fuse_get_featureInd(selFeatures, dataLength)
% Get feature indices

if selFeatures(1) == 1
    featureStartInd = 1;
else
    featureStartInd = sum(dataLength(1:selFeatures(1)-1)) + 1;
end

featureEndInd = sum(dataLength(1:selFeatures(1)));
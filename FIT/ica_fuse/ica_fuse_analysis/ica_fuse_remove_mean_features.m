function [whitesig] = ica_fuse_remove_mean_features(whitesig, featureDataLength)

for nn = 1:length(featureDataLength)
    [featureStartInd, featureEndInd] = ica_fuse_get_featureInd(nn, featureDataLength);
    temp = whitesig(:, featureStartInd:featureEndInd);
    temp = temp';
    temp = detrend(temp, 0);
    whitesig(:, featureStartInd:featureEndInd) = temp';
    clear temp
end
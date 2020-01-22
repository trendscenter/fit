function [V, Lambda, whitesig, whiteM, dewhiteM] = ica_fuse_calculate_pca(dat, numc, type_pca, reference)
% Calculates PCA and does whitening to the data
%
% Input:
% 1. dat - Data
% 2. numc - Number of components
%
% Output:
% 1. V - Eigen vectors
% 2. Lambda - Eigen values in a matrix
% 3. whitesig - Whitened data
% 4. whiteM - Whitening matrix
% 5. dewhiteM - Dewhitening matrix

ica_fuse_defaults;
global RM_PCA;

remove_mean = 'no';
if ~isempty(RM_PCA)
    if RM_PCA
        remove_mean = 'yes';
    end
end

if ~exist('type_pca', 'var')
    type_pca = 'standard';
end

if (~exist('reference', 'var'))
    reference = [];
end

if (strcmpi(remove_mean, 'yes'))
    dat = ica_fuse_remove_mean(dat');
    dat = dat';
end


disp('Doing data reduction ...');

fprintf('\n');
if strcmpi(type_pca, 'standard')
    disp('Running standard pca ...');
    % Calculate PCA
    if (size(dat, 1) > size(dat, 2))
        [V, Lambda] = ica_fuse_svd(dat', numc);
    else
        [V, Lambda] = ica_fuse_v_pca(dat, 1, numc, 0, 'transpose', remove_mean);
    end
else
    disp('Running pca with reference ...');
    % Calculate PCA with reference
    [V, Lambda] = ica_fuse_pca_reference(dat, 1, numc, reference, 0, 'transpose', remove_mean);
end
fprintf('\n');
% Whitening
[whitesig, whiteM, dewhiteM] = ica_fuse_v_whiten(dat, V, Lambda, 'transpose');

disp('Done with data reduction');
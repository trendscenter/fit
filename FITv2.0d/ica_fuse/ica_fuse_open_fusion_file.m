function fusionFile = ica_fuse_open_fusion_file


ica_fuse_defaults;
global FUSION_INFO_MAT_FILE;

fusionFile = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
    ['*', FUSION_INFO_MAT_FILE, '*.mat'], 'title', 'Select joint ICA fusion information file');

load(fusionFile);

if ~exist('fusionInfo', 'var')
    error(['Selected file: ', fusionFile, ' is not a valid joint ICA fusion information file']);
end

% open fusion file
fusionFile = ica_fuse_setup_analysis([], fusionFile);
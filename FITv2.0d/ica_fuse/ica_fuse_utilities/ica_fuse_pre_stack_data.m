function stackInfo = ica_fuse_pre_stack_data(fusionInfo)
%% Purpose: Before loading data, compute mask and dims.
%
% Inputs:
% 1. fusionInfo - Fusion Info data structure
%
% Outputs:
% stackInfo - stackInfo data structure containing necessary fields
%

%% Create mask
stackInfo.mask_ind = ica_fuse_createMask(fusionInfo.setup_analysis.dataInfo, fusionInfo.setup_analysis.maskFile);

%% Get sample size for each feature. Voxels correspond to maximum data size
[stackInfo.dims, stackInfo.voxels] = ica_fuse_getFeatureDIM(stackInfo.mask_ind);

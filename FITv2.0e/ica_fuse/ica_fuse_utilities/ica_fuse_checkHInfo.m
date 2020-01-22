function [HInfo] = ica_fuse_checkHInfo(HInfo, fileName, extn)
% Check header information. HInfo contains V field that will be used for
% writing images. Depending upon the V field, file name of functional data
% and image extension to be written return HInfo.
% 
% Input:
% 1. HInfo - Structure containing volume information
% 2. File name - File name of the functional image from which structure was built.
% 3. extn - Image extension with which the component images will be
%    written.
%
% Output:
% HInfo - Structure containing volume information

% Get the first volume
V = HInfo.V(1);

% Check image extension
if strcmpi(extn, '.nii')
    
    if strcmpi(class(V(1).private), 'struct')
        clear V;
        % create nifti structure
        V = ica_fuse_get_vol_nifti(deblank(fileName));
        HInfo = rmfield(HInfo, 'V');
        HInfo.V(1) = V(1);
        clear V;
    end
    
else
    
    if strcmpi(class(V(1).private), 'ica_fuse_nifti')
        clear V;
        % create analyze structure
        V = ica_fuse_getVol(deblank(fileName), 1);
        HInfo = rmfield(HInfo, 'V');
        HInfo.V(1) = V(1);
        clear V;
    end
    
end
%%%%%%%%% End for checking HInfo %%%%%%%%%

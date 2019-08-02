function V = ica_fuse_get_vol_nifti(file)
% get volume of file by parsing the extension

% Parse extension
[file, fileNum] = ica_fuse_parseExtn(file);

V = ica_fuse_spm_vol_nifti(file, fileNum);
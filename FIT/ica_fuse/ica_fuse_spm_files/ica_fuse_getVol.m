function V = ica_fuse_getVol(files, number)
%% Get volume
%

files = ica_fuse_rename_4d_file(files);

if (exist('number', 'var'))
    if (max(number) > size(files, 1))
        error('Image number requested exceeds the number of images');
    end
    files = deblank(files(number, :));
end

V = ica_fuse_spm_vol(files);

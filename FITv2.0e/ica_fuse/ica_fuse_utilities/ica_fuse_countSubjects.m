function [numSubjects] = ica_fuse_countSubjects(dataInfo)
% Count number of subjects per group

numSubjects = zeros(1, length(dataInfo));
for nGroup = 1:length(dataInfo)
    tmpFile = char(dataInfo(nGroup).feature(1).files.name);
    [~, pp, extn] = fileparts(icatb_parseExtn(deblank(tmpFile(1, :))));
    if (strcmpi(extn, '.img') || strcmpi(extn, '.nii'))
        temp = ica_fuse_rename_4d_file(tmpFile);
        numSubjects(nGroup) = size(temp, 1);
    else
        temp = ica_fuse_loadData(tmpFile);
        siz_tmp = length(size(temp));
        numSubjects(nGroup) = size(temp, siz_tmp);
    end
    
end
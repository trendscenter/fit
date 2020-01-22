function [numSubjects] = ica_fuse_countSubjects(dataInfo)
% Count number of subjects per group

numSubjects = zeros(1, length(dataInfo));
for nGroup = 1:length(dataInfo)
    temp = ica_fuse_rename_4d_file(str2mat(dataInfo(nGroup).feature(1).files.name));
    numSubjects(nGroup) = size(temp, 1);
end
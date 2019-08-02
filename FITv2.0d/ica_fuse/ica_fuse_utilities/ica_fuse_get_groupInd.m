function [groupInd] = ica_fuse_get_groupInd(selGroups, numSubjects)
% Get group indices

if selGroups(1) == 1
    groupInd = (1:numSubjects(1));
else
    groupInd = (sum(numSubjects(1:selGroups(1)-1)) + 1) : (sum(numSubjects(1:selGroups(1))));
end
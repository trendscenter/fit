function [newPrefix, indices_prefix, originalPrefix] = ica_fuse_classifyFiles(files)
% Return summary string from the files where files are in character format.

[prefix, indices_prefix] = ica_fuse_spm_get('fileSummary', files);
originalPrefix = prefix;

% Show the number of entries during the prefix mode
for ii = 1:size(prefix, 1)
    temp = deblank(prefix(ii, :));
    findNumberFiles = find(indices_prefix(ii, :) ~= 0);
    temp = ['[1:', num2str(length(findNumberFiles)), '] ', temp];
    prefixStruct(ii).prefix = temp;    
end
clear prefix;
newPrefix = str2mat(prefixStruct.prefix); % Updated prefix
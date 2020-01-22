function subFolders = ica_fuse_listSubFolders(sourceDir)
% list sub folders excluding the current folder and parent folder

% list sub folders
subFolders = ica_fuse_listDir(sourceDir);
tempVec = [];

index1 = strmatch('.', subFolders, 'exact');
index2 = strmatch('..', subFolders, 'exact');

index = [index1 index2];

CheckOne = ones(size(subFolders, 1), 1);

CheckOne(index) = 0;

tempVec = find(CheckOne ~= 0);

if ~isempty(tempVec)
    subFolders = subFolders(tempVec, :);
    for ii = 1:size(subFolders, 1)
        store(ii).name = deblank(fullfile(sourceDir, subFolders(ii, :)));
    end

    clear  subFolders;
    subFolders = str2mat(store.name);
    clear store;
else
    subFolders = [];
end




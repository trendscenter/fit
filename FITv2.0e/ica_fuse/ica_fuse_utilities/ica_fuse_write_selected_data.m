function ica_fuse_write_selected_data(dataInfo, file_name)
% record information in a file

disp(['Saving data information in file ', file_name]);
fid = fopen(file_name, 'w');

if fid == -1
    error(['File:', file_name, 'cannot be opened']);
end

numGroups = length(dataInfo);

numfeatures = length(dataInfo(1).feature);

% loop over groups
for nGroup = 1:numGroups
    for nfeature = 1:numfeatures        
        titleToPrint = ['Selected data files for ', dataInfo(nGroup).name, ' ', dataInfo(nGroup).feature(nfeature).name, ...
                ' feature are as follows:'];        
        % print the string to a file
        ica_fuse_printString(fid, titleToPrint, str2mat(dataInfo(nGroup).feature(nfeature).files.name));
        if (nGroup*nfeature) ~= (numGroups*numfeatures)
            fprintf(fid, '\n\n');
        end
        
    end
end

fclose(fid);             

disp(['Done saving data information']);

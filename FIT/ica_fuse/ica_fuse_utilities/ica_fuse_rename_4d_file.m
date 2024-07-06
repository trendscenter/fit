function files = ica_fuse_rename_4d_file(files, fileNumber)
%For nifti and analyze files rename the files by adding a number at the
%end

if ~exist('fileNumber', 'var')
    fileNumber = [];
end

if ~isempty(files)

    % Do pattern match to check IMG or NII files
    files = cellstr(files);
    
    % MATCH IMG OR NII
    checkNII = regexpi(files, '\.nii$');
    checkIMG = regexpi(files, '\.img$');
    
    good_inds1 = ica_fuse_good_cells(checkNII);
    good_inds2 = ica_fuse_good_cells(checkIMG);
    
    % Get good cells
    good_inds = (good_inds1 | good_inds2);
    good_inds = find(good_inds ~= 0);
    
    clear good_inds1 good_inds2 checkNII checkIMG;

    % If IMG or NII files exist check the headers of the image files
    if ~isempty(good_inds)
        imFiles = files(good_inds);
        if (~isempty(fileNumber)) && (length(fileNumber) == 1) && (fileNumber == 1)
            %imFiles = [str2mat(imFiles), repmat(',1', length(imFiles), 1)];
            for nIm = 1:length(imFiles)
                imFiles{nIm} = [imFiles{nIm}, ',1'];
            end
            files(good_inds) = imFiles;
        else
            % Loop over img or nii files
            for nIm = 1:length(imFiles)
                currentFile = imFiles{nIm};

                try
                    ni = ica_fuse_nifti(currentFile);
                    numFiles = ni.dat.dim(4);
                catch
                    numFiles = 1;
                end

                if (isempty(fileNumber))
                    tempFileNum = (1:numFiles);
                else
                    tempFileNum = fileNumber;
                end

                tempFileNum(tempFileNum > numFiles) = [];

                if ~isempty(tempFileNum)
                    tempFiles = [repmat(currentFile, length(tempFileNum), 1), repmat(',', length(tempFileNum), 1), num2str(tempFileNum(:))];
                else
                    tempFiles = '';
                end
                files{good_inds(nIm)} = tempFiles;
            end
            % End loop over img or nii files
            files = cellstr(str2mat(files));
            ind = ica_fuse_good_cells(files);
            files = files(ind);
        end
    end
    % End for checking IMG or NII files

    if ~isempty(files)
        files = str2mat(files);
    end

end
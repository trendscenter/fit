function  [outputFiles] = ica_fuse_write_feature_data(icasig, A, groups_icasig, dataFiles, mask_ind, outFile, ...
    outputDir, timeAxis, numSubjects)
% Write corresponding feature images for joint components. The output file
% naming will depend on the global variable JOINT_COMPONENT_NAMING (See
% ica_fuse_defaults.m). This will be added to the fusionInfo variable.
%
% Input:
% 1. icasig - Joint components
% 2. dataFiles - Data files
% 3. mask_ind - Feature mask indices
% 4. outFile - Output file naming
% 5. outputDir - Output Directory
% 6. timeAxis - Time axis
%
% Output:
% outputFiles - Output file names

% Load defaults
ica_fuse_defaults;

global MRI_DATA_FILTER;
global EEG_DATA_FILTER;

% size of icasig
[numComp, data_points] = size(icasig);

numGroups = size(groups_icasig, 3);

% Number of features
numFeatures = length(dataFiles);

% Loop over features
for nFeature = 1:numFeatures
    current_dataFile = deblank(dataFiles(nFeature).name(1, :));
    [pathstr, currentFile, extn] = fileparts(ica_fuse_parseExtn(current_dataFile));

    if (~strcmpi(extn, '.img')) & (~strcmpi(extn, '.nii'))
        % For EEG Data
        xAxis = timeAxis(nFeature).data;
        featureData(nFeature).dim = length(xAxis);
        featureData(nFeature).V = [];
        featureData(nFeature).xAxis = xAxis;
    else
        % For MRI data
        if strcmpi(extn, '.img') & strcmpi(MRI_DATA_FILTER(2:end), '.img')
            V = ica_fuse_getVol(current_dataFile, 1);
        else
            V = ica_fuse_get_vol_nifti(current_dataFile);
        end
        featureData(nFeature).dim = length(find(mask_ind(nFeature).ind ~= 0));
        featureData(nFeature).V = V;
        featureData(nFeature).xAxis = [];
        clear V;
    end
end
% End loop over features

startLength = 1; endLength = 0;

% Loop over features
for nFeature = 1:numFeatures

    endLength = endLength + featureData(nFeature).dim;
    groupFiles = [];

    % loop over components
    for nComp = 1:numComp

        currentData = icasig(nComp, startLength:endLength);

        [fileIndex] = ica_fuse_returnFileIndex(nComp);

        % write ascii data
        if isempty(featureData(nFeature).V)

            currentData = currentData';

            eegDataLength = length(mask_ind(nFeature).ind);
            downSampleFactor = ceil(length(currentData)/eegDataLength);

            fileToWrite = [outFile, 'feature_', num2str(nFeature), '_', fileIndex, EEG_DATA_FILTER(2:end)];
            % Downsample EEG data
            [currentData, timeAxis] = ica_fuse_downsample(currentData, featureData(nFeature).xAxis, ...
                downSampleFactor);
            % Make the data with two columns
            data = zeros(length(currentData), 2);
            data(:, 1) = timeAxis;
            data(:, 2) = currentData;

            clear currentData;
            %%% end for making data to be n by 2
            save(fullfile(outputDir, fileToWrite), 'data', '-ascii');
            clear data;clear timeAxis;

            compFiles(nComp).name = fileToWrite;

            if ~isempty(groups_icasig)
                for nGroup = 1:numGroups

                    groupFileToWrite = [outFile, 'feature_', num2str(nFeature), '_group_', num2str(nGroup), '_', ...
                        fileIndex, EEG_DATA_FILTER(2:end)];

                    % current group data of the feature
                    group_data = squeeze(groups_icasig(startLength:endLength, nComp, nGroup));

                    % Downsample EEG data
                    [group_data, timeAxis] = ica_fuse_downsample(group_data, featureData(nFeature).xAxis, ...
                        downSampleFactor);

                    data = zeros(length(group_data), 2);
                    data(:, 1) = timeAxis;
                    data(:, 2) = group_data;
                    save(fullfile(outputDir, groupFileToWrite), 'data', '-ascii');
                    clear data; clear timeAxis;
                    groupFiles(nGroup).comp(nComp).name = groupFileToWrite;
                end
            end


        else
            % write image data
            if strcmpi(MRI_DATA_FILTER(2:end), '.img')
                fileToWrite = [outFile, 'feature_', num2str(nFeature), '_', fileIndex, MRI_DATA_FILTER(2:end)];
            else
                fileToWrite = [outFile, 'feature_', num2str(nFeature), MRI_DATA_FILTER(2:end)];
            end

            newData = zeros(1, prod(featureData(nFeature).V(1).dim(1:3)));

            ind = mask_ind(nFeature).ind;

            ind = find(ind ~= 0);

            newData(ind) = currentData;

            newData = reshape(newData, featureData(nFeature).V(1).dim(1:3));

            compV = featureData(nFeature).V(1);

            compV.descrip = ['Feature ', dataFiles(nFeature).feature_name];

            compV.fname = fullfile(outputDir, fileToWrite);

            if strcmpi(MRI_DATA_FILTER(2:end), '.img')
                compV.n(1) = 1;
            else
                compV.n(1) = nComp;
            end

            ica_fuse_write_vol(compV, newData);

            clear newData;

            if strcmpi(MRI_DATA_FILTER(2:end), '.img')
                compFiles(nComp).name = fileToWrite;
            else
                compFiles(1).name = fileToWrite;
            end


            if ~isempty(groups_icasig)

                for nGroup = 1:numGroups

                    if strcmpi(MRI_DATA_FILTER(2:end), '.img')
                        groupFileToWrite = [outFile, 'feature_', num2str(nFeature), '_group_', num2str(nGroup), '_', ...
                            fileIndex, MRI_DATA_FILTER(2:end)];
                    else
                        groupFileToWrite = [outFile, 'feature_', num2str(nFeature), '_group_', num2str(nGroup), ...
                            MRI_DATA_FILTER(2:end)];
                    end

                    % current group data of the feature
                    newData = zeros(1, prod(featureData(nFeature).V(1).dim(1:3)));
                    group_data = squeeze(groups_icasig(startLength:endLength, nComp, nGroup));
                    newData(ind) = group_data;
                    newData = reshape(newData, featureData(nFeature).V(1).dim(1:3));
                    if strcmpi(MRI_DATA_FILTER(2:end), '.img')
                        groupFiles(nGroup).comp(nComp).name = groupFileToWrite;
                    else
                        groupFiles(nGroup).comp(1).name = groupFileToWrite;
                    end

                    compV.fname = fullfile(outputDir, groupFileToWrite);
                    ica_fuse_write_vol(compV, newData);
                    clear newData;
                end

            end

            clear compV;

        end


    end

    % Store file names (Don't store full file)
    outputFiles(nFeature).name = str2mat(compFiles.name);
    outputFiles(nFeature).feature_name = dataFiles(nFeature).feature_name;
    outputFiles(nFeature).groupFiles = groupFiles;
    clear compFiles; clear currentData;
    clear groupFiles;

    startLength = endLength + 1;

end
% end loop over features


clear V;

if ~isempty(A)
    %%% Write Mixing coefficients %%%
    V = ica_fuse_getVol(which('nsingle_subj_T1_2_2_5.nii'));
    V.descrip = ['Number of subjects = ', num2str(numSubjects)];

    if ~iscell(A)
        V.dim(1:3) = [size(A, 1), size(A, 2), 1];
        fileToWrite = [outFile, 'mixing_coeff_.img'];
        fileToWrite = fullfile(outputDir, fileToWrite);
        V.fname = fileToWrite;
        ica_fuse_write_vol(V, A);
    else
        for n = 1:length(A)
            tmpA = A{n};
            V.dim(1:3) = [size(tmpA, 1), size(tmpA, 2), 1];
            fileToWrite = [outFile, 'mixing_coeff_feature_', num2str(n), '.img'];
            fileToWrite = fullfile(outputDir, fileToWrite);
            V.fname = fileToWrite;
            ica_fuse_write_vol(V, tmpA);
        end
    end
    %% End for writing mixing coefficients %%%
end
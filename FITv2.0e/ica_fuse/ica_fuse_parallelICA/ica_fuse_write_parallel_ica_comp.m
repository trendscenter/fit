function outputFiles = ica_fuse_write_parallel_ica_comp(aveComp, loadingCoeff, dataInfo, mask_ind, outFile, outputDir, numSubjects)
% Write parallel ICA components based on global variable
% PARALLEL_ICA_COMPONENT_NAMING.
%
% Input:
% 1. aveComp - Parallel ICA components data. It is a cell array of length
% equal to number of modalities involved.
% 2. loadingCoeff - Loading coefficients of each modality.
% 3. dataInfo - dataInfo structure
% 4. mask_ind - Mask indices
% 5. outFile - Output file naming
% 6. outputDir - Output directory to save the images
% 7. numSubjects - Number of subjects for each group
%
% Output:
% outputFiles - Output file names

% Load defaults
ica_fuse_defaults;

global MRI_DATA_FILTER;
global GENE_DATA_FILTER;
global EEG_DATA_FILTER;
global FLIP_SIGN_COMPONENTS;
global BEH_DATA_FILTER;

% Modalities and feature names
modalities = cellstr(str2mat(dataInfo(1).feature.modality));
featureNames = cellstr(str2mat(dataInfo(1).feature.name));


loadingCoeffV = ica_fuse_getVol(which('nsingle_subj_T1_2_2_5.nii'));

% Loop over feature names
for nF = 1:length(featureNames)
    if strcmpi(modalities{nF}, 'fmri') || strcmpi(modalities{nF}, 'smri')
        % Loop over components
        for nComp = 1:size(aveComp{nF}, 1)
            [fileIndex] = ica_fuse_returnFileIndex(nComp);
            % Current component
            currentData = aveComp{nF}(nComp, :);
            % write image data
            if strcmpi(MRI_DATA_FILTER(2:end), '.img')
                fileToWrite = [outFile, 'feature_', num2str(nF), '_', fileIndex, MRI_DATA_FILTER(2:end)];
            else
                fileToWrite = [outFile, 'feature_', num2str(nF), MRI_DATA_FILTER(2:end)];
            end
            
            % Mask indices
            ind = mask_ind(nF).ind;
            
            newData = zeros(size(ind));
            
            newData(ind) = currentData;
            
            % Get the volume of the first file
            compV = ica_fuse_getVol(deblank(dataInfo(1).feature(nF).files(1).name(1, :)), 1);
            
            compV.descrip = ['Feature ', featureNames{nF}];
            
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
            
        end
        % End loop over components
        
    elseif (strcmpi(modalities{nF}, 'gene') || strcmpi(modalities{nF}, 'behavioral'))
        % Loop over components
        for nComp = 1:size(aveComp{nF}, 1)
            [fileIndex] = ica_fuse_returnFileIndex(nComp);
            % Current component
            currentData = aveComp{nF}(nComp, :);
            currentData = currentData(:);
            if (strcmpi(modalities{nF}, 'gene'))
                fileToWrite = [outFile, 'feature_', num2str(nF), '_', fileIndex, GENE_DATA_FILTER(2:end)];
            else
                fileToWrite = [outFile, 'feature_', num2str(nF), '_', fileIndex, BEH_DATA_FILTER(2:end)];
            end
            save(fullfile(outputDir, fileToWrite), 'currentData', '-ascii');
            compFiles(nComp).name = fileToWrite;
        end
        % End loop over components
        
    elseif strcmpi(modalities{nF}, 'eeg')
        
        % Loop over components
        for nComp = 1:size(aveComp{nF}, 1)
            [fileIndex] = ica_fuse_returnFileIndex(nComp);
            % Current component
            currentData = aveComp{nF}(nComp, :)';
            if (nComp == 1)
                timeAxis = ica_fuse_loadData(deblank(dataInfo(1).feature(nF).files(1).name(1, :)));
                timeAxis = squeeze(timeAxis(:, 1));
                timeAxis = timeAxis(mask_ind(nF).ind);
            end
            currentData = [timeAxis, currentData];
            fileToWrite = [outFile, 'feature_', num2str(nF), '_', fileIndex, EEG_DATA_FILTER(2:end)];
            save(fullfile(outputDir, fileToWrite), 'currentData', '-ascii');
            compFiles(nComp).name = fileToWrite;
        end
        % End loop over components
        
    end
    % End for checking modality
    
    % Store file names (Don't store full file)
    outputFiles(nF).name = str2mat(compFiles.name);
    outputFiles(nF).feature_name = featureNames{nF};
    clear compFiles;
    
    %%% Write Mixing coefficients %%
    loadingCoeffV.dim(1:3) = [size(loadingCoeff{nF}), 1];
    loadingCoeffV.descrip = ['Number of subjects = ', num2str(numSubjects)];
    fileToWrite = [outFile, 'feature_', num2str(nF), '_load_coeff_.img'];
    % Store mixing coeffic
    outputFiles(nF).loadingCoeffFile = fileToWrite;
    
    fileToWrite = fullfile(outputDir, fileToWrite);
    loadingCoeffV.fname = fileToWrite;
    ica_fuse_write_vol(loadingCoeffV, loadingCoeff{nF});
    %% End for writing mixing coefficients %%%
    
end

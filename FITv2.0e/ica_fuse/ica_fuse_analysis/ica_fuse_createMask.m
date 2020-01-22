function maskIndices = ica_fuse_createMask(dataInfo, mask_files)
%% Create mask for each feature separately. For fmri, use intersection of
% voxels from different fmri modalities.
%
% Inputs:
% 1. dataInfo - dataInfo structure.
% 2. mask_files - Mask files for each modality in a cell array
%
% Outputs:
% maskIndices - Mask indices data structure.
%

%% Load defaults
ica_fuse_defaults;

global EEG_DATA_INDICES;

if ~exist('mask_files', 'var')
    mask_files = [];
end

featureInfo = ica_fuse_get_feature_info(dataInfo);

%% Number of features
numFeatures = length(featureInfo);

%% Initialise mask indices
maskIndices = repmat(struct('ind', []), 1, numFeatures);

%% Default mask
if isempty(mask_files)
    
    disp(['Creating mask from the data ...']);
    
    %% Loop over features
    for nFeature = 1:numFeatures
        
        if strcmpi(featureInfo(nFeature).modality, 'fmri') || strcmpi(featureInfo(nFeature).modality, 'smri')
            tempV = ica_fuse_spm_vol(featureInfo(nFeature).files);
            % Loop over files
            for nFiles = 1:length(tempV)
                temp = ica_fuse_read_vols(tempV(nFiles));
                temp(isnan(temp)) = 0;
                temp = (temp ~= 0);
                if (nFiles == 1)
                    mask_ind = temp;
                else
                    mask_ind = mask_ind & temp;
                end
            end
            % End loop over files
            
            maskIndices(nFeature).ind = mask_ind;
            clear mask_ind;
            clear currentData;
        elseif strcmpi(featureInfo(nFeature).modality, 'eeg')
            if (isempty(EEG_DATA_INDICES))
                tmp = ica_fuse_loadData(deblank(featureInfo(nFeature).files(1, :)));
                maskIndices(nFeature).ind = (1:size(tmp, 1));
                clear tmp;
            else
                maskIndices(nFeature).ind = EEG_DATA_INDICES;
            end
        elseif (strcmpi(featureInfo(nFeature).modality, 'gene') || strcmpi(featureInfo(nFeature).modality, 'behavioral'))
            [currentData, size_data] = getCurrentData(featureInfo(nFeature));
            clear currentData;
            maskIndices(nFeature).ind = [1:size_data(1)];
        else
            error(['Unknown modality: ', featureInfo(nFeature).modality, ' entered']);
        end
        
    end
    %% End loop over features
    
    matchInd = strmatch('fmri', lower(str2mat(featureInfo.modality)), 'exact');
    
    % Check fmri features
    if ~isempty(matchInd)
        count = 0;
        % collect all fmri features that have same dimensions
        for nn = 1:length(matchInd)
            ind = maskIndices(matchInd(nn)).ind;
            if nn == 1
                [dims] = size(ind);
            end
            
            checkInd = find(dims == size(ind));
            
            if length(checkInd) == length(dims)
                count = count + 1;
                featureVec(count) = matchInd(nn);
            end
        end
        % end for collecting fmri features having same dimensions
        
        if exist('featureVec', 'var')
            
            clear ind
            
            % Do a boolean and of all fmri features
            for nn = 1:length(featureVec)
                
                if nn == 1
                    ind = maskIndices(featureVec(nn)).ind;
                else
                    ind = ind & maskIndices(featureVec(nn)).ind;
                end
                
            end
            % End for doing boolean and of all fmri features
            
            % Assign ind vector
            for nn = 1:length(featureVec)
                maskIndices(featureVec(nn)).ind = ind;
            end
            % End for assiging ind vector
            
        end
        
    end
    % End for checking fmri features
    
else
    for nFeature = 1:length(featureInfo)
        [maskIndices(nFeature).ind] = getMaskData(featureInfo(nFeature), mask_files{nFeature});
    end
    
end
% End for default mask


function [currentData, size_data] = getCurrentData(featureInfo)
%% Load feature data

% Load data
currentData = ica_fuse_loadData(featureInfo.files);

size_data = size(currentData);

if length(size_data) == 3
    size_data(4) = 1;
end


function [maskData] = getMaskData(featureInfo, mask_file)
%% Mask data

modality = featureInfo.modality;

if strcmpi(modality, 'fmri') | strcmpi(modality, 'smri')
    
    dispStr = str2mat(['Using mask from the file ', mask_file], ...
        ['to find the regions that are not zero in the mask for feature ', featureInfo.feature_name]);
    disp(dispStr);
    fprintf('\n');
    files = featureInfo.files;
    fileName = deblank(files(1, :));
    V = ica_fuse_getVol(fileName); V = V(1);
    
    % Calculate non-zero indices for selected mask
    dims = V(1).dim(1:3); % Data dimensions
    [data] = ica_fuse_loadData(deblank(mask_file)); % Mask data
    size_data = size(data);
    
    % Check the dimensions of the mask w.r.t data
    if length(find(size_data == dims)) ~= length(dims)
        error('Error:MaskDim', ['Mask dimensions ([%s]) doesn''t match that of data dimensions ([%s])'], ...
            num2str(size_data), num2str(dims));
    end
    
    maskData = (data ~= 0);
    
    if isempty(find(maskData ~= 0))
        error(['Please check the mask (', deblank(mask_file), ')you have specified for feature ', ...
            featureInfo.feature_name]);
    end
    
elseif strcmpi(modality, 'eeg')
    
    dispStr = ['Using EEG mask for feature ', featureInfo.feature_name];
    disp(dispStr);
    if isempty(mask_file)
        error(['EEG mask for feature ', featureInfo.feature_name, ' is empty']);
    end
    fprintf('\n');
    maskData = str2num(mask_file);
    
elseif (strcmpi(modality, 'gene') || strcmpi(modality, 'behavioral'))
    
    dispStr = ['Using Gene/behavioral mask for feature ', featureInfo.feature_name];
    disp(dispStr);
    if isempty(mask_file)
        error(['Gene/behavioral mask for feature ', featureInfo.feature_name, ' is empty']);
    end
    fprintf('\n');
    maskData = str2num(mask_file);
    
end

function stackInfo = ica_fuse_prepare_data(varargin)
%% Prepare data for joint ICA analysis
%
% Inputs must be in pairs
%
% Outputs:
% stackInfo - data is stacked by features
%

computeMean = 0;
estim_data = 0;

% Loop over number of arguments
for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'dataInfo')
        dataInfo = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'mask_ind')
        mask_ind = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'normalize_scheme')
        normalizeVal = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'voxels')
        voxels = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'sel_groups')
        selGroups = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'sel_features')
        selFeatures = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'num_subjects')
        numSubjects = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'feature_norm_para')
        featureNormPara = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'compute_mean')
        computeMean = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'estim_data')
        estim_data = varargin{ii + 1};
    end
end
% end loop over number of arguments

if ~exist('dataInfo', 'var')
    error('dataInfo variable is not passed');
end

%% Number of groups and features
numGroups = length(dataInfo);
numFeatures = length(dataInfo(1).feature);

if ~exist('selGroups', 'var')
    selGroups = (1:numGroups);
end

if ~exist('selFeatures', 'var')
    selFeatures = (1:numFeatures);
end

%% Truncate Feature info
featureInfo = ica_fuse_get_feature_info(dataInfo, selGroups, selFeatures);

%% Truncate mask indices
mask_ind = mask_ind(selFeatures);

%% Selected modalities
selectedModalities = cellstr(str2mat(featureInfo.modality));
selectedFeatures = cellstr(str2mat(featureInfo.feature_name));

combName = ica_fuse_formatStr(selectedFeatures, ' & ');

disp(['Loading data for ', combName, ' ...']);


if (~estim_data)
    
    % Apply mask
    dataN = ica_fuse_applyMask(featureInfo, mask_ind);
    
else
    
    %dataN = repmat(struct('data', [], 'xAxis', [], 'files', [], 'feature_name', []), 1, length(featureInfo));
    
    estimationInfo = repmat(struct('comp', [], 'mdl', [], 'aic', [], 'feature_name', ''), 1, length(featureInfo));
    for nF = 1:length(featureInfo)
        
        disp(['Doing dimensionality estimation for feature ', featureInfo(nF).feature_name, ' ...']);
        
        % Apply mask
        dataN = ica_fuse_applyMask(featureInfo(nF), mask_ind(nF));
        
        % Do iid sampling for fmri, smri data
        doSampling = (strcmpi(featureInfo(nF).modality, 'fmri') || strcmpi(featureInfo(nF).modality, 'smri'));
        
        estimationInfo(nF).feature_name = featureInfo(nF).feature_name;
        
        [estimationInfo(nF).comp, estimationInfo(nF).mdl, estimationInfo(nF).aic] = ica_fuse_estim_dim(dataN.data, 'maskvec', mask_ind(nF).ind, 'doSampling', doSampling);
        
        clear dataN;
        
        fprintf('\n');
        
        disp(['Components estimated to be ', num2str(estimationInfo(nF).comp)]);
        
        fprintf('\n');
        
        
    end
    
    if (length(featureInfo) > 1)
        
        disp('Using PCA-CCA based estimation ..');
        allCombinations = nchoosek(1:numFeatures, 2);
        allFeatureNames = cellstr(char(featureInfo.feature_name));
        ccaEstimationInfo = repmat(struct('name', [], 'comp', []), 1, size(allCombinations, 1));
        for nF = 1:size(allCombinations, 1)
            featureA = allCombinations(nF, 1);
            featureB = allCombinations(nF, 2);
            % feature 1 and 2
            xData = ica_fuse_applyMask(featureInfo(featureA), mask_ind(featureA));
            yData = ica_fuse_applyMask(featureInfo(featureB), mask_ind(featureB));
            ccaEstimationInfo(nF).name = [allFeatureNames{featureA}, ' & ', allFeatureNames{featureB}];
            disp(['Loading features ', ccaEstimationInfo(nF).name, '...']);
            ccaEstimationInfo(nF).comp = ica_fuse_model_est_cca(xData.data', yData.data');
            
        end
        
        stackInfo.ccaEstimationInfo = ccaEstimationInfo;
        
    end
    
    stackInfo.estimationInfo = estimationInfo;
    
    return;
    
end

disp('Done');

% Get the flag for standardizing subjects
standardize_sub = ica_fuse_flagStandardizeSub(selectedModalities);

checkEEG = strmatch('eeg', selectedModalities, 'exact');

if ~isempty(checkEEG)
    
    checkEEG = checkEEG(:)';
    
    % Interpolate the EEG Data
    for nn = checkEEG
        % EEG Data length
        eegDataLength = size(dataN(nn).data, 2);
        if (eegDataLength < voxels)
            
            disp(['Resampling data of feature ', featureInfo(nn).feature_name]);
            % Loop over number of subjects
            for ii = 1:size(dataN(nn).data, 1)
                % Resample data
                temp = ica_fuse_resample(dataN(nn).data(ii, :), voxels, eegDataLength);
                if (ii == 1)
                    new_data = zeros(size(dataN(nn).data, 1), length(temp));
                end
                new_data(ii, :) = temp;
                clear temp;
            end
            % end loop over number of subjects
            
            % Resample x axis
            xAxis = ica_fuse_resample(dataN(nn).xAxis, voxels, eegDataLength);
            
            dataN(nn).data = new_data;
            dataN(nn).xAxis = xAxis;
            clear xAxis; clear new_data;
            
        end
        
    end
    % End for interpolating data
    
end

% Standard deviation parameters
stdParameters = repmat(struct('stdValues', []), 1, length(dataN));
dataLength = zeros(1, length(dataN));

for nD = 1:length(dataLength)
    dataLength(nD) = size(dataN(nD).data, 2);
    stdParameters(nD).stdValues = ones(1, size(dataN(nD).data, 1));
end

%% Compute mean before applying normalization parameters
meanData = [];
meanDataGroups = [];
if (computeMean)
    %% Truncate subjects
    numSubjects = numSubjects(selGroups);
    [meanData, meanDataGroups] = ica_fuse_compute_mean_features(dataN, dataLength, numSubjects);
end


if ~exist('featureNormPara', 'var')
    
    % Do normalization
    [dataN, featureNormPara] = ica_fuse_featureNormalize(dataN, normalizeVal);
    
else
    
    %% Truncate normalization parameters
    featureNormPara = featureNormPara(selFeatures);
    
    %%%% Do normalization %%%%%%
    % Loop over features
    for nF = 1:length(dataN)
        currentData = dataN(nF).data;
        currentData = currentData / featureNormPara(nF);
        dataN(nF).data = currentData;
        clear currentData;
    end
    % End loop over features
    %%% End for feature normalization %%%
    
end

if strcmpi(standardize_sub, 'yes')
    disp('Converting each subject''s data to z-scores ...');
    % Normalize data for each subject by its standard deviation (Memory
    % constrained way)
    for nFeatures = 1:length(dataN)
        stdValues = ones(1, size(dataN(nFeatures).data, 1));
        % Convert subject data to z scores
        for nn = 1:length(stdValues)
            stdValues(nn) = std(dataN(nFeatures).data(nn, :));
            dataN(nFeatures).data(nn, :) = dataN(nFeatures).data(nn, :) ./ stdValues(nn);
        end
        stdParameters(nFeatures).stdValues = stdValues;
        clear stdValues;
    end
    % end loop over features
    disp('Done converting each subject''s data to z-scores');
end

fprintf('\n');
disp('Stacking data ...');
fprintf('\n');

timeAxis = repmat(struct('data', []), 1, length(dataN));
fileInfo = repmat(struct('name', [], 'feature_name', []), 1, length(dataN));
for nF = 1:length(timeAxis)
    timeAxis(nF).data = dataN(nF).xAxis;
    fileInfo(nF).name = deblank(featureInfo(nF).files(1, :));
    fileInfo(nF).feature_name = deblank(featureInfo(nF).feature_name);
end

%% Stack information
stackInfo.combName = combName;
stackInfo.featureNormPara = featureNormPara;
stackInfo.stdParameters = stdParameters;
stackInfo.timeAxis = timeAxis;
clear timeAxis;
stackInfo.data = [dataN.data];
clear dataN;
stackInfo.dataLength = dataLength;
stackInfo.meanData = meanData;
stackInfo.meanDataGroups = meanDataGroups;
stackInfo.fileInfo = fileInfo;

disp('Done stacking data');
fprintf('\n');


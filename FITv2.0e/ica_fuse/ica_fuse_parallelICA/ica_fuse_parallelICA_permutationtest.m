function varargout = ica_fuse_parallelICA_permutationtest(parallel_ica_fusion_file, num_permutations, numICARuns)
% Run parallel ICA fusion


% Global vars
ica_fuse_defaults;
global PARALLEL_ICA_INFO_MAT_FILE;

% Open file selection dialog box
if ~exist('parallel_ica_fusion_file', 'var')
    parallel_ica_fusion_file = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        ['*', PARALLEL_ICA_INFO_MAT_FILE, '*.mat'], 'title', 'Select parallel ICA information file');
end

if isempty(parallel_ica_fusion_file)
    error(['Parallel ICA fusion information file is not selected']);
end

outputDir = fileparts(parallel_ica_fusion_file);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

load(parallel_ica_fusion_file);

try
    maxCorr = max(abs(paraICAInfo.run_analysis.average_correlation));
catch
    error('You need to run parallel ICA in order to do permutation test');
end


if (~exist('numICARuns', 'var'))
    numICARuns = 15;
end

if (~exist('num_permutations', 'var'))
    
    answer = ica_fuse_inputdlg2({'Enter no. of permutations', 'Enter no. of ICA runs'}, 'No. of permutations & ICA runs', 1, {'1000', num2str(numICARuns)});
    
    drawnow;
    
    if (isempty(answer) || isempty(answer{1}))
        error('No. of permutations is not entered');
    end
    
    num_permutations = str2num(answer{1});
    
    if isempty(answer{1})
        error('No. of ICA runs is not entered');
    end
    
    numICARuns = str2num(answer{2});
    
end

tic;

fprintf('\n');

disp('Doing permutation test');

fprintf('\n');

% Get information from paraICAInfo file
numSubjects = paraICAInfo.setup_analysis.numSubjects;
dataInfo = paraICAInfo.setup_analysis.dataInfo;
modalities = cellstr(char(dataInfo(1).feature.modality));

try
    numComp = paraICAInfo.setup_analysis.numComp;
catch
    modality1_numComp = paraICAInfo.setup_analysis.modality1_numComp;
    modality2_numComp = paraICAInfo.setup_analysis.modality2_numComp;
    numComp = [modality1_numComp, modality2_numComp];
end

% modality1_numComp = paraICAInfo.setup_analysis.modality1_numComp;
% modality2_numComp = paraICAInfo.setup_analysis.modality2_numComp;
% numComp = [modality1_numComp, modality2_numComp];


try
    type_pca = paraICAInfo.setup_analysis.type_pca;
catch
    type_pca = 'standard';
end


try
    type_parallel_ica = paraICAInfo.setup_analysis.type_parallel_ica;
catch
    type_parallel_ica = 'aa';
end

%numICARuns = paraICAInfo.setup_analysis.numICARuns;

% Create mask
%mask_ind = ica_fuse_createMask(dataInfo, maskFile);
mask_ind = paraICAInfo.run_analysis.mask_ind;

% Get featureInfo
featureInfo = ica_fuse_get_feature_info(dataInfo);

% Apply mask
featureData = ica_fuse_applyMask(featureInfo, mask_ind);

% Initialise reference vector
reference = ones(sum(numSubjects), 1);

% Print to command window
if strcmpi(paraICAInfo.setup_analysis.type_pca, 'standard')
    disp('Using standard PCA ...');
else
    disp('Using PCA with reference ...');
    
    if isfield(paraICAInfo.setup_analysis, 'reference') && ~isempty(paraICAInfo.setup_analysis.reference)
        reference = paraICAInfo.setup_analysis.reference;
    else
        
        if (length(numSubjects) == 2)
            reference(1:numSubjects(1)) = -1;
        else
            for nS = 1:length(numSubjects)
                groupInd = ica_fuse_get_groupInd(nS, numSubjects);
                reference(groupInd) = 1/numSubjects(nS);
            end
        end
    end
end

corrVals = zeros(1, num_permutations);
for nPerm = 1:num_permutations
    
    fprintf('\n');
    disp(['Permutation ', num2str(nPerm), '/', num2str(num_permutations)]);
    
    fprintf('\n');
    
    
    % Data reduction
    pca_output = parallelICA_dataReduction(featureData, modalities, numComp, reference, type_pca, ...
        type_parallel_ica);
    
    data = {pca_output(1).whitesig, pca_output(2).whitesig}; % Data
    dewhiteM = {pca_output(1).dewhiteM, pca_output(2).dewhiteM}; % Dewhitening matrix
    whiteM = {pca_output(1).whiteM, pca_output(2).whiteM}; % Whitening matrix
    
    clear pca_output;
    
    [aveComp, loadingCoeff, avecorr] = ica_fuse_calculate_parallelICA(data, dewhiteM, whiteM, ...
        numICARuns, type_parallel_ica, modalities, paraICAInfo.setup_analysis.ICA_Options, 'icassso', featureData);
    
    [dd, dd_inds] = max(abs(avecorr));
    
    corrVals(nPerm) = avecorr(dd_inds);
    
    fprintf('\n');
    
end

fprintf('\n');

p = (length(find(abs(corrVals) >= maxCorr)) / num_permutations);

disp(['Ratio of correlation values above the maximally linked component is ', num2str(p)]);
fprintf('\n');

t_end = toc;

disp(['Time taken to run the permutation analysis is ', num2str(t_end), ' seconds']);
fprintf('\n');

if (nargout > 0)
    
    varargout{1} = corrVals;
    varargout{2} = p;
    
else
    
    %% Save info
    paraICAInfo.run_analysis.permutation_test_corr_vals = corrVals;
    ica_fuse_save(parallel_ica_fusion_file, 'paraICAInfo');
    
end

fprintf('\n');


function [maxcorr, maxrow, maxcol] = computeCorr(loadingCoeff, aveComp, type_parallel_ica)
%% Compute max correlation
%

maxcorr = zeros(1, size(aveComp{1}, 1));
maxcol = []; maxrow = [];

% Match components
for j=1:size(loadingCoeff{1}, 2)
    for m = 1:size(aveComp{2}, 1)
        a = loadingCoeff{1}(:, j);
        if strcmpi(type_parallel_ica, 'as')
            b = aveComp{2}(m, :)';
        else
            b = loadingCoeff{2}(:, m);
        end
        temp = ica_fuse_corr(a, b); % correlation
        if abs(temp) > abs(maxcorr(j))
            maxcol(j)=j;
            maxrow(j)=m;
            maxcorr(j)=temp;
        end
    end
end


function pca_output = parallelICA_dataReduction(featureData, modalities, numComp, reference, pca_type, type_parallel_ica)
% Data reduction step

pca_output = repmat(struct('dewhiteM', [], 'whiteM', [], 'whitesig', []), 1, length(modalities));


nSubs = size(featureData(1).data, 1);

% Remove the mean only for fmri modality
for nModality = 1:length(modalities)
    
    ia = randperm(nSubs);
    
    fprintf('\n');
    currentData = featureData(nModality).data(ia, :);
    if strcmpi(modalities{nModality}, 'fmri')
        currentData = currentData';
        disp(['Removing mean for feature ', featureData(nModality).feature_name, '...']);
        % Remove mean of fmri data for each subject
        try
            tempVar = ones(size(currentData, 1), 1)*mean(currentData);
            currentData = currentData - tempVar;
            clear tempVar;
        catch
            clear tempVar;
            % --slower way to zero mean data
            disp('Using slower but less memory constrained way to zero mean data');
            
            for j = 1:size(currentData, 2)
                currentData(:, j) = currentData(:,  j) - mean(currentData(:, j));
            end
        end
        currentData = currentData';
    elseif strcmpi(modalities{nModality}, 'gene')
        if strcmpi(type_parallel_ica, 'as')
            % Arrange SNP data as SNPS by subjects
            currentData = currentData';
        end
    end
    
    featureData(nModality).data = currentData;
    clear currentData;
    
    disp(['Doing data reduction for feature ', featureData(nModality).feature_name]);
    %%%%%%%%%%%%%%% PCA Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(pca_type, 'standard')
        
        % Calculate PCA
        [V, Lambda] = ica_fuse_v_pca(featureData(nModality).data, 1, numComp(nModality), 0, 'transpose', 'no');
        
    else
        % Calculate PCA with reference
        [V, Lambda] = ica_fuse_pca_reference(featureData(nModality).data, 1, numComp(nModality), reference(ia), 0, 'transpose', 'no');
    end
    
    % Whitening
    [whitesig, whiteM, dewhiteM] = ica_fuse_v_whiten(featureData(nModality).data, V, Lambda, 'transpose');
    disp('Done');
    
    
    fprintf('\n');
    
    % Store PCA output
    pca_output(nModality).dewhiteM = dewhiteM;
    pca_output(nModality).whiteM = whiteM;
    pca_output(nModality).whitesig = whitesig;
    
end
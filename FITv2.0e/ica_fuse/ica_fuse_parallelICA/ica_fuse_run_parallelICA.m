function ica_fuse_run_parallelICA(parallel_ica_fusion_file)
% Run parallel ICA fusion


% Global vars
ica_fuse_defaults;
global PARALLEL_ICA_INFO_MAT_FILE;
global PARALLEL_ICA_COMPONENT_NAMING;
global MRI_DATA_FILTER;
global FLIP_ANALYZE_IM;
global PARALLEL_ICA_PCA_FILE;
global PARALLEL_ICA_ICA_FILE;
global OPEN_DISPLAY_WINDOW;

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

output_prefix = paraICAInfo.setup_analysis.prefix;

% Results log file
output_LogFile = fullfile(outputDir, [output_prefix, '_para_ica_results.log']);

diary(output_LogFile);

tic;

disp('Starting parallel ICA fusion analysis');

try
    
    drawnow;
    
    % Get information from paraICAInfo file
    maskFile = paraICAInfo.setup_analysis.maskFile;
    numGroups = paraICAInfo.setup_analysis.numGroups;
    numFeatures = paraICAInfo.setup_analysis.numFeatures;
    numSubjects = paraICAInfo.setup_analysis.numSubjects;
    dataInfo = paraICAInfo.setup_analysis.dataInfo;
    modalities = cellstr(str2mat(dataInfo(1).feature.modality));
    
    preproc_type = 'none';
    try
        preproc_type = paraICAInfo.setup_analysis.preproc_type;
    catch
    end
    
    try
        numComp = paraICAInfo.setup_analysis.numComp;
    catch
        modality1_numComp = paraICAInfo.setup_analysis.modality1_numComp;
        modality2_numComp = paraICAInfo.setup_analysis.modality2_numComp;
        numComp = [modality1_numComp, modality2_numComp];
    end
    
    numICARuns = paraICAInfo.setup_analysis.numICARuns;
    
    % Store dataInfo, numGroups, numFeatures and numSubjects
    paraICAInfo.run_analysis.prefix = output_prefix;
    paraICAInfo.run_analysis.maskFile = maskFile;
    paraICAInfo.run_analysis.numGroups = numGroups;
    paraICAInfo.run_analysis.numFeatures = numFeatures;
    paraICAInfo.run_analysis.numSubjects = numSubjects;
    paraICAInfo.run_analysis.modalities = modalities;
    paraICAInfo.run_analysis.numComp = numComp;
    if ~isfield(paraICAInfo.setup_analysis, 'ICA_Options')
        paraICAInfo.setup_analysis.ICA_Options = ica_fuse_paraICAOptions(min(numComp), 'off');
    end
    
    paraICAInfo.run_analysis.ICA_Options = paraICAInfo.setup_analysis.ICA_Options;
    
    % data Info
    paraICAInfo.run_analysis.dataInfo = dataInfo;
    
    % Store flip parameter
    paraICAInfo.run_analysis.flip_analyze_images = FLIP_ANALYZE_IM;
    
    
    % When using ica_fuse_runica_parallelicaMul_AA function arrange SNP data
    % subjects by SNPS
    if ~isfield(paraICAInfo.setup_analysis, 'type_parallel_ica')
        type_parallel_ica = 'as';
    else
        type_parallel_ica = paraICAInfo.setup_analysis.type_parallel_ica;
    end
    
    
    if (paraICAInfo.run_analysis.numFeatures == 3)
        type_parallel_ica = 'aa';
    end
    
    
    paraICAInfo.run_analysis.type_parallel_ica = type_parallel_ica;
    
    if ~isfield(paraICAInfo.setup_analysis, 'type_pca')
        paraICAInfo.setup_analysis.type_pca = 'standard';
    end
    
    if strcmpi(type_parallel_ica, 'as')
        %% Arrange SNP data as SNPS by subjects
        %SNPData = SNPData';
        paraICAInfo.setup_analysis.type_pca = 'standard';
    end
    
    paraICAInfo.run_analysis.type_pca = paraICAInfo.setup_analysis.type_pca;
    
    
    type_ica = 'average';
    try
        type_ica = lower(paraICAInfo.setup_analysis.type_ica);
    catch
    end
    
    if strcmpi(type_ica, 'icasso') && numICARuns < 2
        error('You need to select atleast 2 runs inorder to run ICASSO');
    end
    
    
    paraICAInfo.setup_analysis.type_ica = type_ica;
    paraICAInfo.run_analysis.type_ica = type_ica;
    
    % Get featureInfo
    featureInfo = ica_fuse_get_feature_info(dataInfo);
    
    % Create mask
    mask_ind = ica_fuse_createMask(dataInfo, maskFile);
    
    paraICAInfo.run_analysis.mask_ind = mask_ind;
    
    % Apply mask
    featureData = ica_fuse_applyMask(featureInfo, mask_ind);
    
    fprintf('\n');
    disp('-----------------------------------------------------------------------------------------------');
    disp(['Calculating PCA ']);
    disp('-----------------------------------------------------------------------------------------------');
    
    fprintf('\n');
    
    % Initialise reference vector
    reference = ones(sum(numSubjects), 1);
    
    % Print to command window
    if strcmpi(paraICAInfo.setup_analysis.type_pca, 'standard')
        disp('Running standard PCA ...');
    else
        disp('Running PCA with reference ...');
        
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
    
    
    % Preprocess data
    featureData = preprocessData(featureData, preproc_type);
    
    
    % Data reduction
    pca_output = ica_fuse_parallelICA_dataReduction(featureData, modalities, numComp, reference, paraICAInfo.setup_analysis.type_pca, ...
        type_parallel_ica);
    
    pcaFile = [output_prefix, '_', PARALLEL_ICA_PCA_FILE, '.mat'];
    ica_fuse_save(fullfile(outputDir, pcaFile), 'pca_output');
    %clear featureData;
    
    fprintf('\n');
    
    disp('-----------------------------------------------------------------------------------------------');
    disp(['Done calculating PCA']);
    disp('-----------------------------------------------------------------------------------------------');
    
    paraICAInfo.run_analysis.pcaFile = pcaFile;
    
    % Run parallel ICA
    %ica_output = repmat(struct('aveComp', [], 'loadingCoeff', []), 1, length(modalities));
    
    data = cell(1, length(pca_output));
    dewhiteM = data;
    whiteM = data;
    
    for nM = 1:length(data)
        data{nM} = pca_output(nM).whitesig; % Data
        dewhiteM{nM} = pca_output(nM).dewhiteM; % Dewhitening matrix
        whiteM{nM} = pca_output(nM).whiteM; % Whitening matrix
    end
    
    clear pca_output;
    
    fprintf('\n');
    
    disp('-----------------------------------------------------------------------------------------------');
    disp(['Calculating ICA ']);
    disp('-----------------------------------------------------------------------------------------------');
    
    fprintf('\n');
    
    if (length(data) == 2)
        [aveComp, loadingCoeff, avecorr, corrIndices, sR_modality1, sR_modality2] = ica_fuse_calculate_parallelICA(data, dewhiteM, whiteM, ...
            numICARuns, type_parallel_ica, modalities, paraICAInfo.run_analysis.ICA_Options, paraICAInfo.run_analysis.type_ica, featureData);
    else
        [aveComp, loadingCoeff, avecorr, corrIndices, sR_modality1, sR_modality2, sR_modality3] = ica_fuse_calculate_parallelICA(data, dewhiteM, whiteM, ...
            numICARuns, type_parallel_ica, modalities, paraICAInfo.run_analysis.ICA_Options, paraICAInfo.run_analysis.type_ica, featureData);
    end
    clear data;
    
    if (strcmpi(paraICAInfo.run_analysis.type_ica, 'icasso'))
        icasso_file = fullfile(outputDir, [output_prefix, '_para_ica_icasso_results.mat']);
        ica_fuse_save(icasso_file, 'sR_modality1', 'sR_modality2');
        if (exist('sR_modality3', 'var'))
            ica_fuse_save(icasso_file, 'sR_modality3', '-append');
        end
    end
    
    % ICA Output file
    icaFile = [output_prefix, '_', PARALLEL_ICA_ICA_FILE, '.mat'];
    ica_fuse_save(fullfile(outputDir, icaFile), 'aveComp', 'loadingCoeff', 'avecorr', 'corrIndices');
    
    fprintf('\n');
    
    disp('-----------------------------------------------------------------------------------------------');
    disp(['Done calculating ICA ']);
    disp('-----------------------------------------------------------------------------------------------');
    
    fprintf('\n');
    
    % Scale EEG components to data units
    [aveComp, loadingCoeff] = scaleComps(aveComp, loadingCoeff, dataInfo, mask_ind);
    
    if (length(loadingCoeff) == 2)
        % Compute correlations part
        [avecorr, corrIndices] = computeCorr(loadingCoeff, aveComp, type_parallel_ica);
    else
        [bestIndex, pair1_corr, pair2_corr, pair3_corr] = ica_fuse_mostCorrelatedTrio(loadingCoeff{1}, loadingCoeff{2}, loadingCoeff{3});
        avecorr = [pair1_corr(1), pair2_corr(1), pair3_corr(1)];
        corrIndices = bestIndex(1, 1:3);
    end
    
    paraICAInfo.run_analysis.icaFile = icaFile;
    
    paraICAInfo.run_analysis.average_correlation = avecorr;
    
    paraICAInfo.run_analysis.corrIndices = corrIndices;
    
    
    % Save parallel ICA fusion information file
    ica_fuse_save(parallel_ica_fusion_file, 'paraICAInfo');
    
    imExtn = '.img';
    if strcmpi(MRI_DATA_FILTER, '*.nii')
        imExtn = MRI_DATA_FILTER(2:end);
    end
    
    fprintf('\n');
    %%%%%%%%%%%% Write Results %%%%%%%%%%%%%
    disp('-----------------------------------------------------------------------------------------------');
    disp(['Writing component results ']);
    disp('-----------------------------------------------------------------------------------------------');
    fprintf('\n');
    
    nFeature = 1;
    outFile = [output_prefix, '_', PARALLEL_ICA_COMPONENT_NAMING];
    disp(['Writing data for component set ', outFile]);
    
    % Write data as images or ascii files
    outputFiles = ica_fuse_write_parallel_ica_comp(aveComp, loadingCoeff, dataInfo, mask_ind, outFile, outputDir, numSubjects);
    
    paraICAInfo.run_analysis.outputFiles = outputFiles;
    
    %%%%%%%%% End for writing results %%%%%%%%
    fprintf('\n');
    
    disp('-----------------------------------------------------------------------------------------------');
    disp(['Done writing component results']);
    disp('-----------------------------------------------------------------------------------------------');
    
    
    fprintf('\n');
    ica_fuse_save(parallel_ica_fusion_file, 'paraICAInfo');
    
    disp(['Parallel ICA fusion information is saved in MAT file: ', parallel_ica_fusion_file]);
    
    t_end = toc;
    
    % Suspend diary
    diary('off');
    
    disp(['Time taken to run the analysis is ', num2str(t_end), ' seconds']);
    fprintf('\n');
    disp(['All the analysis information is stored in the file ', output_LogFile]);
    
    disp('Finished with Analysis');
    fprintf('\n');
    
    if OPEN_DISPLAY_WINDOW
        disp('Opening display GUI ...');
        ica_fuse_parallelICA_displayGUI(parallel_ica_fusion_file);
    end
    
catch
    diary('off');
    % Display error message
    ica_fuse_displayErrorMsg;
end

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

function [aveComp, loadingCoeff] = scaleComps(aveComp, loadingCoeff, dataInfo, mask_ind)

ica_fuse_defaults;
global FLIP_SIGN_COMPONENTS;

modalities = cellstr(char(dataInfo(1).feature.modality));

eegInds = strmatch('eeg', lower(modalities), 'exact');
eegInds = eegInds(:)';

for nM = eegInds
    
    disp(['Scaling feature ', dataInfo(1).feature(nM).name, ' components w.r.t the mean of the data ...']);
    
    comps = detrend(aveComp{nM}', 0);
    featureInfo = ica_fuse_get_feature_info(dataInfo, 1:length(dataInfo), nM);
    meanData = ica_fuse_loadData(char(featureInfo.files));
    clear featureInfo;
    meanData = mean(squeeze(meanData(:, 2, :)), 2);
    meanData = meanData(mask_ind(nM).ind);
    meanData = detrend(meanData(:), 0);
    betas = ica_fuse_regress(meanData, comps);
    
    A = loadingCoeff{nM};
    
    % Loop over components
    for nComp = 1:size(comps, 2)
        if (sign(betas(nComp)) == -1 && FLIP_SIGN_COMPONENTS)
            disp(['Flipping sign for component ', num2str(nComp), ' of feature ', dataInfo(1).feature(nM).name]);
        end
        [fileIndex] = ica_fuse_returnFileIndex(nComp);
        % Current component
        currentData = comps(:, nComp);
        if (FLIP_SIGN_COMPONENTS)
            currentData = betas(nComp)*currentData;
            A(:, nComp) = sign(betas(nComp))*A(:, nComp);
        else
            currentData = abs(betas(nComp))*currentData;
        end
        comps(:, nComp) = currentData;
    end
    
    aveComp{nM} = comps';
    loadingCoeff{nM} = A;
    
    fprintf('Done\n\n');
end



function featureData = preprocessData(featureData, preproc_type)
%% Scale data to z-scores
%

if (~strcmpi(preproc_type, 'none'))
    for nModality = 1:length(featureData)
        disp(['Converting data to z-scores for feature ', featureData(nModality).feature_name, '...']);
        currentData = (featureData(nModality).data)';
        if (~strcmpi(featureData(nModality).modality, 'gene'))
            currentData = ica_fuse_remove_mean(currentData);
        end
        currentData = currentData*diag(1./std(currentData));
        featureData(nModality).data = currentData';
    end
end
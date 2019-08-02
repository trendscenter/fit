function ica_fuse_leave_one_out(parallel_ica_fusion_file)
% Function to do leave one out evaluation

%% Load defaults
ica_fuse_defaults;

global PARALLEL_ICA_INFO_MAT_FILE;
global UI_FG_COLOR;
% FONT DEFAULTS
global UI_FONT_NAME; % font name
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size

% Open file selection dialog box
if ~exist('parallel_ica_fusion_file', 'var')
    parallel_ica_fusion_file = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        ['*', PARALLEL_ICA_INFO_MAT_FILE, '*.mat'], 'title', 'Select parallel ICA information file');
end

if isempty(parallel_ica_fusion_file)
    error('Parallel ICA fusion information file is not selected');
end

outputDir = fileparts(parallel_ica_fusion_file);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

drawnow;

fprintf('\n');
disp('Using leave one out evaluation to test the stability of components.');
disp('This may take lot of time depending upon the number of subjects and number of ICA runs.');
fprintf('\n');

load(parallel_ica_fusion_file);

if ~isfield(paraICAInfo, 'run_analysis')
    error('Please run the parallel ICA analysis in order to do leave one out evaluation');
end

output_prefix = paraICAInfo.run_analysis.prefix;

%% Get information from paraICAInfo file
[d, fN, extn] = fileparts(paraICAInfo.run_analysis.icaFile);
icaFile = fullfile(outputDir, [fN, extn]);
load(icaFile);
[sorted_values, inds] = sort(abs(avecorr));
inds = inds(end:-1:1);
bestComp = {aveComp{1}(inds(1), :), aveComp{2}(corrIndices(inds(1)), :)};
maxcorr = sorted_values(end);
clear aveComp avecorr;

numSubjects = paraICAInfo.run_analysis.numSubjects;
dataInfo = paraICAInfo.run_analysis.dataInfo;
modalities = cellstr(str2mat(dataInfo(1).feature.modality));
numComp = paraICAInfo.run_analysis.numComp;
featureNames = cellstr(str2mat(dataInfo(1).feature.name));

numICARuns = paraICAInfo.setup_analysis.numICARuns;

if ~isfield(paraICAInfo.run_analysis, 'ICA_Options')
    ICA_Options = ica_fuse_paraICAOptions(min(numComp), 'off');
else
    ICA_Options = paraICAInfo.run_analysis.ICA_Options;
end

%% When using ica_fuse_runica_parallelicaMul_AA function arrange SNP data subjects by SNPS
if ~isfield(paraICAInfo.run_analysis, 'type_parallel_ica')
    type_parallel_ica = 'as';
else
    type_parallel_ica = paraICAInfo.run_analysis.type_parallel_ica;
end

if ~isfield(paraICAInfo.run_analysis, 'type_pca')
    type_pca = 'standard';
end

if strcmpi(type_parallel_ica, 'as')
    type_pca = 'standard';
end

if isfield(paraICAInfo.run_analysis, 'type_pca')
    type_pca = paraICAInfo.run_analysis.type_pca;
end


type_ica = 'average';
try
    type_ica = lower(paraICAInfo.run_analysis.type_ica);
catch
end


% Mask indices
mask_ind = paraICAInfo.run_analysis.mask_ind;

% Get featureInfo
featureInfo = ica_fuse_get_feature_info(dataInfo);

% Apply mask
featureData = ica_fuse_applyMask(featureInfo, mask_ind);

%[featureData] = ica_fuse_getFeatureData(dataInfo);

% Apply mask
%[featureData] = ica_fuse_applyMask(featureData, mask_ind);

%% Initialise reference vector
reference = ones(sum(numSubjects), 1);

% Print to command window
if strcmpi(paraICAInfo.run_analysis.type_pca, 'standard')
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

origData = featureData;
origReference = reference;
numOfSub = length(reference);

% Initialise average components
modality1_aveComp = zeros(length(bestComp{1}), numOfSub);
modality2_aveComp = zeros(length(bestComp{2}), numOfSub);

corr_between_modalities = zeros(1, numOfSub);

%% Loop over subjects
for nSub = 1:numOfSub
    
    sub_included = find((1:numOfSub) ~= nSub);
    
    for nM = 1:length(origData)
        featureData(nM).data = origData(nM).data(sub_included, :);
        featureData(nM).files = origData(nM).files(sub_included, :);
    end
    
    reference = origReference(sub_included);
    
    % PCA output
    pca_output = ica_fuse_parallelICA_dataReduction(featureData, modalities, numComp, reference, type_pca, type_parallel_ica);
    
    % Data
    data = {pca_output(1).whitesig, pca_output(2).whitesig};
    % Dewhitening matrix
    dewhiteM = {pca_output(1).dewhiteM, pca_output(2).dewhiteM};
    % Whitening matrix
    whiteM = {pca_output(1).whiteM, pca_output(2).whiteM};
    
    % Parallel ICA
    [aveComponents, loadingCoeff] = ica_fuse_calculate_parallelICA(data, dewhiteM, whiteM, numICARuns, type_parallel_ica, ...
        modalities, ICA_Options, type_ica, featureData);
    
    % Modality 1 and 2 components
    [modality1_comp, best_comp_corr1] = best_corr_comp(aveComponents{1}, bestComp{1});
    [modality2_comp, best_comp_corr2] = best_corr_comp(aveComponents{2}, bestComp{2});
    
    % Store the average components
    modality1_aveComp(:, nSub) = sign(best_comp_corr1)*(aveComponents{1}(modality1_comp, :))';
    modality2_aveComp(:, nSub) = sign(best_comp_corr2)*(aveComponents{2}(modality2_comp, :))';
    
    % Compute correlation between linked between modality 1 and modality 2 %
    a = loadingCoeff{1}(:, modality1_comp);
    
    if strcmpi(type_parallel_ica, 'as')
        b = aveComponents{2}(modality2_comp, :)';
    else
        b = loadingCoeff{2}(:, modality2_comp);
    end
    % End for computing correlation between linked fMRI and SNP component %
    %
    
    a = ica_fuse_corr(a(:), b(:));
    corr_between_modalities(nSub) = a;
    
    clear a b;
    
end
%% End loop over subjects

fprintf('\n');

%% Correlation matrix of modalities
corr_modality = zeros(numOfSub, numOfSub, 2);

disp('Correlation less than 0.6 is set to zero in correlation matrix (corr_modality)');

% Any correlation value less than 0.6 is not taken into account
temp1 = ica_fuse_corr(modality1_aveComp);
temp1(abs(temp1) < 0.6) = 0;
corr_modality(:, :, 1) = temp1;

temp2 = ica_fuse_corr(modality2_aveComp);
temp2(abs(temp2) < 0.6) = 0;
corr_modality(:, :, 2) = temp2;


fprintf('\n');

%% Save results
resultsFile = fullfile(outputDir, [output_prefix, '_leave_one_out_results.mat']);

ica_fuse_save(resultsFile, 'corr_modality', 'corr_between_modalities', 'modality1_aveComp', 'modality2_aveComp');

clear modality1_aveComp modality2_aveComp;

disp(['Leave one out evaluation results are saved in file ', resultsFile]);

fprintf('\n');

disp('The following are the variables in the results file:');
disp(['1. corr_modality - Correlation matrix of size ', num2str(numOfSub), ' X ', num2str(numOfSub), ' X 2 where third dimension refers to modalities']);
disp(['2. corr_between_modalities - Correlation between modalities of size 1 X ', num2str(numOfSub)]);
disp(['3. modality1_aveComp - Average component of modality 1 of size ', num2str(length(bestComp{1})), ' X ', num2str(numOfSub)]);
disp(['4. modality2_aveComp - Average component of modality 2 of size ', num2str(length(bestComp{2})), ' X ', num2str(numOfSub)]);

fprintf('\n');

%% Loop over modalities
for nM = 1:length(featureNames)
    drawnow;
    %% Plot correlation matrix of modality
    plotTitle = ['Correlation Matrix Of ', featureNames{nM}];
    graphicsHandle = ica_fuse_getGraphics(plotTitle, 'normal', ['corr_matrix_modality', num2str(nM)], 'on');
    axesPos = [0.125 0.125 0.75 0.75];
    axesH = axes('parent', graphicsHandle, 'units', 'normalized', 'position', axesPos, 'fontunits', UI_FONT_UNITS, 'fontname', UI_FONT_NAME, 'fontSize', UI_FONT_SIZE);
    image(abs(squeeze(corr_modality(:, :, nM))), 'parent', axesH, 'CDataMapping', 'scaled');
    title(plotTitle, 'color',  UI_FG_COLOR, 'HorizontalAlignment', 'center', 'fontunits', UI_FONT_UNITS, 'fontname', UI_FONT_NAME, 'fontSize', UI_FONT_SIZE);
    xlabel('Subjects', 'parent', axesH);
    ylabel('Subjects', 'parent', axesH);
    set(axesH, 'YColor', UI_FG_COLOR, 'XColor', UI_FG_COLOR);
    axis(axesH, 'tight');
    
    %% Colorbar handle
    matlab_version = ica_fuse_get_matlab_version;
    if (matlab_version < 14)
        colorbarHandle = colorbar('peer', axesH);
        axesPos = get(axesH, 'position');
        %% Colorbar position
        colorbarPosition = [axesPos(1) + axesPos(3) + 0.05, axesPos(2), 0.03, axesPos(4)];
        set(colorbarHandle, 'position', colorbarPosition);
    else
        colorbarHandle = colorbar('peer', axesH, 'location', 'eastoutside');
    end
    
    set(colorbarHandle, 'YColor', UI_FG_COLOR);
    set(colorbarHandle, 'XColor', [0, 0, 0]);
    set(colorbarHandle, 'XTick', []);
    
end
%% End loop over modalities

function [best_comp_ind, best_comp_corr] = best_corr_comp(aveComponents, bestComp)
% Best correlated component which is obtained by correlating the maximal linked component
% between the modalities with the components

compCorr = zeros(1, size(aveComponents, 1));
for nC = 1:length(compCorr)
    compCorr(nC) = ica_fuse_correlateFunctions(aveComponents(nC, :), bestComp);
end

[sorted_comp, inds] = sort(abs(compCorr));
best_comp_ind = inds(end);
best_comp_corr = compCorr(best_comp_ind);
function ica_fuse_run_pgica_ica(param_file)
%% Run pgica ica
%

ica_fuse_defaults;
global PGICA_ICA_DEFAULTS;


if (isempty(which('gift.m')))
    error('Requires GIFT toolbox on path to run the back-reconstructed output from pgica');
end

if (~exist('param_file', 'var'))
    param_file = ica_fuse_selectEntry('typeEntity', 'file', 'title', ...
        'Select parallel group ICA + ICA fusion param file', 'filter', '*pgica*ica*mat');
    drawnow;
end

if (ischar(param_file))
    load(param_file);
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
else
    pgicaInfo = param_file;
    outputDir = pgicaInfo.outputDir;
end

clear param_file;

cd(outputDir);

pgicaInfo.outputDir = outputDir;
%param_file = pgicaInfo.param_file;

if (~isfield(pgicaInfo, 'featuresInfo'))
    featuresInfo = repmat(struct('name', '', 'modality_type', '', 'components', '', 'files', '', 'mask', ''), 1, num_features);
    
    for nF = 1:num_features
        featuresInfo(nF).name = ['Feature ', num2str(nF)];
        featuresInfo(nF).modality_type = pgicaInfo.modalities{nF};
        featuresInfo(nF).components = pgicaInfo.(['numComp', num2str(nF)]);
        featuresInfo(nF).mask = pgicaInfo.(['maskFile_modality', num2str(nF)]);
        featuresInfo(nF).files = pgicaInfo.files{nF};
    end
    pgicaInfo.featuresInfo = featuresInfo;
end

if (~isfield(pgicaInfo, 'ica_options'))
    ica_options = ica_fuse_paraICAOptions(min([pgicaInfo.featuresInfo.components]), 'off');
    pgicaInfo.ica_options = ica_options;
end

ica_options = pgicaInfo.ica_options;
featuresInfo = pgicaInfo.featuresInfo;

results_file = fullfile(outputDir, [pgicaInfo.prefix, '_pgica_ica.log']);
diary(results_file);

%% Step 1: Do pca and get the whitening dewhitening information
% mask1 = pgicaInfo.maskFile_modality1;
% mask2 = pgicaInfo.maskFile_modality2;

modalities = {'smri', 'fmri'};
try
    modalities = pgicaInfo.modalities;
catch
end

VStruct = ica_fuse_getVol(which('nsingle_subj_T1_2_2_5.nii'));
VStruct = VStruct(1);

structV = cell(1, length(modalities));
maskIndices = cell(1, length(modalities));
for nM = 1:length(modalities)
    mask = featuresInfo(nM).mask;
    currentFiles = featuresInfo(nM).files;
    
    if (isempty(mask))
        % Default mask
        if (strcmpi(modalities{nM}, 'smri') || strcmpi(modalities{nM}, 'dti'))
            mask = ica_fuse_generateMask(currentFiles, 'multiplier', PGICA_ICA_DEFAULTS.mask.smri_mult, 'threshold', ...
                PGICA_ICA_DEFAULTS.mask.smri_threshold, 'prefix', [pgicaInfo.prefix, '_smri'], 'corr_threshold', ...
                PGICA_ICA_DEFAULTS.mask.smri_corr_threshold, 'outputdir', outputDir, 'disp_corr', 0);
        elseif (strcmpi(modalities{nM}, 'fmri'))
            mask = ica_fuse_generateMask(currentFiles, 'multiplier', PGICA_ICA_DEFAULTS.mask.fmri_mult, 'threshold', ...
                PGICA_ICA_DEFAULTS.mask.fmri_threshold, 'prefix', [pgicaInfo.prefix, '_fmri'], 'corr_threshold', ...
                PGICA_ICA_DEFAULTS.mask.fmri_corr_threshold, 'outputdir', outputDir, 'disp_corr', 0);
        else
            tmp_dat = ica_fuse_loadData(deblank(currentFiles(1, :)));
            mask = (1:length(tmp_dat));
        end
    else
        % Custom mask
        mask = ica_fuse_loadData(mask);
        if (~strcmpi(modalities{nM}, 'smri') && ~strcmpi(modalities{nM}, 'fmri') && ~strcmpi(modalities{nM}, 'dti'))
            mask = squeeze(mask(:, 2));
        end
    end
    
    if (strcmpi(modalities{nM}, 'smri') || strcmpi(modalities{nM}, 'dti') || strcmpi(modalities{nM}, 'fmri'))
        tmpV = ica_fuse_spm_vol(deblank(currentFiles(1,:)));
        tmpV = tmpV(1);
        tmpV.fname = fullfile(outputDir, [pgicaInfo.prefix, '_', modalities{nM}, 'Mask.nii']);
        ica_fuse_write_vol(tmpV, mask);
        structV{nM} = tmpV;
    else
        structV{nM} = VStruct;
    end
    
    mask = find(abs(mask) > eps);
    maskIndices{nM} = mask;
    
end


numComp2 = pgicaInfo.featuresInfo(2).components;
if (length(numComp2) == 1)
    numComp2(2) = numComp2(1);
end
pgicaInfo.featuresInfo(2).components = numComp2;
featuresInfo(2).components = numComp2;

fmri_files = pgicaInfo.featuresInfo(2).files;
pgicaInfo.num_subjects = size(fmri_files, 1);

Timepoints = zeros(1, size(fmri_files, 1));
for nT = 1:length(Timepoints)
    Timepoints(nT) = size(ica_fuse_rename_4d_file(deblank(fmri_files(nT, :))), 1);
end

%% PCA Reduction
pgicaInfo.mask_indices = maskIndices;

data = cell(1, length(modalities));
whiteM = data;
dewhiteM = data;
for nM = 1:length(modalities)
    [data{nM}, whiteM{nM}, dewhiteM{nM}]  = getPCARedFiles(pgicaInfo, nM);
end

dewhiteM_fmri = dewhiteM{2}{1};
dewhiteM_subjects = dewhiteM{2}{2};
dewhiteM{2}{2} = cat(2, dewhiteM_subjects{:});

%% Run parallel ica
if (length(data) == 2)
    dewhiteM = cat(2, dewhiteM{:});
    whiteM = cat(2, whiteM{:});
    [icasig, loading, W] = ica_fuse_calculate_parallelICA(data, dewhiteM, whiteM, 1, 'AT', [], ica_options);
else
    dewhiteM = dewhiteM([1, 3, 2]);
    dewhiteM = cat(2, dewhiteM{:});
    whiteM = whiteM([1, 3, 2]);
    whiteM = cat(2, whiteM{:});
    [icasig, loading, W] = ica_fuse_calculate_parallelICA(data([2, 1, 3]), dewhiteM, whiteM, 1, ...
        'TAA', [], ica_options);
    icasig = icasig([2, 1, 3]);
    loading = loading([2, 1, 3]);
    W = W([2, 1, 3]);
end
fileName = [pgicaInfo.prefix, '_pgica_ica_ica.mat'];
save(fullfile(outputDir, fileName), 'icasig',  'W', 'icasig');


S_group = icasig{2};
num_sub = size(featuresInfo(2).files, 1);
num_pc_fmri_group = featuresInfo(2).components(end);

clear icasig

back_reconstructed_s_fmri_sub_stacked = whiteM{end};
S_sub = zeros(num_sub, num_pc_fmri_group, size(data{2}, 2));
for sub = 1:num_sub%%calculate source for all subjects, ie back reconstruction
    temp = back_reconstructed_s_fmri_sub_stacked(sub,:,:);
    temp_whiteM = reshape(temp,num_pc_fmri_group, size(data{2}, 2));
    S_sub(sub,:,:) = W{2}*temp_whiteM;
end
D2 = zeros(num_sub,num_pc_fmri_group);%D2 is the constructed A for data2 like mx for data1
for i = 1:num_pc_fmri_group
    for j = 1:num_sub
        temp = S_sub(j,:,:);
        temp_S_sub = reshape(temp,num_pc_fmri_group, size(data{2}, 2));
        D2(j,i) = sum((S_group(i,:)-temp_S_sub(i,:)).^2);
    end
end

mixing_coeff{1} = loading{1};
mixing_coeff{2} = D2;
if (length(modalities) == 3)
    mixing_coeff{3} = loading{3};
end

if (length(loading) == 2)
    [corr_matrix, p_matrix] = corr(loading{1}, D2);
    
    [dd, corr_inds] = max(abs(corr_matrix)');
    
    corr_vals = diag(corr(loading{1}, D2(:, corr_inds)));
    
else
    [bestIndex, pair1_corr, pair2_corr, pair3_corr] = ica_fuse_mostCorrelatedTrio(mixing_coeff{1}, mixing_coeff{2}, mixing_coeff{3});
    corr_vals = [pair1_corr(1), pair2_corr(1), pair3_corr(1)];
    corr_inds = bestIndex(1, 1:3);
    
end


%% write feature nifti files
loadingFiles = cell(1, length(featuresInfo));
outputFiles = loadingFiles;
for nM = 1:length(featuresInfo)
    %files = featuresInfo(nM).files;
    V = structV{nM};
    %     if (strcmpi(modalities{nM}, 'smri') || strcmpi(modalities{nM}, 'fmri') || strcmpi(modalities{nM}, 'dti'))
    %         V = ica_fuse_getVol(deblank(files(1,:)));
    %     else
    %         V = VStruct;
    %     end
    V = V(1);
    V.n(1)=1;
    V.dt(1) = 4;
    
    if  (strcmpi(modalities{nM}, 'fmri'))
        VV = repmat(V, num_pc_fmri_group, 1);
    end
    
    [outputFiles{nM}, loadingFiles{nM}] = writeNiftis(V, data{nM}, mixing_coeff{nM}, pgicaInfo.mask_indices{nM}, ...
        [pgicaInfo.prefix, '_feature', num2str(nM)], outputDir, modalities{nM});
    
end

backReconDir = [pgicaInfo.prefix, '_backrecon'];
if (exist(fullfile(outputDir, backReconDir), 'dir') ~= 7)
    mkdir(outputDir, backReconDir);
end

backReconFiles = cell(size(fmri_files, 1), 2);


num_pc_fmri_sub = featuresInfo(2).components(1);
num_pc_fmri_group = featuresInfo(2).components(2);
GA = dewhiteM_fmri*pinv(W{2});

disp('Computing subject back-reconstructed fmri components  ...');

%% write back-reconstructed nifti files for fmri
endT = 0;
minTp = min(Timepoints);
for n = 1:size(fmri_files, 1)
    
    ic = squeeze(S_sub(n, :,:));
    
    if (n == 1)
        meanIC = zeros(size(ic));
    end
    
    meanIC = ic + meanIC;
    
    ic = ic';
    
    tmp = zeros(prod(VV(1).dim(1:3)), size(ic, 2));
    tmp(pgicaInfo.mask_indices{2}, :) = ic;
    
    tmp = reshape(tmp, [VV(1).dim(1:3), length(VV)]);
    
    outputFileName = fullfile(backReconDir, [pgicaInfo.prefix, '_sub', ica_fuse_returnFileIndex(n), '_component_ica_s1_.nii']);
    backReconFiles{n, 1} = outputFileName;
    ica_fuse_write_nifti_data(outputFileName, VV, tmp, '4d nifti');
    clear tmp;
    
    startT = endT  + 1;
    endT = endT + num_pc_fmri_sub;
    
    TC = dewhiteM_subjects{n}*GA(startT:endT, :);
    outputFileName = fullfile(backReconDir, [pgicaInfo.prefix, '_sub', ica_fuse_returnFileIndex(n), '_timecourses_ica_s1_.nii']);
    V.fname = fullfile(outputDir, outputFileName);
    V.dim(1:3) = [size(TC, 1), size(TC, 2), 1];
    
    if (n == 1)
        meanTC = zeros(minTp, size(TC, 2));
    end
    
    meanTC = meanTC + TC(1:minTp, :);
    
    backReconFiles{n, 2} = outputFileName;
    ica_fuse_write_vol(V, TC);
    
end

meanTC = meanTC./size(fmri_files, 1);
meanIC = meanIC./size(fmri_files, 1);


tmp = zeros(prod(VV(1).dim(1:3)), size(meanIC, 1));
tmp(pgicaInfo.mask_indices{2}, :) = meanIC';

tmp = reshape(tmp, [VV(1).dim(1:3), length(VV)]);

outputFileName = fullfile(backReconDir, [pgicaInfo.prefix, '_mean_component_ica_s_all_.nii']);
ica_fuse_write_nifti_data(outputFileName, VV, tmp, '4d nifti');


outputFileName = fullfile(backReconDir, [pgicaInfo.prefix, '_mean_timecourses_ica_s_all_.nii']);
V.fname = fullfile(outputDir, outputFileName);
V.dim(1:3) = [size(meanTC, 1), size(meanTC, 2), 1];
ica_fuse_write_vol(V, meanTC);

pgicaInfo.outputFiles = outputFiles;
pgicaInfo.loadingFiles = loadingFiles;
pgicaInfo.backReconFiles = backReconFiles;
pgicaInfo.corr_vals = corr_vals;
pgicaInfo.corr_inds = corr_inds;

fileName = fullfile(outputDir, pgicaInfo.param_file);
disp(['Saving PGICA fusion parameters in file ', fileName]);
save(fileName, 'pgicaInfo');
disp('Done');
fprintf('\n');


diary('off');


% Perform some additional steps to mimic GIFT parameter settings to run
% post-processing tools like mancovan, dfnc, etc on the reconstructed subject components.
disp('Performing some additional steps to use subject back-reconstructed fmri components in the GIFT toolbox ...');


global PICA_INPUT_FMRI_FILES;
global PICA_OUTPUT_FILE_DIR;
global PICA_INPUT_PREFIX;
global PICA_COMPS;
global PICA_INPUT_MASK;

PICA_INPUT_FMRI_FILES = cellstr(fmri_files);
PICA_OUTPUT_FILE_DIR = backReconDir;
PICA_INPUT_PREFIX = pgicaInfo.prefix;
PICA_COMPS = pgicaInfo.featuresInfo(2).components;
PICA_INPUT_MASK = fullfile(outputDir, [pgicaInfo.prefix, '_fmriMask.nii']);
copyfile(fullfile(fileparts(which('fusion.m')), 'ica_fuse_batch_files', 'input_gica_pica.m'), backReconDir);
ica_param_file = icatb_read_batch_file(fullfile(backReconDir, 'input_gica_pica.m'));
load(ica_param_file);
icatb_runAnalysis(sesInfo, 2);



cd(outputDir);


function [outputFiles, loadingFiles] = writeNiftis(V, S, A, mask1, prefix, outputDir, modality)
%% Write Niftis


outputFiles = cell(1, size(S, 1));
for nS = 1:size(S, 1)
    if (strcmpi(modality, 'smri') || strcmpi(modality, 'fmri'))
        fileName = [prefix, '_comp_', ica_fuse_returnFileIndex(nS), '.nii'];
        tmp = zeros(V.dim(1:3));
        tmp(mask1) = squeeze(S(nS, :));
        V.fname = fullfile(outputDir, fileName);
        ica_fuse_write_vol(V, tmp);
    else
        tmp = squeeze(S(nS, :));
        fileName = [prefix, '_comp_', ica_fuse_returnFileIndex(nS), '.asc'];
        save(fullfile(outputDir, fileName), 'tmp', '-ascii');
    end
    
    outputFiles{nS} = fileName;
end


fileName = [prefix, '_loading_coeff.nii'];
V.dim(1:3) = [size(A, 1), size(A, 2), 1];
loadingFiles = fileName;
V.fname = fileName;
ica_fuse_write_vol(V, A);


function varargout = getPCARedFiles(pgicaInfo, nM)
%% Get PCA Reduction files

modality = lower(pgicaInfo.modalities{nM});
featuresInfo = pgicaInfo.featuresInfo(nM);
featureName = featuresInfo.name;
feaInfo.files = featuresInfo.files;
feaInfo.modality = modality;
feaInfo.feature_name = featureName;
mask(1).ind = pgicaInfo.mask_indices{nM};
outputDir = pgicaInfo.outputDir;
components = featuresInfo.components;
pcaDir = [pgicaInfo.prefix, '_pca'];
if (exist(fullfile(outputDir, pcaDir), 'dir') ~= 7)
    mkdir(outputDir, pcaDir);
end
disp(['Loading ', featureName, ' data ...']);

if (~strcmpi(modality, 'fmri'))
    %% Load smri data
    smri_data = ica_fuse_applyMask(feaInfo, mask);
    smri_data = smri_data.data;
    
    if (~strcmpi(modality, 'gene'))
        % don't remove mean for gene data
        smri_data = ica_fuse_remove_mean(smri_data', 0)';
    end
    
    [V, Lambda, whitesig, whiteM, dewhiteM] = ica_fuse_calculate_pca(smri_data, components(1));
    fileName = fullfile(pcaDir, [pgicaInfo.prefix, '_pgica_ica_pca_feature', num2str(nM), '.mat']);
    save(fullfile(outputDir, fileName), 'V', 'Lambda', 'whitesig', 'whiteM', 'dewhiteM');
    
    varargout{1} = whitesig;
    varargout{2}  = whiteM;
    varargout{3}  = dewhiteM;
    
else
    % fmri
    
    whitesig_g = zeros(components(1)*size(feaInfo.files, 1), length(mask.ind));
    dewhiteM_g = cell(size(feaInfo.files, 1), 1);
    whiteM_g = dewhiteM_g;
    endT = 0;
    %mask(1).ind = mask2(:);
    Timepoints = zeros(1, size(feaInfo.files, 1));
    for nF = 1:length(Timepoints)
        disp(['Computing pca on subject ', num2str(nF), ' ...']);
        feaInfo.files =  deblank(featuresInfo.files(nF, :));
        tmp = ica_fuse_applyMask(feaInfo, mask);
        tmp = tmp.data;
        Timepoints(nF) = size(tmp, 1);
        [V, Lambda, whitesig, whiteM, dewhiteM] = ica_fuse_calculate_pca(ica_fuse_remove_mean(tmp', 0)', components(1));
        fileName = fullfile(pcaDir, [pgicaInfo.prefix, '_pgica_ica_pca_feature', num2str(nM), '_sub', ica_fuse_returnFileIndex(nF), '.mat']);
        save(fullfile(outputDir, fileName), 'V', 'Lambda', 'whitesig', 'whiteM', 'dewhiteM');
        startT = endT + 1;
        endT = endT + components(1);
        whitesig_g(startT:endT, :) = whitesig;
        dewhiteM_g{nF} = dewhiteM;
        whiteM_g{nF} = whiteM;
        clear V Lambda whiteM dewhiteM tmp whitesig;
    end
    
    fprintf('\n');
    disp(['Computing pca on stacked components across subjects ...']);
    [V, Lambda, whitesig, whiteM, dewhiteM] = ica_fuse_calculate_pca(whitesig_g, components(2));
    fileName = fullfile(pcaDir, [pgicaInfo.prefix, '_pgica_ica_pca_feature2', num2str(nM), '.mat']);
    save(fullfile(outputDir, fileName), 'V', 'Lambda', 'whitesig', 'whiteM', 'dewhiteM');
    
    whitesig_fmri = whitesig;
    whiteM_fmri = whiteM;
    dewhiteM_fmri = dewhiteM;
    
    %dewhiteM_subjects = dewhiteM_g;
    %dewhiteM_g = cat(2, dewhiteM_g{:});
    whiteM_g = cat(1, whiteM_g{:});
    
    num_pc_fmri_sub = components(1);
    num_pc_fmri_group = components(2);
    %num_sub = length(Timepoints);
    
    dewhiteM_fmri_group_cell = mat2cell(dewhiteM, components(1)*ones(1, length(Timepoints)), components(2));
    clear V Lambda whiteM dewhiteM tmp whitesig;
    
    back_reconstructed_s_fmri_sub_stacked = zeros(length(Timepoints), num_pc_fmri_group, size(whitesig_fmri, 2));
    for sub = 1:length(Timepoints)
        temp_M_sub = pinv(dewhiteM_fmri_group_cell{sub,1})*whitesig_g((sub-1)*num_pc_fmri_sub+1:sub*num_pc_fmri_sub, :);
        back_reconstructed_s_fmri_sub_stacked(sub,:,:) = temp_M_sub;
    end
    
    %pgicaInfo.featuresInfo(nM).timepoints = Timepoints;
    varargout{1} = whitesig_fmri;
    varargout{2} = {whiteM_fmri, whiteM_g, back_reconstructed_s_fmri_sub_stacked};
    varargout{3} = {dewhiteM_fmri, dewhiteM_g};
    %varargout{4} = Timepoints;
end

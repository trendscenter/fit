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



numComp1 = pgicaInfo.numComp1;
numComp2 = pgicaInfo.numComp2;


if (~isfield(pgicaInfo, 'ica_options'))
    ica_options = ica_fuse_paraICAOptions(min([numComp1, numComp2]), 'off');
    pgicaInfo.ica_options = ica_options;
end

ica_options = pgicaInfo.ica_options;

results_file = fullfile(outputDir, [pgicaInfo.prefix, '_pgica_ica.log']);
diary(results_file);

%% Step 1: Do pca and get the whitening dewhitening information
mask1 = pgicaInfo.maskFile_modality1;
mask2 = pgicaInfo.maskFile_modality2;


if (isempty(mask1))
    disp('Calculating mask for modality 1 ...');
    mask1 = ica_fuse_generateMask(pgicaInfo.files{1}, 'multiplier', PGICA_ICA_DEFAULTS.mask.smri_mult, 'threshold', ...
        PGICA_ICA_DEFAULTS.mask.smri_threshold, 'prefix', [pgicaInfo.prefix, '_smri'], 'corr_threshold', ...
        PGICA_ICA_DEFAULTS.mask.smri_corr_threshold, 'outputdir', outputDir, 'disp_corr', 0);
    
    disp('Calculating mask for modality 2 ...');
    mask2 = ica_fuse_generateMask(pgicaInfo.files{2}, 'multiplier', PGICA_ICA_DEFAULTS.mask.fmri_mult, 'threshold', ...
        PGICA_ICA_DEFAULTS.mask.fmri_threshold, 'prefix', [pgicaInfo.prefix, '_smri'], 'corr_threshold', ...
        PGICA_ICA_DEFAULTS.mask.fmri_corr_threshold, 'outputdir', outputDir, 'disp_corr', 0);
    
else
    
    mask1 = ica_fuse_loadData(mask1);
    mask2 = ica_fuse_loadData(mask2);
    
end


VStruct = ica_fuse_spm_vol(deblank(pgicaInfo.files{1}(1,:)));
VStruct = VStruct(1);
VStruct.fname = fullfile(outputDir, [pgicaInfo.prefix, '_smriMask.nii']);
VFunc = ica_fuse_spm_vol(deblank(pgicaInfo.files{2}(1,:)));
VFunc = VFunc(1);
VFunc.fname = fullfile(outputDir, [pgicaInfo.prefix, '_fmriMask.nii']);
VStruct.dt(1) = 4;
VFunc.dt(1) = 4;

ica_fuse_write_vol(VStruct, mask1);
ica_fuse_write_vol(VFunc, mask2);

mask1 = find(abs(mask1) > eps);
mask2 = find(abs(mask2) > eps);

pcaDir = [pgicaInfo.prefix, '_pca'];

if (exist(fullfile(outputDir, pcaDir), 'dir') ~= 7)
    mkdir(outputDir, pcaDir);
end


%% Load smri data
disp('Loading modality 1 data ...');
featureInfo.files = pgicaInfo.files{1};
featureInfo.modality = 'smri';
featureInfo.feature_name = 'Feature';
mask(1).ind = mask1(:);
smri_data = ica_fuse_applyMask(featureInfo, mask);
smri_data = smri_data.data;
[V, Lambda, whitesig, whiteM, dewhiteM] = ica_fuse_calculate_pca(smri_data, numComp1);
fileName = fullfile(pcaDir, [pgicaInfo.prefix, '_pgica_ica_pca_feature1.mat']);
save(fullfile(outputDir, fileName), 'V', 'Lambda', 'whitesig', 'whiteM', 'dewhiteM');

whitesig_smri = whitesig;
dewhiteM_smri = dewhiteM;
whiteM_smri = whiteM;
clear V Lambda whiteM dewhiteM smri_data featureInfo;


disp(['Loading modality 2 data ...']);
files = pgicaInfo.files{2};
featureInfo.modality = 'fmri';

whitesig_g = zeros(numComp2(1)*size(files, 1), length(mask2));
dewhiteM_g = cell(size(files, 1), 1);
whiteM_g = dewhiteM_g;
endT = 0;
mask(1).ind = mask2(:);
Timepoints = zeros(1, size(files, 1));
for nF = 1:size(files, 1)
    disp(['Computing pca on subject ', num2str(nF), ' ...']);
    featureInfo.files = deblank(files(nF, :));
    featureInfo.feature_name = 'Feature';
    tmp = ica_fuse_applyMask(featureInfo, mask);
    tmp = tmp.data;
    Timepoints(nF) = size(tmp, 1);
    [V, Lambda, whitesig, whiteM, dewhiteM] = ica_fuse_calculate_pca(tmp, numComp2(1));
    fileName = fullfile(pcaDir, [pgicaInfo.prefix, '_pgica_ica_pca_feature2_sub', ica_fuse_returnFileIndex(nF), '.mat']);
    %     if (exist(fullfile(outputDir, 'sub'), 'dir') ~a= 7)
    %         mkdir(outputDir, [pgicaInfo.prefix, '_pca']);
    %     end
    save(fullfile(outputDir, fileName), 'V', 'Lambda', 'whitesig', 'whiteM', 'dewhiteM');
    startT = endT + 1;
    endT = endT + numComp2(1);
    whitesig_g(startT:endT, :) = whitesig;
    dewhiteM_g{nF} = dewhiteM;
    whiteM_g{nF} = whiteM;
    clear V Lambda whiteM dewhiteM tmp whitesig;
end

fprintf('\n');
disp(['Computing pca on stacked components across subjects ...']);
[V, Lambda, whitesig, whiteM, dewhiteM] = ica_fuse_calculate_pca(whitesig_g, numComp2(2));
fileName = fullfile(pcaDir, [pgicaInfo.prefix, '_pgica_ica_pca_feature2.mat']);
save(fullfile(outputDir, fileName), 'V', 'Lambda', 'whitesig', 'whiteM', 'dewhiteM');

whitesig_fmri = whitesig;
whiteM_fmri = whiteM;
dewhiteM_fmri = dewhiteM;

dewhiteM_subjects = dewhiteM_g;
dewhiteM_g = cat(2, dewhiteM_g{:});
whiteM_g = cat(1, whiteM_g{:});

num_pc_fmri_sub = numComp2(1);
num_pc_fmri_group = numComp2(2);
num_sub = size(files, 1);


pgicaInfo.num_subjects = num_sub;

dewhiteM_fmri_group_cell = mat2cell(dewhiteM, numComp2(1)*ones(1, size(files, 1)), numComp2(2));
clear V Lambda whiteM dewhiteM tmp whitesig;

back_reconstructed_s_fmri_sub_stacked = zeros(size(files, 1), num_pc_fmri_group, size(whitesig_fmri, 2));
for sub = 1:size(files, 1)
    temp_M_sub = pinv(dewhiteM_fmri_group_cell{sub,1})*whitesig_g((sub-1)*num_pc_fmri_sub+1:sub*num_pc_fmri_sub, :);
    back_reconstructed_s_fmri_sub_stacked(sub,:,:) = temp_M_sub;
end

data = {whitesig_smri, whitesig_fmri};
whiteM = {whiteM_smri, whiteM_fmri, whiteM_g, back_reconstructed_s_fmri_sub_stacked};
dewhiteM = {dewhiteM_smri, dewhiteM_fmri, dewhiteM_g};

%% Run parallel ica
[icasig, loading, W] = ica_fuse_calculate_parallelICA(data, dewhiteM, whiteM, 1, 'AT', [], ica_options);
fileName = [pgicaInfo.prefix, '_pgica_ica_ica.mat'];
save(fullfile(outputDir, fileName), 'icasig',  'W', 'icasig');


S_group = icasig{2};

clear icasig

back_reconstructed_s_fmri_sub_stacked = whiteM{4};
S_sub = zeros(num_sub, num_pc_fmri_group, size(whitesig_fmri, 2));
for sub = 1:num_sub%%calculate source for all subjects, ie back reconstruction
    temp = back_reconstructed_s_fmri_sub_stacked(sub,:,:);
    temp_whiteM = reshape(temp,num_pc_fmri_group, size(whitesig_fmri, 2));
    S_sub(sub,:,:) = W{2}*temp_whiteM;
end
D2 = zeros(num_sub,num_pc_fmri_group);%D2 is the constructed A for data2 like mx for data1
for i = 1:num_pc_fmri_group
    for j = 1:num_sub
        temp = S_sub(j,:,:);
        temp_S_sub = reshape(temp,num_pc_fmri_group, size(whitesig_fmri, 2));
        D2(j,i) = sum((S_group(i,:)-temp_S_sub(i,:)).^2);
    end
end

[corr_matrix, p_matrix] = corr(loading{1}, D2);

[dd, corr_inds] = max(abs(corr_matrix)');


corr_vals = diag(corr(loading{1}, D2(:, corr_inds)));

mixing_coeff{1} = loading{1};
mixing_coeff{2} = D2;


%% write feature nifti files
smri_files = pgicaInfo.files{1};
fmri_files = pgicaInfo.files{2};
V = ica_fuse_getVol(deblank(smri_files(1,:)));
V = V(1);
V.n(1)=1;
V.dt(1) = 4;

[fileSet1, loadingFiles{1}] = writeNiftis(V, whitesig_smri, mixing_coeff{1}, mask1, [pgicaInfo.prefix, '_feature1'], outputDir);

V = ica_fuse_getVol(deblank(fmri_files(1,:)));
V = V(1);
V.dt(1) = 4;
[fileSet2, loadingFiles{2}] = writeNiftis(V, S_group, mixing_coeff{2}, mask2, [pgicaInfo.prefix, '_feature2'], outputDir);

VV = repmat(V, num_pc_fmri_group, 1);


outputFiles{1} = fileSet1;
outputFiles{2} = fileSet2;


backReconDir = [pgicaInfo.prefix, '_backrecon'];
if (exist(fullfile(outputDir, backReconDir), 'dir') ~= 7)
    mkdir(outputDir, backReconDir);
end

backReconFiles = cell(size(fmri_files, 1), 2);


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
    tmp(mask2, :) = ic;
    
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
tmp(mask2, :) = meanIC';

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
PICA_COMPS = numComp2;
PICA_INPUT_MASK = fullfile(outputDir, [pgicaInfo.prefix, '_fmriMask.nii']);
copyfile(fullfile(fileparts(which('fusion.m')), 'ica_fuse_batch_files', 'input_gica_pica.m'), backReconDir);
ica_param_file = icatb_read_batch_file(fullfile(backReconDir, 'input_gica_pica.m'));
load(ica_param_file);
icatb_runAnalysis(sesInfo, 2);



cd(outputDir);


function [outputFiles, loadingFiles] = writeNiftis(V, S, A, mask1, prefix, outputDir)
%% Write Niftis


outputFiles = cell(1, size(S, 1));
for nS = 1:size(S, 1)
    fileName = [prefix, '_comp_', ica_fuse_returnFileIndex(nS), '.nii'];
    tmp = zeros(V.dim(1:3));
    tmp(mask1) = squeeze(S(nS, :));
    V.fname = fullfile(outputDir, fileName);
    outputFiles{nS} = fileName;
    ica_fuse_write_vol(V, tmp);
end


fileName = [prefix, '_loading_coeff.nii'];
V.dim(1:3) = [size(A, 1), size(A, 2), 1];
loadingFiles = fileName;
V.fname = fileName;
ica_fuse_write_vol(V, A);
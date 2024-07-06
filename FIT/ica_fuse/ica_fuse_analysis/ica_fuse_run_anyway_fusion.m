function ica_fuse_run_anyway_fusion(fusionInfo)
%% Run any way fusion
%
% Inputs:
% fusionInfo - Data structure used in fusion anaysis

getFeatureInfo = 0;
if (~exist('fusionInfo', 'var'))
    param_file = ica_fuse_selectEntry('title', 'Select anyway fusion Info file', 'filter', '*anyway_fusion.mat');
    drawnow;
    load(param_file);
    getFeatureInfo = 1;
end


outputDir = fusionInfo.setup_analysis.outputDir;
prefix = fusionInfo.setup_analysis.prefix;
out_file = fullfile(outputDir, [prefix, '_anyway_fusion.mat']);

scale_components = 'z-scores';
try
    scale_components = fusionInfo.setup_analysis.scale_components;
catch
end

% record analysis information in a log file
diaryFile = fullfile(outputDir, [prefix, '_anyway_results.log']);
diary(diaryFile);

tic;

fprintf('\n');
disp('Starting Analysis ...');

if (getFeatureInfo)
    
    feature_info = repmat(struct('feature_name', '', 'modality', '', 'files', ''), 1, length(fusionInfo.setup_analysis.dataInfo(1).feature));
    for nFeature = 1:length(fusionInfo.setup_analysis.dataInfo(1).feature)
        feature_info(nFeature).feature_name = fusionInfo.setup_analysis.dataInfo(1).feature(nFeature).name;
        feature_info(nFeature).modality = fusionInfo.setup_analysis.dataInfo(1).feature(nFeature).modality;
        tmp_files = cell(length(fusionInfo.setup_analysis.dataInfo), 1);
        for nG = 1:length(tmp_files)
            tmp_files{nG} = char(fusionInfo.setup_analysis.dataInfo(nG).feature(nFeature).files.name);
        end
        feature_info(nFeature).files = ica_fuse_rename_4d_file(char(tmp_files));
        feature_info(nFeature).comp = fusionInfo.setup_analysis.numComp(nFeature);
    end
    
else
    
    feature_info = fusionInfo.setup_analysis.feature_info;
    
end

fusionInfo.run_analysis.feature_info = feature_info;


num_comps = [feature_info.comp];

fusionInfo.run_analysis.num_comps = num_comps;

lambda_iva = 0.2;
try
    lambda_iva = fusionInfo.setup_analysis.algorithm_opts.lambda_iva;
catch
end


numFeatures = length(feature_info);
mask_ind = repmat(struct('ind', []), 1, numFeatures);
timeAxis = repmat(struct('data', []), 1, numFeatures);
%dataInfo = repmat(struct('name', [], 'feature_name', []), 1, numFeatures);

%% Mask
if ~isfield(fusionInfo.setup_analysis, 'mask_ind')
    disp('Creating mask ...');
    for nF = 1:length(feature_info)
        maskFile = '';
        try
            maskFile = feature_info(nF).mask;
        catch
        end
        dataInfo.feature.files.name = feature_info(nF).files;
        dataInfo.feature.name = feature_info(nF).feature_name;
        dataInfo.feature.modality = feature_info(nF).modality;
        maskIndices = ica_fuse_createMask(dataInfo, maskFile);
        mask_ind(nF).ind = maskIndices.ind;
    end
    disp('Done creating mask');
else
    mask_ind = fusionInfo.setup_analysis.mask_ind;
end

%% Load features
data = cell(1, length(feature_info));
%dataInfo = [];
for nF = 1:length(feature_info)
    disp(['Loading feature ', num2str(feature_info(nF).feature_name), ' ...']);
    tmp = ica_fuse_applyMask(feature_info(nF), mask_ind(nF));
    tmp.data = ica_fuse_normalize_data(tmp.data, 1);
    if (strcmpi(feature_info(nF).modality, 'eeg'))
        tmp.data = diag(1./std(tmp.data, [], 2))*tmp.data;
    end
    if (~strcmpi(feature_info(nF).modality, 'gene'))
        tmp.data = tmp.data - repmat(mean(tmp.data, 2), 1, size(tmp.data, 2));
    end
    data{nF} = tmp.data;
    timeAxis(nF).data = tmp.xAxis;
    dataInfo.feature(nF).name = deblank(feature_info(nF).feature_name);
    dataInfo.feature(nF).files.name = deblank(feature_info(nF).files(1, :));
    dataInfo.feature(nF).modality =  feature_info(nF).modality;
end

disp('Computing anyway fusion ...');

[~, ~, icasig, A] = ica_fuse_anywayica_nmod(data, 'ncomps', num_comps, 'iva_lambda', lambda_iva);

fusionInfo.run_analysis.numSubjects = size(data{1}, 1);

% ICA file
ica_file = fullfile(outputDir, [prefix, '_anyway_ica.mat']);
save(ica_file, 'icasig', 'A');

disp('Done');


fusionInfo.run_analysis.mask_ind = mask_ind;


fusionInfo.run_analysis.scale_components = scale_components;

if (scale_components)
    %% Scale components to z-scores
    for nF = 1:length(feature_info)
        tmp_ic = icasig{nF};
        tmp_tc = A{nF}';
        tmp_tc = ica_fuse_convertToZScores(tmp_tc)';
        tmp_ic = ica_fuse_convertToZScores(tmp_ic);
        A{nF} = tmp_tc;
        icasig{nF} = tmp_ic;
    end
end

%% Write components to file
outFileNaming = [prefix, '_anyway_'];
outputFiles = ica_fuse_write_parallel_ica_comp(icasig, A, dataInfo, mask_ind, outFileNaming, outputDir,  size(A{1}, 1));

fusionInfo.run_analysis.outputFiles = outputFiles;

fusionInfo.run_analysis.isInitialized = 1;


% Save fusion file
save(out_file, 'fusionInfo');

t_end = toc;

disp(['Time taken to run the analysis is ', num2str(t_end), ' seconds']);
fprintf('\n');
disp(['All the analysis information is stored in log file: ', diaryFile]);
fprintf('\n');

diary('off');
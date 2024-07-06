function ica_fuse_anyway_fusion_batch(inputFile)
%% Batch script for running any way fusion
%


disp(['Reading batch file ', inputFile]);

inputData = ica_fuse_eval_script(inputFile);

outputDir = pwd;

try
    outputDir = inputData.outputDir;
catch
end

prefix = inputData.prefix;

modality = cellstr(inputData.modality);

featureNames = strcat('Feature ', num2str((1:length(modality))'));

try
    featureNames = inputData.featureNames;
catch
end

featureNames = cellstr(featureNames);

files = inputData.files;

if (length(files) ~= length(modality))
    error('Input files cell array length doesn''t match the number of modalities');
end

if (length(featureNames) ~= length(modality))
    error('featureNames cell array length doesn''t match the number of modalities');
end

numComp = inputData.numComp;
if (length(numComp) ~= length(modality))
    error('numComp vector length doesn''t match the number of modalities');
end

maskFile = [];
try
    maskFile = inputData.maskFile;
catch
end

if (~isempty(maskFile))
    if (length(maskFile) ~= length(modality))
        error('maskFile cell array length doesn''t match the number of modalities');
    end
end


feature_info = repmat(struct('feature_name', '', 'files', [], 'modality', '', 'comp', '', 'mask', ''), 1, length(modality));
for nModality = 1:length(modality)
    feature_info(nModality).feature_name = featureNames{nModality};
    feature_info(nModality).comp = numComp(nModality);
    feature_info(nModality).modality = modality{nModality};
    if (~isempty(maskFile))
        feature_info(nModality).mask = maskFile{nModality};
    end
    
    tmp = files{nModality};
    if (size(tmp, 1) == 1)
        [tmp_path, tmpF, extn] = fileparts(tmp);
        tmp = ica_fuse_listFiles_inDir(tmp_path, [tmpF, extn]);
        if isempty(tmp)
            error(['Files not found for (', fullfile(tmp, [tmpF, extn]), ')']);
        end
        tmp = ica_fuse_fullFile('files', tmp, 'directory', tmp_path);
    end
    
    tmp = ica_fuse_rename_4d_file(tmp);
    feature_info(nModality).files = tmp;
    
end

scale_components = 'z-scores';
try
    scale_components = inputData.scale_components;
catch
end

fusionInfo.setup_analysis.scale_components = scale_components;

try
    fusionInfo.setup_analysis.algorithm_opts = inputData.algorithm_opts;
catch
end

fusionInfo.setup_analysis.feature_info = feature_info;
fusionInfo.setup_analysis.prefix = prefix;
fusionInfo.setup_analysis.outputDir = outputDir;
out_file = fullfile(outputDir, [prefix, '_anyway_fusion.mat']);
save(out_file, 'fusionInfo');

disp(['Fusion information is saved in ', out_file]);
fprintf('\n');


%% Run analysis
ica_fuse_run_anyway_fusion(fusionInfo);


%% Display results

try
    display_results = inputData.display_results;
catch
end

if (exist('display_results', 'var'))
    ica_fuse_anyway_fusion_summary(out_file, display_results);
end
function ica_fuse_setup_neural_net_fusion(multi_param_file)
%% Neural network based approach is used to fuse smri and fmri modalities
% This toolbox assumes you have run ICA analysis using SBM and dfnc
% analysis using GIFT

if (~exist('multi_param_file', 'var'))
    multi_param_file = ica_fuse_selectEntry('typeEntity', 'file', 'filter', '*ica*param*mat;*neural*net*fusion*mat', 'title', 'Select SBM ICA Parameter file/Neural net fusion file');
    drawnow;
end

if (isempty(multi_param_file))
    error('Select neural net fusion file/sbm ica parameter file');
end

load(multi_param_file);

if (exist('sesInfo', 'var'))
    if (~strcmpi(sesInfo.modality, 'smri'))
        error('You need to use SBM ica parameter file');
    end
    neuralNetFusionInfo.smri_param_file = multi_param_file;
    dfnc_param_file = ica_fuse_selectEntry('typeEntity', 'file', 'filter', '*dfnc*mat', 'title', 'Select dFNC parameter file');
    drawnow;
    if (isempty(dfnc_param_file))
        error('dFNC parameter file is not selected');
    end
    neuralNetFusionInfo.dfnc_param_file = dfnc_param_file;
else
    if (~exist('neuralNetFusionInfo', 'var'))
        error('Selected file is not a valid neural network fusion file');
    end
end


% select output analysis directory
outputDir = ica_fuse_selectEntry('typeEntity', 'directory', 'title', 'Select neural net fusion output analysis directory');
if (isempty(outputDir))
    error('Output directory is not selected');
end

neuralNetFusionInfo.outputDir = outputDir;

prefix = 'nnet';
try
    prefix = neuralNetFusionInfo.prefix;
catch
end

lrate = 0.01;
try
    lrate = neuralNetFusionInfo.lrate;
catch
end

max_epochs = 100;
try
    max_epochs = neuralNetFusionInfo.max_epochs;
catch
end

num_runs = 100;
try
    num_runs = neuralNetFusionInfo.num_runs;
catch
end

numParameters = 1;

inputText(numParameters).promptString = 'Enter output files prefix';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = prefix;
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'prefix';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Enter learning rate';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(lrate);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'lrate';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;


numParameters = numParameters + 1;


inputText(numParameters).promptString = 'No. of runs';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(num_runs);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'num_runs';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;


inputText(numParameters).promptString = 'No. of epochs';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(max_epochs);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'max_epochs';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;


answer = ica_fuse_inputDialog('inputtext', inputText, 'Title', 'Select nnet fusion params', 'handle_visibility',  'on', 'windowStyle', 'modal');

if (~isempty(answer))
    prefix = answer{1};
    lrate = answer{2};
    num_runs = answer{3};
    max_epochs = answer{4};
else
    error('Input parameters are not selected');
end

neuralNetFusionInfo.prefix = prefix;
neuralNetFusionInfo.lrate = lrate;
neuralNetFusionInfo.num_runs = num_runs;
neuralNetFusionInfo.max_epochs = max_epochs;

%% Save info
fusionFile = fullfile(outputDir, [prefix, '_neural_net_fusion.mat']);
save(fusionFile, 'neuralNetFusionInfo');


%% Run
cd(outputDir);

load(neuralNetFusionInfo.smri_param_file);

%% Writing smri loadings to csv file
load(fullfile(fileparts(neuralNetFusionInfo.smri_param_file), [sesInfo.calibrate_components_mat_file, '1-1.mat']), 'tc');
neuralNetFusionInfo.features(1).file = [prefix, '_neural_net_fusion_smri_file.csv'];
fileN1 = fullfile(outputDir, neuralNetFusionInfo.features(1).file);
csvwrite(fileN1, tc);

%% Writing dfnc states to csv file
load(neuralNetFusionInfo.dfnc_param_file);
post_process_file = fullfile(fileparts(neuralNetFusionInfo.dfnc_param_file),  [dfncInfo.prefix, '_post_process.mat']);
load(post_process_file);
Nwin = length(clusterInfo.IDXall) / length(dfncInfo.outputFiles);
states = reshape(clusterInfo.IDXall, length(dfncInfo.outputFiles), Nwin);
neuralNetFusionInfo.features(2).file = [prefix, '_neural_net_fusion_dfnc_file.csv'];
fileN2 = fullfile(outputDir, neuralNetFusionInfo.features(2).file);
csvwrite(fileN2, states);

%% write centroids
neuralNetFusionInfo.features(2).centroids = [prefix, '_neural_net_fusion_dfnc_centroids.mat'];
C_all = clusterInfo.Call;
fileN3 = fullfile(outputDir, neuralNetFusionInfo.features(2).centroids);
save(fileN3, 'C_all');

alignment_outputDir = [prefix, '_alignments'];

neuralNetFusionInfo.alignment_outputDir = alignment_outputDir;

fusionFile = [prefix, '_neural_net_fusion.mat'];
disp(['.... saving neural net fusion information in file ', fusionFile]);
save(fullfile(outputDir, fusionFile), 'neuralNetFusionInfo');
disp('Done');
fprintf('\n');

%dist_file = fullfile(fileparts(which('fusion.m')), 'dist', 'nmt_dfnc_correlation_alldata');
disp('Computing alignment scores ...');
dist_file = fullfile(fileparts(which('fusion.m')), 'ica_fuse_bin', 'src', 'nmt_dfnc_correlation_alldata.py');
commandStr = ['python "', dist_file, '" --smri "', fileN1, '" --dfnc "', fileN2, '" --states "', fileN3, '" --lrate ', num2str(lrate), ' --max_epochs ', num2str(max_epochs), ' --num_runs ', ...
    num2str(num_runs), ' --outdir "', fullfile(outputDir, alignment_outputDir), '"'];

disp(commandStr);
[status, message ] = system(commandStr);

if (status)
    error(message);
end

alignmentInfo = cell(1, neuralNetFusionInfo.num_runs);
for nRun = 1:length(alignmentInfo)
    alignments = ica_fuse_read_alignment_scores(fullfile(neuralNetFusionInfo.outputDir, neuralNetFusionInfo.alignment_outputDir), nRun);
    alignmentInfo{nRun} = alignments;
    subj_align = subject_alignment(alignments, states);
    S = get_state_alignments(alignments, states);
    S = cell2mat(cellfun(@(x) mean(x,1), S, 'UniformOutput', false));
    if (nRun == 1)
        states_all = nan(size(S, 1), size(S, 2), length(alignmentInfo));
        subject_alignments = nan(size(subj_align, 1), size(subj_align, 2), size(subj_align, 3), length(alignmentInfo));
    end
    states_all(:, :, nRun) = S;
    subject_alignments(:, :, :, nRun) = subj_align;
end

states_all = squeeze(ica_fuse_nanmean(states_all, 3));
subject_alignments = squeeze(ica_fuse_nanmean(subject_alignments, 4));
subject_alignments = permute(subject_alignments, [3, 1, 2]);
alignmentFile = [prefix, '_alignments.mat'];
neuralNetFusionInfo.alignment_file = alignmentFile;
disp(['.... saving alignment scores in file ', alignmentFile]);
save(fullfile(outputDir, alignmentFile), 'alignmentInfo', 'subject_alignments', 'states_all');

save(fullfile(outputDir, fusionFile), 'neuralNetFusionInfo');

fprintf('Done\n\n');

function A = subject_alignment(Algn, Dfnc)

nSub = size(Dfnc,1);
%nState = 5;
nState = max(Dfnc(:));
nComp = size(Algn,2);
A = nan(nState, nComp, nSub);
for i = 1:nSub
    for j = 1:nState
        idx = find(Dfnc(i,:)==j);
        if ~isempty(idx)
            A(j, :, i) = mean(Algn(idx, :,i),1);
        end
    end
end

function S = get_state_alignments(Algn, Dfnc)
nSub = size(Algn,3);
nTime = size(Dfnc, 2);
nState = max(Dfnc(:));

S = cell(nState,1);

for i = 1:nSub
    for j = 1:nTime
        S{Dfnc(i,j)}(end+1,:) = Algn(j,:, i);
    end
end

function fusionInfo = ica_fuse_calculateICA(fusionInfo)
%% Calculate ICA based on the ICA algorithm. The information will be stored
% in ICA MAT files. The MAT file name is dependent on global variable
% ICA_FILE (See ica_fuse_defaults.m) and combination number.
%
% Inputs:
% fusionInfo - Fusion Info data structure
%
% Outputs:
% fusionInfo - Fusion Info data structure with additional field icaFiles
%

%% Load defaults
ica_fuse_defaults;
global ICA_FILE;
global NUM_RUNS_ICA;
global DATA_REDUCTION_FILE;
global OPEN_DISPLAY_WINDOW;


if (~isfield(fusionInfo.run_analysis, 'numICARuns'))
    if NUM_RUNS_ICA < 1
        numRuns = 1;
    else
        numRuns = ceil(NUM_RUNS_ICA);
    end
else
    numRuns = fusionInfo.run_analysis.numICARuns;
end


% Output directory
outputDir = fusionInfo.run_analysis.outputDir;
numSubjects = fusionInfo.run_analysis.numSubjects;
ica_algorithm = fusionInfo.run_analysis.algorithm;
combinationName = fusionInfo.run_analysis.currentCombName;
comb_number = fusionInfo.run_analysis.currentComb;

ica_algo = ica_fuse_icaAlgorithm;
ica_algo = cellstr(ica_algo);
algorithmName = ica_algo{fusionInfo.run_analysis.algorithm};

try
    analysisType = fusionInfo.run_analysis.type_ica;
catch
    analysisType = 'average';
end


if (strcmpi(algorithmName, 'iva-g') || strcmpi(algorithmName, 'iva-ggd'))
    % For now use single run. MST is required for stable run estimates
    analysisType = 'average';
    numRuns = 1;
end

%% Load PCA file
pcaFile = fullfile(outputDir, fusionInfo.run_analysis.pcaFiles(comb_number).name);
load(pcaFile, 'whitesig', 'dewhiteM');

%% ICA Options
ICA_Options = {};
if isfield(fusionInfo.run_analysis, 'ICA_Options')
    ICA_Options = fusionInfo.run_analysis.ICA_Options;
end

numOfIC = size(whitesig, 1);

disp('-----------------------------------------------------------------------------------------------');
disp(['Doing ICA/IVA for ',  combinationName]);
disp('-----------------------------------------------------------------------------------------------');

fprintf('\n');
disp(['Number of times ICA/IVA will run is ', num2str(numRuns)]);

%% For Fast ICA use eigen vectors and eigen values from PCA step
if (strcmpi(algorithmName, 'fast ica'))
    
    [V, Lambda, whiteSignal, whiteM, dewhiteM] = ica_fuse_calculate_pca(whitesig, numOfIC);
    
    % Eigen vectors
    ICA_Options{length(ICA_Options) + 1} = 'pcaE';
    ICA_Options{length(ICA_Options) + 1} = V;
    
    % Eigen values
    ICA_Options{length(ICA_Options) + 1} = 'pcaD';
    ICA_Options{length(ICA_Options) + 1} = Lambda;
    
    % Whitened signal
    ICA_Options{length(ICA_Options) + 1} = 'whiteSig';
    ICA_Options{length(ICA_Options) + 1} = whiteSignal;
    
    % Whitening matrix
    ICA_Options{length(ICA_Options) + 1} = 'whiteMat';
    ICA_Options{length(ICA_Options) + 1} = whiteM;
    
    % Dewhitening matrix
    ICA_Options{length(ICA_Options) + 1} = 'dewhiteMat';
    ICA_Options{length(ICA_Options) + 1} = dewhiteM;
    
    clear V Lambda whiteSignal whiteM dewhiteM;
    
end

%% For CCICA pass dewhitening matrix and number of subjects of each group as parameters
if (strcmpi(algorithmName, 'ccica'))
    pcaFile = [fusionInfo.run_analysis.prefix, DATA_REDUCTION_FILE, 'comb_', num2str(comb_number), '.mat'];
    load(pcaFile, 'dewhiteM');
    % Dewhitening matrix
    ICA_Options{length(ICA_Options) + 1} = 'dewhiteM';
    ICA_Options{length(ICA_Options) + 1} = dewhiteM;
    
    % No. of subjects in each group
    ICA_Options{length(ICA_Options) + 1} = 'num_subjects';
    ICA_Options{length(ICA_Options) + 1} = numSubjects;
    
    clear dewhiteM;
end

if comb_number ~= 1
    %ICA_Options = {};
    
    if length(ICA_Options) > 1
        tags = cell(1, ceil(length(ICA_Options)/2));
        for nn = 1:length(tags)
            tags{nn} = ICA_Options{2*nn-1};
        end
        matchIndex = strmatch('block', tags, 'exact');
        if ~isempty(matchIndex)
            ICA_Options(2*matchIndex - 1 : 2*matchIndex) = [];
        end
    end
    
end

%icaAlgo = ica_fuse_icaAlgorithm; % available ICA algorithms
% selected ICA algorithm
%algorithmName = deblank(icaAlgo(ica_algorithm, :));

if (strcmpi(analysisType, 'average'))
    
    disp('Using average method ...');
    disp('');
    
    icasig2 = zeros(numRuns, size(whitesig, 1), size(whitesig, 2));
    A2 = zeros(numRuns, size(whitesig, 1), size(whitesig, 1));
    
    % Loop over number of runs
    for nRun = 1:numRuns
        
        fprintf('\n');
        
        % Run ICA
        [icaAlgo, W, A, icasig] = ica_fuse_icaAlgorithm(ica_algorithm, whitesig, ICA_Options);
        
        %A = dewhiteM*pinv(W);
        
        if (nRun > 1),
            rho = zeros(size(icasig, 1), size(icasig, 1));
            for k1 = 1:size(icasig, 1)
                for k2 = 1:size(icasig, 1)
                    rho(k1,k2) = ica_fuse_corr2(flatrow(icasig(k1, :)), flatrow(icasig2(1, k2, :)));
                end
            end
            % rho = (rho1+rho2)/2;
            
            Y = zeros(1, size(icasig, 1));
            I = zeros(1, size(icasig, 1));
            Ys = zeros(1, size(icasig, 1));
            
            for k = 1:size(icasig, 1)
                [Y(k) I(k)] = max(abs(rho(:,k)));
                Ys(k) = rho(I(k),k);%get signed correlation
                rho(I(k),k) = 0;
            end;
            
            %reorder and force to be positively correlated
            icasig = sign(repmat(Ys',1,size(icasig,2))).*icasig(I,:);
            A = sign(repmat(Ys,size(A,1),1)).*A(:,I);
            
        end
        
        % store icasig and A
        if (~iscell(icasig))
            icasig2(nRun, :, :) = icasig;
            A2(nRun, :, :) = A;
        else
            icasig2 = icasig;
            A2 = A;
        end
        
    end
    % end loop over number of runs
    
    if numRuns > 1
        icasig = squeeze(mean(icasig2));
        A = squeeze(mean(A2));
        clear W;
        %W = pinv(A);
    end
    
    clear icasig2;
    clear A2;
    
else
    
    disp('Using ICASSO. Stable ICA run estimates will be used ...');
    disp('');
    
    %%%%% Calculate PCA and Whitening matrix %%%%%
    % PCA
    [V, Lambda] = ica_fuse_v_pca(whitesig, 1, numOfIC, 0, 'transpose', 'yes');
    
    % Whiten matrix
    [w, White, deWhite] = ica_fuse_v_whiten(whitesig, V, Lambda, 'transpose');
    
    clear V Lambda;
    
    sR = ica_fuse_icassoEst('randinit', whitesig, numRuns, 'numOfPC', numOfIC, 'algoIndex', ica_algorithm, ...
        'dewhiteM', deWhite, 'whiteM', White, 'whitesig', w, 'icaOptions', ICA_Options);
    
    clear data w deWhite White;
    
    %%%% Visualization %%%%%%
    
    sR = icassoExp(sR);
    
    %%% Visualization & returning results
    %%% Allow to disable visualization
    if OPEN_DISPLAY_WINDOW
        disp(['Launch Icasso visualization supposing ', num2str(numOfIC), ' estimate-clusters.']);
        disp('Show demixing matrix rows.');
        icassoShow(sR, 'L', numOfIC, 'estimate', 'demixing');
        
        disp(['Launch Icasso visualization supposing ', num2str(numOfIC), ' estimate-clusters.']);
        disp('Show IC source estimates (default), reduce number of lines');
        disp('Collect results.');
        iq = icassoShow(sR, 'L', numOfIC, 'colorlimit', [.8 .9]);
    else
        iq = icassoResult(sR, numOfIC);
    end
    
    minClusterSize = ceil(numRuns*0.8);
    maxClusterSize = numRuns;
    
    [metric_Q, A, W, icasig] = getStableRunEstimates(sR, minClusterSize, maxClusterSize);
    
    icassoResultsFile = [fusionInfo.run_analysis.prefix, ICA_FILE, 'comb_', num2str(comb_number), '_icasso_results.mat'];
    
    icassoResultsFile = fullfile(outputDir, icassoResultsFile);
    ica_fuse_save(icassoResultsFile, 'iq', 'A', 'W', 'icasig', 'sR', 'algorithmName', 'metric_Q');
    clear sR;
    
end

clear whitesig;


if isfield(fusionInfo.run_analysis, 'type_pca') && (strcmpi(fusionInfo.run_analysis.type_pca, 'cca') || strcmpi(fusionInfo.run_analysis.type_pca, 'mcca')) ...
        && (fusionInfo.run_analysis.numFeatures == 2)
    
    corr_modalities = zeros(1, size(icasig, 1));
    if (strcmpi(fusionInfo.run_analysis.type_pca, 'cca') && all(fusionInfo.run_analysis.newDims == fusionInfo.run_analysis.newDims(1)))
        %if (all(fusionInfo.run_analysis.newDims == fusionInfo.run_analysis.newDims(1)))
        
        disp('Sorting components based on the correlation between features ...');
        
        for nC = 1:length(corr_modalities)
            corr_modalities(nC) = ica_fuse_corr2(icasig(nC, 1:fusionInfo.run_analysis.newDims(1)), icasig(nC, fusionInfo.run_analysis.newDims(1) + 1:end));
        end
        
    else
        
        load(pcaFile, 'dewhiteM');
        disp('Sorting components based on the correlation between mixing coefficients ...');
        
        GA1 = dewhiteM{1}*A;
        GA2 = dewhiteM{2}*A;
        
        for nC = 1:length(corr_modalities)
            corr_modalities(nC) = ica_fuse_corr2(GA1(:, nC), GA2(:, nC));
        end
        
        clear GA1 GA2;
        
    end
    
    [dd, corr_inds] = sort(abs(corr_modalities));
    corr_inds = corr_inds(end:-1:1);
    icasig = icasig(corr_inds, :);
    A = A(:, corr_inds);
    corr_modalities = corr_modalities(corr_inds);
    
end

if (~iscell(A))
    W = pinv(A);
end

%%%%%%%% Saving ICA Information %%%%%%%%%%%%%%%
icaFile = [fusionInfo.run_analysis.prefix, ICA_FILE, 'comb_', num2str(comb_number), '.mat'];
fusionInfo.run_analysis.icaFiles(comb_number).name = icaFile;
fusionInfo.run_analysis.icaFiles(comb_number).combinationName = combinationName;
icaFile = fullfile(outputDir, icaFile);

if (exist('corr_modalities', 'var'))
    ica_fuse_save(icaFile, 'W', 'icasig', 'A', 'combinationName', 'corr_modalities');
else
    ica_fuse_save(icaFile, 'W', 'icasig', 'A', 'combinationName');
end
disp(['ICA information for ', combinationName, ' is saved in ', icaFile]);

fprintf('\n');

%% Save fusion file
fusionFile = fusionInfo.run_analysis.fusionFile;
ica_fuse_save(fusionFile, 'fusionInfo');

disp('-----------------------------------------------------------------------------------------------');
disp(['Done ICA/IVA for ',  combinationName]);
disp('-----------------------------------------------------------------------------------------------');

fprintf('\n');


function [data] = flatrow(data)

data = data(:);

function [metric_Q, A, W, icasig] = getStableRunEstimates(sR, minClusterSize, maxClusterSize)
%% Get stable run based on code by Sai Ma. Stable run estimates will be used instead of centrotype
%

% number of runs and ICs
numOfRun = length(sR.W);
numOfIC = size(sR.W{1},1);

% Get the centrotype for each cluster and Iq
index2centrotypes = icassoIdx2Centrotype(sR,'partition', sR.cluster.partition(numOfIC,:));
Iq = icassoStability(sR, numOfIC, 'none');

% Find IC index  within each cluster
partition = sR.cluster.partition(numOfIC, :);
clusterindex = cell(1, numOfIC);
for i = 1:numOfIC
    temp = (partition == i);
    clusterindex{i} = sR.index(temp, :);
    clear temp;
end

% Compute stability metric for each run within each cluster
eachRun = zeros(numOfRun, numOfIC);
qc = 0; % num of qualified clusters
for i = 1:numOfIC
    thisCluster = (clusterindex{i}(:,1))';
    clusterSize = length(clusterindex{i});
    if ((clusterSize >= minClusterSize) && (clusterSize <= maxClusterSize) && (Iq(i)>=0.7))
        qc = qc + 1;
        for k = 1:numOfRun
            thisRun = find(thisCluster == k);
            ICindex = (clusterindex{i}(thisRun,1)-1)*numOfIC + clusterindex{i}(thisRun,2);
            if ~isempty(thisRun)
                eachRun(k,i) = max(sR.cluster.similarity(index2centrotypes(i),ICindex'));
            end
            clear thisRun ICindex;
        end
    end
    clear thisCluster clusterSize;
end

%% Find stable run
metric_Q = sum(eachRun,2)/qc;
[dd, stableRun] = max(metric_Q);

%% Get stable run estimates
W = sR.W{stableRun};
clusters_stablerun = partition((stableRun - 1)*numOfIC + 1 : stableRun*numOfIC);
[dd, inds] = sort(clusters_stablerun);
W = W(inds, :);
A = pinv(W);
icasig = W*sR.signal;
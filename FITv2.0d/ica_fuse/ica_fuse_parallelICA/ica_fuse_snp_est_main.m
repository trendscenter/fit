% Input:
% filepath: directory to save the results
% data_snp: sample*SNP
% ncomp_snp_min: min testing component number, default = 2;
% ncomp_snp_max: max testing component number, default = min([No. of subjects - 1, No. of SNps - 1]);
% estimation method: % to estimate orders based on average consistencies, set type1 = 'column', type2 = [];
% to estimate orders using reference selected from consistency map, set type1 = 'row', type2 = 'blind';
% to estimate orders using reference selected based on a given phenotype, set type1 = 'row', type2 = 'pheno', and provide a variable pheno;
% pheno: phenotype used for reference selection, requried when using method row-pheno
% winsize: sliding window size, default = 10
% testN: number of selected references, default = 3

% Output:
% ncomp_snp_est: the final estimated component number

function ica_fuse_snp_est_main(data_snp, numruns, ncomp_snp_max, verbose)


if (~exist('verbose', 'var'))
    verbose = 'on';
end

if (~exist('data_snp', 'var'))
    files = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.asc', 'title', 'Select SNP files of subjects');
    if (isempty(files))
        error('SNP files are not selected');
    end
    
    if (size(files, 1) < 15)
        error('At least 15 subjects need to be selected in order to use this utility');
    end
    
    data_snp = ica_fuse_loadData(files);
    data_snp = detrend(squeeze(data_snp(:, 2, :)), 0)';
end

if (~exist('numruns', 'var'))
    answer = ica_fuse_inputdlg2({'Enter no. of ICA runs'}, 'No. of ICA runs', 1, {num2str(15)});
    if (isempty(answer) || isempty(answer{1}))
        error('No. of ICA runs is not entered');
    end
    numruns = str2num(answer{1});
end


if (~exist('ncomp_snp_max', 'var'))
    answer = ica_fuse_inputdlg2({'Enter upper limit for components to be estimated'}, 'Max No. of Components', 1, {num2str(min([min(size(data_snp)) - 1, 100]))});
    if (isempty(answer) || isempty(answer{1}))
        error('Upper limit for components to be estimated is not entered');
    end
    ncomp_snp_max = str2num(answer{1});
end


%% ICA on SNP data
type1 = 'column';
ncomp_snp_min = 2;
%ncomp_snp_max = size(data_snp, 1);
[subN, snpN] = size(data_snp);
ncomp_snp_max = min([ncomp_snp_max,subN-1,snpN-1]);
filepath = pwd;
type2='blind';
testN=3;
if (size(data_snp, 1) > 31)
    winsize = 5;
else
    winsize = 3;
end
pheno = [];

clear s_ica;
ica_options = ica_fuse_icaOptions([8,subN],1,'off');
ica_options{10} = 256;
%numruns = 10;
for ncomp_snp = ncomp_snp_min:ncomp_snp_max
    [A_info, W, icasig_info] = getIcassoEstimates(data_snp, ncomp_snp, numruns, verbose);
    s_ica(ncomp_snp).icasig = icasig_info;
    s_ica(ncomp_snp).loading = A_info;
    %     if ~rem(ncomp_snp,5)
    %         save([filepath,'ica.mat'],'s_ica');
    %     end
end
%save([filepath,'ica.mat'],'s_ica');

%% component number estimation

%=== consistency map construction ======================================
[consmap_s, consmap_a, indmap_s, indmap_a] = ica_fuse_snp_est_consmap(filepath, s_ica, ncomp_snp_min, ncomp_snp_max);

[ncomp_snp_est, ncomp_snp_candid] = ica_fuse_snp_est_order(filepath, s_ica, consmap_s, consmap_a, [], [], winsize, testN, ncomp_snp_min, ncomp_snp_max, type1, type2);

disp(['Estimated components is found out to be ', num2str(ncomp_snp_est)]);
fprintf('\n');
fprintf('\n');

%=== plot ====================================
ica_fuse_snp_est_plot(consmap_s,consmap_a,ncomp_snp_est);

function [A, W, icasig] = getIcassoEstimates(data, numOfIC, numRuns, verbose)

[V, Lambda] = ica_fuse_v_pca(data, 1, numOfIC, 0, 'transpose', 'yes');

% Whiten matrix
[w, White, deWhite] = ica_fuse_v_whiten(data, V, Lambda, 'transpose');

clear V Lambda;

sR = ica_fuse_icassoEst('randinit', data, numRuns, 'numOfPC', numOfIC, 'algoIndex', 1, ...
    'dewhiteM', deWhite, 'whiteM', White, 'whitesig', w, 'icaOptions', {'verbose', verbose});

clear data w deWhite White;

%%%% Visualization %%%%%%

sR = icassoExp(sR);

%iq = icassoResult(sR, numOfIC);

[metric_Q, A, W, icasig] = getStableRunEstimates(sR, 2, numRuns);


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


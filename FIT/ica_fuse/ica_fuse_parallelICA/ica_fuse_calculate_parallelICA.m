function varargout = ica_fuse_calculate_parallelICA(data, dewhiteM, whiteM, ...
    numICARuns, type_parallel_ica, modalities, ICA_Options, analysisType, featureData)
% Runs parallel ICA multiple times
%
% Input:
% 1. data - Data is a cell array of length 2
% 2. dewhiteM - de-whitening matrix is a cell array of length 2
% 3. whiteM - whitening matrix is a cell array of length 2
% 4. numICARuns - Number of ICA runs
% 5. type_parallel_ica - Parallel ICA type
% 6. modalities - Modalities used
%
% Output
% 1. aveComp - Average component is of the same dimensions as inpu data
% 2. loadingCoeff - Loading coefficients is a cell array of length 2
% 3. W - weights come out of PICA, it is needed to  back recontrct subject
% specific fmri component
% 4. maxrow - SNP component numbers
% 5. SNPWeight - Indices surpassing SNP_Z_


ica_fuse_defaults;

if (~exist('ICA_Options', 'var'))
    ICA_Options = ica_fuse_paraICAOptions(3, 'off');
end

if (~exist('analysisType', 'var'))
    analysisType = 'average';
end

numComp1 = size(data{1}, 1);
numComp2 = size(data{2}, 1);

% if (length(data) == 3)
%     type_parallel_ica = 'aa';
% end

if ~strcmpi(type_parallel_ica, 'spica')
    chkOpts = strmatch('hoyerindex', ICA_Options(1:2:end), 'exact');
    if (~isempty(chkOpts))
        ICA_Options(2*chkOpts-1:2*chkOpts) = [];
    end
end


if (strcmpi(analysisType, 'average'))
    
    sR1 = [];
    sR2 = [];
    sR3 = [];
    
    disp('Using average method ...');
    disp('');
    
    % Run ICA multiple times
    for run = 1:numICARuns
        
        % Parallel ICA
        if strcmpi(type_parallel_ica, 'aa')
            disp('Using parallel ICA algorithm AA ..');
            if (length(data) == 2)
                [W, sphere, icasig_tmp] = ica_fuse_runica_parallelicaMul_AA(data, 'dewhitem', ...
                    dewhiteM, 'whitem', whiteM, ICA_Options{:});
            else
                % Three way parallel ica
                [W, sphere] = ica_fuse_run_three_way_parallel_ica(data, 'dewhitem', dewhiteM,'whitem', whiteM, ICA_Options{:});
            end
        elseif strcmpi(type_parallel_ica, 'aa-ref')
            disp('Using parallel ICA algorithm AA-ref ..');
            [W, sphere, icasig_tmp] = ica_fuse_runica_picar(data, 'dewhitem', ...
                dewhiteM, 'whitem', whiteM, ICA_Options{:});
        elseif strcmpi(type_parallel_ica, 'TAA')
            disp('Using 3way temporal included parallel ICA algorithm ..');
            [W, sphere] = ica_fuse_run_three_way_parallel_ica_TAA(data, 'dewhitem', dewhiteM,'whitem', whiteM, ICA_Options{:});
        elseif strcmpi(type_parallel_ica,'att')
            disp('Using parallel ICA algorithm ATT')
            [W,sphere,icasig_tem] = ica_fuse_runica_parallelicaMul_ATT(data, 'dewhitem', ...
                dewhiteM, 'whitem', whiteM, ICA_Options{:});
        elseif strcmpi(type_parallel_ica,'at')
            disp('Using parallel ICA algorithm AT')
            [W,sphere,icasig_tem] = ica_fuse_runica_parallelicaMul_AT(data, 'dewhitem', ...
                dewhiteM, 'whitem', whiteM, ICA_Options{:});
            
        elseif strcmpi(type_parallel_ica, 'spica')
            disp('Using parallel ICA algorithm sPICA ..');
            [W, sphere,icasig_tmp,laststep,sparse_sign2,max_corr_corropt_tmp] = ica_fuse_runica_SparseParallelICAMul_AA(data, 'dewhitem', dewhiteM, ...
                'whitem', whiteM, ICA_Options{:});
            
        else
            disp('Using parallel ICA algorithm AS ..');
            [W, sphere, icasig_tmp] = ica_fuse_runica_parallelica_AS(data, 'dewhitem', ...
                dewhiteM, 'whitem', whiteM, ICA_Options{:});
        end
        
        % Calculate the A matrix and icasig for each
        % modality separately
        A = cell(1, length(W));
        %if (~strcmpi(type_parallel_ica, 'taa'))
        for nn = 1:length(W)
            W{nn} = W{nn}*sphere{nn};
            icasig_tmp{nn} = W{nn}*data{nn};
            A{nn} = dewhiteM{nn}*pinv(W{nn});
        end
        %else
        if (strcmpi(type_parallel_ica, 'taa'))
            % TIAA (requires first modality as fmri)
            A{1} = dewhiteM{3}*pinv(W{1});
            A{2} = dewhiteM{1}*pinv(W{2});
            A{3} = dewhiteM{2}*pinv(W{3});
        end
        
        
        gene_modalities = find(strcmpi('gene', modalities));
        
        if  strcmpi(type_parallel_ica, 'spica')
            if (~isempty(gene_modalities))
                for current_gene_index = gene_modalities
                    [A{current_gene_index},T_tmp,C_tmp,alpha_tmp] = AA_est_RegularizedLeastSquare(icasig_tmp{current_gene_index}, featureData(current_gene_index).data);
                end
            end
        end
        
        %end
        %% Divide components by z-scores
        
        % Loop over modalities
        for nM = 1:length(icasig_tmp)
            tempS = icasig_tmp{nM};
            tempA = A{nM};
            % Loop over components
            for nC = 1:size(tempS, 1)
                s = std(tempS(nC, :));
                tempS(nC, :) = tempS(nC, :)./s;
                tempA(:, nC) = tempA(:, nC).*s;
            end
            % End of loop over components
            icasig_tmp{nM} = tempS;
            A{nM} = tempA;
        end
        % End of loop over modalities
        
        %  AA = A{1};
        %  TT = A{2};
        
        %  clear A;
        
        % Record the resultd from different runs
        if run == 1
            
            aveComp = icasig_tmp;
            loadingCoeff = A; %{loadingCoeff1, loadingCoeff2};
            
            %aveComp1 = icasig_tmp{1};
            %aveComp2 = icasig_tmp{2};
            
            % loadingCoeff1 = AA;
            % loadingCoeff2 = TT;
            
        else
            
            for nM = 1:length(aveComp)
                [aveComp{nM}, loadingCoeff{nM}] = sortCompRuns(aveComp{nM}, loadingCoeff{nM}, A{nM}, icasig_tmp{nM});
            end
            %[aveComp2, loadingCoeff2] = sortCompRuns(aveComp2, loadingCoeff2, icasig_tmp{2});
            
            %             % sort the components based on fMRi spacial correlation
            %             temp = (ica_fuse_corr(aveComp1', icasig_tmp{1}'));
            %             rowt = zeros(1, numComp1);
            %             colt = rowt; tcorr = rowt;
            %             for i=1:numComp1
            %                 [t,loc]=max(abs(temp(:)));
            %                 [tempRow, tempCol] = find(abs(temp) == t(1));
            %                 rowt(i) = tempRow(1);
            %                 colt(i) = tempCol(1);
            %                 tcorr(i)=temp(loc(1));
            %                 temp(rowt(i),:)=0;
            %                 temp(:,colt(i))=0;
            %
            %             end
            %             atemp=icasig_tmp{1}(colt,:);
            %             ltemp = AA(:, colt);
            %             ind=find(tcorr<0);
            %             if ~isempty(ind)
            %                 atemp(ind, :)=atemp(ind, :).*(-1);
            %                 ltemp(:, ind) = ltemp(:, ind)*-1;
            %             end
            %
            %             aveComp1 = aveComp1(rowt, :) + atemp;
            %             loadingCoeff1 = loadingCoeff1(:, rowt) + ltemp;
            
            %             % SNP decomposition consistancy check
            %             temp=(ica_fuse_corr(aveComp2',icasig_tmp{2}'));
            %             srowt = zeros(1, numComp2);
            %             scolt = srowt; stcorr = srowt;
            %             for i=1:numComp2
            %                 [t,loc]=max(abs(temp(:)));
            %                 [tempRow, tempCol] = find(abs(temp)==t(1));
            %                 srowt(i) = tempRow(1);
            %                 scolt(i) = tempCol(1);
            %                 stcorr(i)=temp(loc(1));
            %                 temp(srowt(i),:)=0;  temp(:,scolt(i))=0;
            %             end
            %             astemp = icasig_tmp{2}(scolt,:);
            %             lstemp = TT(:, scolt);
            %             ind = find(stcorr<0);
            %             if ~isempty(ind)
            %                 astemp(ind, :) = astemp(ind, :).*(-1);
            %                 lstemp(:, ind) = lstemp(:, ind)*-1;
            %             end
            %             aveComp2=aveComp2(srowt, :) + astemp;
            %             loadingCoeff2 = loadingCoeff2(:, srowt) + lstemp;
            
        end
        % End for number of ICA runs
        
    end
    
    for nRun = 1:length(aveComp)
        aveComp{nM} = aveComp{nM}/run;
        loadingCoeff{nM} = loadingCoeff{nM}/run;
    end
    
    %aveComp1 = aveComp1 / run;
    %aveComp2 = aveComp2 / run;
    
    
    % Average components
    %aveComp = {aveComp1, aveComp2};
    
    %    clear aveComp1 aveComp2;
    
    %     loadingCoeff1 = loadingCoeff1 / run;
    %     loadingCoeff2 = loadingCoeff2 / run;
    
    % Loading coefficients
    %loadingCoeff = {loadingCoeff1, loadingCoeff2};
    
    %    clear loadingCoeff1 loadingCoeff2
    
    
else
    
    
    disp('Using ICASSO. Stable ICA run estimates will be used ...');
    disp('');
    
    
    sR1 = icassoStruct(featureData(1).data);
    sR2 = icassoStruct(featureData(2).data);
    sR3 = [];
    if (length(featureData) == 3)
        sR3 = icassoStruct(featureData(3).data);
    end
    
    k1 = 0; k2 = 0; k3 = 0;
    % Run ICA multiple times
    for i = 1:numICARuns
        
        fprintf('\n\n%s\n\n', ['Randomization round ', num2str(i), '/', num2str(numICARuns)]);
        
        % Parallel ICA
        if strcmpi(type_parallel_ica, 'aa')
            disp('Using parallel ICA algorithm AA ..');
            if (length(data) == 2)
                [W, sphere] = ica_fuse_runica_parallelicaMul_AA(data, 'dewhitem', dewhiteM, 'whitem', whiteM, ICA_Options{:});
            else
                [W, sphere] = ica_fuse_run_three_way_parallel_ica(data, 'dewhitem', dewhiteM, 'whitem', whiteM, ICA_Options{:});
            end
        elseif strcmpi(type_parallel_ica, 'aa-ref')
            disp('Using parallel ICA algorithm AA-ref ..');
            [W, sphere] = ica_fuse_runica_picar(data, 'dewhitem', dewhiteM, 'whitem', whiteM, ICA_Options{:});
            
        elseif strcmpi(type_parallel_ica, 'spica')
            disp('Using parallel ICA algorithm sPICA ..');
            [W, sphere,icasig_tmp,laststep,sparse_sign2,max_corr_corropt_tmp] = ica_fuse_runica_SparseParallelICAMul_AA(data, 'dewhitem', dewhiteM, ...
                'whitem', whiteM, ICA_Options{:});
        else
            disp('Using parallel ICA algorithm AS ..');
            [W, sphere] = ica_fuse_runica_parallelica_AS(data, 'dewhitem', dewhiteM, 'whitem', whiteM, ICA_Options{:});
        end
        
        % Calculate the A matrix and icasig for each
        % modality separately
        A = cell(1, length(W));
        for nn = 1:length(W)
            W{nn} = W{nn}*sphere{nn};
            A{nn} = dewhiteM{nn}*pinv(W{nn});
        end
        
        
        [sR1, k1] = appendMixingMat(sR1, A{1}, k1);
        [sR2, k2] = appendMixingMat(sR2, A{2}, k2);
        if (length(A) == 3)
            [sR3, k3] = appendMixingMat(sR3, A{3}, k3);
        end
        
        
        %         A1 = A{1};
        %         W1 = pinv(A1);
        %
        %         % Store results in icasso
        %         n = size(A1, 2);
        %         if n >0
        %             k1 = k1 + 1;
        %             sR1.index(end+1:end + n, :)= [repmat(k1, n, 1), (1:n)'];
        %             sR1.A{k1} = A1;
        %             sR1.W{k1} = W1;
        %         end
        %
        %
        %         A2 = A{2};
        %         W2 = pinv(A2);
        %         n = size(A2, 2);
        %         if n > 0
        %             k2 = k2 + 1;
        %             sR2.index(end+1:end + n, :)= [repmat(k2, n, 1), (1:n)'];
        %             sR2.A{k2} = A2;
        %             sR2.W{k2} = W2;
        %         end
        %
        
    end
    
    
    minClusterSize = ceil(numICARuns*0.8);
    maxClusterSize = numICARuns;
    
    sR1.mode = 'randinit';
    sR1.dewhiteningMatrix = dewhiteM{1};
    sR1.whiteningMatrix = whiteM{1};
    
    sR2.mode = 'randinit';
    sR2.dewhiteningMatrix = dewhiteM{2};
    sR2.whiteningMatrix = whiteM{2};
    
    
    %% ICASSO on modality 1
    sR1 = icassoExp(sR1);
    %iq1 = icassoResult(sR1, size(dewhiteM{1}, 2));
    %iq1 = icassoShow(sR1, 'L', size(dewhiteM{1}, 2), 'colorlimit', [.8 .9]);
    [metric_Q, loadingCoeff{1}, W1, aveComp{1}] = getStableRunEstimates(sR1, minClusterSize, maxClusterSize);
    
    
    %% ICASSO on modality 2
    sR2 = icassoExp(sR2);
    %iq2 = icassoResult(sR2, size(dewhiteM{2}, 2));
    %iq2 = icassoShow(sR2, 'L', size(dewhiteM{2}, 2), 'colorlimit', [.8 .9]);
    [metric_Q, loadingCoeff{2}, W2, aveComp{2}] = getStableRunEstimates(sR2, minClusterSize, maxClusterSize);
    
    if (length(A) == 3)
        sR3.mode = 'randinit';
        sR3.dewhiteningMatrix = dewhiteM{3};
        sR3.whiteningMatrix = whiteM{3};
        %% ICASSO on modality 3
        sR3 = icassoExp(sR3);
        %iq2 = icassoResult(sR2, size(dewhiteM{2}, 2));
        %iq2 = icassoShow(sR2, 'L', size(dewhiteM{2}, 2), 'colorlimit', [.8 .9]);
        [metric_Q, loadingCoeff{3}, W3, aveComp{3}] = getStableRunEstimates(sR3, minClusterSize, maxClusterSize);
    end
    
    
    gene_modalities = find(strcmpi('gene', modalities));
    if  strcmpi(type_parallel_ica, 'spica')
        if (~isempty(gene_modalities))
            for current_gene_index = gene_modalities
                [loadingCoeff{current_gene_index},T_tmp,C_tmp,alpha_tmp] = AA_est_RegularizedLeastSquare(aveComp{current_gene_index}, featureData(current_gene_index).data);
            end
        end
    end
    
    for nM = 1:length(loadingCoeff)
        tempS = aveComp{nM};
        tempA = loadingCoeff{nM};
        % Loop over components
        for nC = 1:size(tempS, 1)
            s = std(tempS(nC, :));
            tempS(nC, :) = tempS(nC, :)./s;
            tempA(:, nC) = tempA(:, nC).*s;
        end
        % End of loop over components
        aveComp{nM} = tempS;
        loadingCoeff{nM} = tempA;
    end
    % End of loop over modalities
    
    
end

maxcol = []; maxrow = [];
maxcorr = [];
try
    if (length(data) == 2)
        maxcorr = zeros(1, numComp1);
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
    end
catch
end

if ~(strcmpi(type_parallel_ica, 'at') || strcmpi(type_parallel_ica, 'taa'))
    varargout = {aveComp, loadingCoeff, maxcorr, maxrow, sR1, sR2, sR3};
else
    varargout = {aveComp, loadingCoeff, W};
end


% Average correlation
%avecorr = avecorr / run;


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


function [aveComp1, loadingCoeff1] = sortCompRuns(aveComp1, loadingCoeff1, AA, icasig_tmp)
%% Sort components between the runs

numComp1 = size(aveComp1, 1);
% sort the components based on fMRi spacial correlation
temp = (ica_fuse_corr(aveComp1', icasig_tmp'));
rowt = zeros(1, numComp1);
colt = rowt; tcorr = rowt;
for i=1:numComp1
    [t,loc]=max(abs(temp(:)));
    [tempRow, tempCol] = find(abs(temp) == t(1));
    rowt(i) = tempRow(1);
    colt(i) = tempCol(1);
    tcorr(i)=temp(loc(1));
    temp(rowt(i),:)=0;
    temp(:,colt(i))=0;
    
end
atemp=icasig_tmp(colt,:);
ltemp = AA(:, colt);
ind=find(tcorr<0);
if ~isempty(ind)
    atemp(ind, :)=atemp(ind, :).*(-1);
    ltemp(:, ind) = ltemp(:, ind)*-1;
end

aveComp1 = aveComp1(rowt, :) + atemp;
loadingCoeff1 = loadingCoeff1(:, rowt) + ltemp;

function [sR1, k1] = appendMixingMat(sR1, A1, k1)
%% Append mixing matrix to results structure
%

W1 = pinv(A1);

% Store results in icasso
n = size(A1, 2);
if n >0
    k1 = k1 + 1;
    sR1.index(end+1:end + n, :) = [repmat(k1, n, 1), (1:n)'];
    sR1.A{k1} = A1;
    sR1.W{k1} = W1;
end


function [A_alpha_n_T,T,C,alpha,scale_IC] = AA_est_RegularizedLeastSquare(icasig,X_noise,noise_thr,signal_thr)
%%%Reconstruct the loading matrix using Tychonov-regularized least squares
if ~exist('noise_thr','var') %%if user does not input the noise threshold
    noise_thr = 1;
end

if ~exist('signal_thr','var') %%if user does not input the signal threshold
    signal_thr = 2.5;
end

noise_region_grudtrth = find(abs(zscore(icasig(1,:)))<noise_thr);
signal_region_potential = [];
noise_region_potential = [];
for i = 1:size(icasig,1)
    icasig_z(i,:) = zscore(icasig(i,:));
    signal_region_potential = [signal_region_potential,find(abs(icasig_z(i,:))>signal_thr)];
    noise_region_potential_tmp = intersect(noise_region_grudtrth,find(abs(icasig_z(i,:))<noise_thr));
    noise_region_potential = [noise_region_potential,noise_region_potential_tmp];
end
signal_region_potential = unique(signal_region_potential);
noise_region_potential = unique(noise_region_potential);
scale_IC = max(icasig_z(:))-min(icasig_z(:));
icasig_z_n = icasig_z./(max(icasig_z(:))-min(icasig_z(:)));
%%%%%SVD and PCA and weighted least square are sensentive to the scale, but the scale of components
%%%%%from ICA is arbitratry, so need to normalize the scale of S, then do SVD
[U,S,V] = svd(icasig_z_n,'econ');%%SVD
C = 1/((min(diag(S)))^3);  %%%get the minimum singular value and use it to compute C
%%%%using weighted least square to compute A matrix, one subject at once
for i = 1:size(X_noise,1)
    A_alpha0_T(i,:) = inv(icasig_z_n*icasig_z_n')*icasig_z_n*X_noise(i,:)';
    alpha(i) = (norm(X_noise(i,noise_region_potential),2)/(2*C*norm(X_noise(i,signal_region_potential),2)))^(2/3);%%compute the optimal alpha for each subject
    A_alpha_n_T(i,:) = inv(icasig_z_n*icasig_z_n'+alpha(i).*eye(size(icasig_z_n,1)))*icasig_z_n*X_noise(i,:)';
end
%%%test if the error bounds condition can be satisfied by using the suggested alpha.
for i = 1:size(X_noise,1)
    diff_t1 = norm((A_alpha0_T(i,:)-A_alpha_n_T(i,:)),2);
    cmp_val = (2*C*norm(X_noise(i,signal_region_potential),2)*(norm(X_noise(i,noise_region_potential),2)^2))^1/3;
    T(i) = (diff_t1<=cmp_val);
end

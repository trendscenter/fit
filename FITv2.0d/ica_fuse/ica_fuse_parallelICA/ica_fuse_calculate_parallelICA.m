function [aveComp, loadingCoeff, maxcorr, maxrow, sR1, sR2, sR3] = ica_fuse_calculate_parallelICA(data, dewhiteM, whiteM, ...
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
% 3. maxrow - SNP component numbers
% 4. SNPWeight - Indices surpassing SNP_Z_


ica_fuse_defaults;

if (~exist('ICA_Options', 'var'))
    ICA_Options = ica_fuse_paraICAOptions(3, 'off');
end

if (~exist('analysisType', 'var'))
    analysisType = 'average';
end

numComp1 = size(data{1}, 1);
numComp2 = size(data{2}, 1);

if (length(data) == 3)
    type_parallel_ica = 'aa';
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
        else
            disp('Using parallel ICA algorithm AS ..');
            [W, sphere, icasig_tmp] = ica_fuse_runica_parallelica_AS(data, 'dewhitem', ...
                dewhiteM, 'whitem', whiteM, ICA_Options{:});
        end
        
        % Calculate the A matrix and icasig for each
        % modality separately
        A = cell(1, length(W));
        for nn = 1:length(W)
            W{nn} = W{nn}*sphere{nn};
            icasig_tmp{nn} = W{nn}*data{nn};
            A{nn} = dewhiteM{nn}*pinv(W{nn});
        end
        
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
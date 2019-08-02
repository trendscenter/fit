function [joint_ica_sig, dewhiteM, whiteM] = ica_fuse_mccaN(data, dims, numPC, numComp, featureNames)
%% Do N way MCCA
%
%

ica_fuse_defaults;
global MCCA_COST_FUNC;

mcca_costfun = MCCA_COST_FUNC;
if (isempty(mcca_costfun))
    mcca_costfun = 'ssqcor';
end

if (~exist('featureNames', 'var'))
    featureNames = cellstr(num2str(1:length(dims))');
end

pcasig = cell(1, length(dims));

endT = 0;
for n = 1:length(dims)
    startT = endT + 1;
    endT = endT + dims(n);
    pcasig{n} = pcawhiten(data(:, startT:endT)', numPC(n), featureNames{n});
end

W = mcca(pcasig, numComp, mcca_costfun);
dewhiteM = cell(1, length(W));
whiteM = dewhiteM;

endT = 0;
for n = 1:length(W)
    startT = endT + 1;
    endT = endT + dims(n);
    dewhiteM{n} = (W{n}*pcasig{n})';
    pcasig{n} = dewhiteM{n}'*data(:, startT:endT);
    whiteM{n} = pinv(dewhiteM{n});
end

joint_ica_sig = cat(2, pcasig{:});

function [pcasig, dw, L] = pcawhiten(data, M, featureName)
%% Do whitening
%

disp(['Doing SVD on feature ', featureName]);
[U1, S1, V1] = svd(data, 0);
L1 = (inv(S1(1:M, 1:M)));
L = S1(1:M, 1:M).^2;
Y1 = L1*U1(:,1:M)'*data;
pcasig = Y1;
dw = L1*U1(:,1:M)';
temp = diag(S1).^2;
reminfo = sum(temp(1:M))/sum(temp);
disp([num2str(reminfo*100), '% of (non-zero) eigenvalues retained in feature ', featureName]);
fprintf('\n');

function W=mcca(X,numOfCV,costType)
% W = mcca(X,numOfCV,costType)
%
% Multiset Canonical Correlation Analysis using specified cost function.
% Cost functions implemented are {'genvar'}, 'maxvar', 'minvar', 'ssqcor'
%
% This version is from an original implementation by Yiou Li in Dec 2007
% based on [Kettenring, Biometrica 1971].
%
% References:
% [1] J. R. Kettenring, "Canonical analysis of several sets of variables Biometrika," 1971, 58, 433-451
% [2] Y.-O. Li, T. Adali, W. Wang, V. D. Calhoun, "Joint Blind Source Separation by Multiset Canonical Correlation Analysis," IEEE Trans. Signal Process., 2009, 57, 3918-3929


if nargin==0
    help mcca
    demo_mcca
    return
end

if nargin==1
    numOfCV=size(X,1);
end
if nargin==2
    costType='genvar';
end

switch lower(costType)
    case 'genvar'
        W = mcca_genvar(X,numOfCV);
    case 'maxvar'
        W = mcca_maxminvar(X,numOfCV,1);
    case 'minvar'
        W = mcca_maxminvar(X,numOfCV,0);
    case 'ssqcor'
        W = mcca_ssqcor(X,numOfCV);
    otherwise
        error('Unknown MCCA cost type option')
end


function B = mcca_genvar(y, numOfCV)
% W = mcca_genvar(X,numOfCV)
%
% Multiset Canonical Correlation Analysis using GENVAR
%
% This version is from an original implementation by Yiou Li in Dec 2007 based on
% [Kettenring, Biometrica 1971].

[N, T] = size(y{1});
K = length(y);

%[N,T,K]=size(X);
%y=cell(K,1);
%V=cell(K,1);
% for k=1:K
%    %[y{k},V{k}]=whiten(squeeze(X(:,:,k)));
%    y{k} = squeeze(X(:,:,k));
% end

R=cell(K, K);
for k1=1:K
    for k2=1:K
        R{k1,k2}=y{k1}*y{k2}'/T;
    end
end

numMaxIter = 1000;
M = size(R,1);
p = size(R{1,1},1);
eps = 0.0001;

clear y;
B = cell(M,1);

for i = 1:M
    B0{i} = eye(numOfCV,p);
end

for s = 1:numOfCV
    
    %% Iterations to solve the s-th stage canonical vectors for 1--M
    %% datasets
    theta_old = 0; % zeros(M,1);
    theta = 0; % zeros(M,1);
    for n = 1:numMaxIter
        
        if n == 1
            %% Initialize B{1--M}(s,:) by B0;
            for j = 1:M
                % B{j}(s,:) = B0{j}(s,:);
                % B{j}(s,:) = randn(1,p)/norm(randn(1,p));
                B{j}(s,:) = ones(1,p)/norm(ones(1,p));
            end
            %% Calculate the cost funtion at the initial step
            for j = 1:M
                for k = 1:M
                    R_hat(j,k) = B0{j}(s,:)*R{j,k}*B0{k}(s,:)';
                end
            end
            theta_0 = trace(R_hat*R_hat');
        end
        
        %% Solve the current canonical vector for the j-th dataset
        for j = 1:M
            
            %             %% Calculate the cost function
            %             jtheta_old(j) = 0;
            %             for k = 1:M
            %                 jtheta_old(j) = jtheta_old(j) + B{k}(s,:)*R{k,j}*B{j}(s,:)';
            %             end
            
            %% Calculate the terms for updating jbn
            jC = B{j}(1:s-1,:)';
            if s ~= 1
                jA = eye(p) - jC*jC'; % *inv(jC'*jC)
            else
                jA = eye(p);
            end
            
            jPhi = R_hat(([1:M]~=j),([1:M]~=j));
            for k = 1:M
                jN(:,k) = R{j,k}*B{k}(s,:)';
            end
            jN = jN(:,([1:M]~=j));
            jQ = jN*inv(jPhi)*jN';
            %% update jbn
            [Ev Dv] = eig(jA*jQ);
            DD = diag(Dv);
            [maxv maxi] = max(DD);
            B{j}(s,:) = Ev(:,maxi)';
            tmp(j) = det(jPhi)*(2-Dv(maxi,maxi)); % should = jtheta(j)
            
            %             %% Calculate the cost function
            %             jtheta(j) = 0;
            %             for k = 1:M
            %                 jtheta(j) = jtheta(j) + B{k}(s,:)*R{k,j}*B{j}(s,:)';
            %             end
            %             chec(j) = tmp(j) - jtheta(j);
            %             delta(j) = jtheta(j) - jtheta_old(j);
        end
        
        %% Calculate the cost funtion at the current step
        for j = 1:M
            for k = 1:M
                R_hat(j,k) = B{j}(s,:)*R{j,k}*B{k}(s,:)';
            end
        end
        theta(n) = det(R_hat);
        
        %% Check termination condition
        if sum(abs(theta(n) - theta_old)) < eps | n == numMaxIter
            break;
        end
        theta_0;
        theta_old = theta(n);
    end
    
    fprintf('\n Component #%d is estimated, in %d iterations',s, n);
    
end
fprintf('\n')

%% Sorting is not necessary -- just prefer to order sources according to degree of correlation
detSCV=zeros(1,numOfCV);
for n=1:numOfCV
    % Efficient version of Sigma_n=Yn*Yn'/T;
    Sigma_n=eye(K); %
    for k1=1:K
        for k2=(1+k1):K
            Sigma_n(k1,k2)=B{k1}(n,:)*R{k1,k2}*B{k2}(n,:)';
            Sigma_n(k2,k1)=Sigma_n(k1,k2)';
        end % k2
    end %k3
    detSCV(n)=det(Sigma_n);
end

[sortedDetSCV,isort]=sort(detSCV);

% Sort W for unwhitened data.
for k=1:K
    B{k}=B{k}(isort,:);
end


function B = mcca_maxminvar(y, numOfCV, cost_type)
% W=mcca_maxminvar(X,numOfCV,cost_type)
%
% Multiset Canonical Correlation Analysis using MAXVAR or MINVAR
%
% This version is from an original implementation by Yiou Li in Dec 2007 based on
% [Kettenring, Biometrica 1971].


MAXVAR=1;
MINVAR=0;
if nargin<3
    cost_type=MAXVAR;
end

%% Input: y: cell array containing datasets in M_i(number of mixture)xN(number of sample) form

%%        numOfCV: the number of canonical variates to be estimated, should be less than min{M_1, M_2,...}

%% Output: B: cell array containing the canonical transformations for each dataset



N = size(y{1},2);

M = length(y);



%% Formulate the covariance matrix of the concatenated PC's

for i = 1:M;
    
    numPCv(i) = size(y{i},1);
    
    if i == 1
        
        yc = y{i};
        
    else
        
        yc = vertcat(yc,y{i});
        
    end
    
end

R = (yc*yc')/N;

p = sum(numPCv);



%% Performing multigroup CCA using MAXVAR and MINVAR criteria [Kettenring 1971]

s = numOfCV;

% s = min(numPCv); % number of stages



B = cell(M,1);

for i = 1:s  % stage index
    %% Formulate the constraint matrix Dc for the current stage
    
    Dc = cell(M,M);
    
    for j = 1:M  % j: index of dataset
        if i ~= 1
            for k = 1:i-1  % k: index of CVs  ~ Constraint 1* in Biometrika paper
                if k == 1
                    Dc{j,j} = b{j}{k};
                else
                    Dc{j,j} = horzcat(Dc{j,j},b{j}{k});
                end
            end
        end
    end
    
    %% Pad zero-submatrices into Dc
    for j1 = 1:M
        for j2 = 1:M
            if j2 ~= j1
                Dc{j1,j2} = zeros(size(Dc{j1,j1}));
            end
        end
    end
    
    Dc = cell2mat(Dc);
    
    %% Calculate the canonical coefficient vectors
    
    if i == 1
        A = eye(p);
    else
        A = eye(p) - Dc*inv(Dc'*Dc)*Dc';
    end
    
    if cost_type==MAXVAR
        [ev, lamb] = eigs(A*R*A,1);
    else % MINVAR
        [ev, lamb] = eigs(A*R*A,M);
        ev=ev(:,end);
    end
    
    %% Partition the canonical coefficient vectors for each dataset
    
    for j = 1:M  % j: index of dataset
        
        idx_eigvect = 1;
        
        if j == 1
            tmp = ev(1:numPCv(1),idx_eigvect);
        else
            tmp = ev(sum(numPCv(1:j-1))+1:sum(numPCv(1:j)),idx_eigvect);
        end
        
        tmp = tmp/sqrt(tmp'*tmp);
        
        b{j}{i} = tmp;
        
        if i == 1  % i: index of CVs
            B{j} = tmp';
        else
            B{j} = vertcat(B{j},tmp');
        end
    end
    
    clear Dc;
    
end



function B = mcca_ssqcor(y, numOfCV, B0)
% W = mcca_ssqcor(X,numOfCV,B0)
%
% Multiset Canonical Correlation Analysis using SSQCOR
%
% This version is from an original implementation by Yiou Li in Mar 2009 based on
% [Kettenring, Biometrica 1971].


%% This code implement the M-CCA algorithm based on the SSQCOR cost
%% Reference:
%% J. R. Kettenring, "Canonical analysis of several sets of variables,"
%% Biometrika, vol. 58, pp. 433-51, 1971.

%% Input:
%% y: M by 1 cell array containing the *prewhitened* group datasets
%% numOfCV: Number conanical vectors to be estimated
%% B0: M by 1 cell array containing the initial guess of the demixing
%% matrices for the group dataset: default is identity matrix
%% Output:
%% B: M by 1 cell array containing the estimated demixing matrices
%% theta_opt: Vector containing cost function values at the optimal
%% solutions

%% Yiou (Leo) Li Mar. 2009

[N, T] = size(y{1});
K = length(y);

% [N,T,K]=size(X);
% y=cell(K,1);
% V=cell(K,1);
% for k=1:K
%    [y{k},V{k}]=whiten(squeeze(X(:,:,k)));
% end

numMaxIter = 1000;
M = length(y);
[p, T] = size(y{1});
eps = 0.0001;

%% Calculate all pairwise correlation matrices
R=cell(M,M);
for i = 1:M
    for j = i:M
        R{i,j} = y{i}*y{j}'/T;
        R{j,i} = R{i,j}';
    end
    R{i,i} = eye(p);
end

%% Obtain a prilimary estimate by MAXVAR algorithm
if nargin < 3
    %     B0 = maxminvar_cca(y, numOfCV);
    B0 = cell(M,1);
    for i = 1:M
        B0{i} = eye(numOfCV,p);
    end
end
clear y;
B = cell(M,1);

for s = 1:numOfCV
    
    %% Iterations to solve the s-th stage canonical vectors for 1--M
    %% datasets
    theta_old = 0; % zeros(M,1);
    theta = 0; % zeros(M,1);
    for n = 1:numMaxIter
        
        if n == 1
            %% Initialize B{1--M}(s,:) by B0;
            for j = 1:M
                B{j}(s,:) = B0{j}(s,:)/norm(B0{j}(s,:));  %% Use normalized B0
            end
            %% Calculate the cost funtion at the initial step
            for j = 1:M
                for k = 1:M
                    R_hat(j,k) = B0{j}(s,:)*R{j,k}*B0{k}(s,:)';
                end
            end
            theta_0 = trace(R_hat*R_hat');
        end
        
        %% Solve the current canonical vector for the j-th dataset
        for j = 1:M
            
            %             %% Calculate the cost function
            %             jtheta_old(j) = 0;
            %             for k = 1:M
            %                 jtheta_old(j) = jtheta_old(j) + B{k}(s,:)*R{k,j}*B{j}(s,:)';
            %             end
            
            %% Calculate the terms for updating jbn
            jC = B{j}(1:s-1,:)';
            if s ~= 1
                jA = eye(p) - jC*jC'; % *inv(jC'*jC)
            else
                jA = eye(p);
            end
            jP = 0;
            for k = 1:M
                if k ~= j
                    jP = jP + R{j,k}*B{k}(s,:)'*B{k}(s,:)*R{k,j};
                end
            end
            %% update jbn
            [Ev Dv] = eig(jA*jP);
            DD = diag(Dv);
            [maxv maxi] = max(DD);
            B{j}(s,:) = Ev(:,maxi)';
            tmp(j) = Dv(maxi,maxi) + 1; % should = jtheta(j)
            
            %             %% Calculate the cost function
            %             jtheta(j) = 0;
            %             for k = 1:M
            %                 jtheta(j) = jtheta(j) + B{k}(s,:)*R{k,j}*B{j}(s,:)';
            %             end
            %             chec(j) = tmp(j) - jtheta(j);
            %             delta(j) = jtheta(j) - jtheta_old(j);
        end
        
        %% Calculate the cost funtion at the current step
        for j = 1:M
            for k = 1:M
                R_hat(j,k) = B{j}(s,:)*R{j,k}*B{k}(s,:)';
            end
        end
        theta(n) = trace(R_hat*R_hat');
        
        %% Check termination condition
        if sum(abs(theta(n) - theta_old)) < eps | n == numMaxIter
            theta_opt(s) = theta(n);
            break;
        end
        theta_0;
        theta_old = theta(n);
    end
    
    fprintf('\n Component #%d is estimated, in %d iterations',s, n);
    
end
fprintf('\n ')

%% Sorting is not necessary -- just prefer to order sources according to degree of correlation
detSCV=zeros(1, numOfCV);
for n=1:numOfCV
    % Efficient version of Sigma_n=Yn*Yn'/T;
    Sigma_n=eye(K); %
    for k1=1:K
        for k2=(1+k1):K
            Sigma_n(k1,k2)=B{k1}(n,:)*R{k1,k2}*B{k2}(n,:)';
            Sigma_n(k2,k1)=Sigma_n(k1,k2)';
        end % k2
    end %k3
    detSCV(n)=det(Sigma_n);
end

[sortedDetSCV,isort]=sort(detSCV);

% Sort W for unwhitened data.
for k=1:K
    B{k}=B{k}(isort,:);
end






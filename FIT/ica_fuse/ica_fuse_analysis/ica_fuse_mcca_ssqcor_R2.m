function [W,costfunction,cf11,cf22] = ica_fuse_mcca_ssqcor_R2(X,numOfCV,ref,lamda,B0)
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

[N,T,K]=size(X);
one = ones(1,N);
y=cell(K,1);
V=cell(K,1);

for k=1:K
   [y{k},V{k}]=whiten(squeeze(X(:,:,k)));
end
Ref = cell(K,1);

if (length(ref) == numel(ref(:)))
    
    ref_matrix = repmat(ref,1,N);
    
    for k = 1:K
        Ref{k} = corr(y{k}',ref_matrix,'rows','complete');
    end
    
else
    ref_matrix = detrend(ref, 0);
    for k = 1:K
        tmpyy = y{k}';
        tmp_ref = zeros(size(tmpyy, 2), size(tmpyy, 2));
        for nYY = 1:size(tmpyy, 2)
            [~, tmp_ref(nYY, :)] = ica_fuse_regress(detrend(tmpyy(:, nYY), 0), ref_matrix);
        end
        Ref{k} = tmp_ref;
    end
end


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
if nargin < 5
%     B0 = maxminvar_cca(y, numOfCV);
   B0 = cell(M,1);
    for i = 1:M
        B0{i} = eye(numOfCV,p);
    end
end
clear y;
B = cell(M,1);
costfunction = cell(numOfCV,1);
cf11 = cell(numOfCV,1);
cf22 = cell(numOfCV,1);
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
                f2(j) = B{j} (s,:)*Ref{j}*((1/5)*one');
            end
            theta_0 = trace(R_hat*R_hat')+50*lamda*sum(f2.^2);
%             theta_0 = trace(R_hat*R_hat')+50*lamda*trace(f2*f2');
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
            MATRIX = jA*jP+jA*(lamda*Ref{j}*one'*one*Ref{j}');
            [Ev Dv] = eig(MATRIX);
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
            f2(j) = B{j} (s,:)*Ref{j}*((1/5)*one');
        end
        theta(n) = trace(R_hat*R_hat')+50*50*lamda*sum(f2.^2);
        cf1(n) = trace(R_hat*R_hat');
        cf2(n) = 50*sum(f2.^2);
        %% Check termination condition
        if sum(abs(theta(n) - theta_old)) < eps | n == numMaxIter
            theta_opt(s) = theta(n);
            break;
        end
        theta_0;
        theta_old = theta(n);
    end
    costfunction{s,1} = theta;
    cf11{s,1} = cf1;
    cf22{s,1} = cf2;
    fprintf('\n Component #%d is estimated, in %d iterations',s, n);
   
end
fprintf('\n ')

%% Sorting is not necessary -- just prefer to order sources according to degree of correlation
% detSCV=zeros(1,N);
% for n=1:N
%    % Efficient version of Sigma_n=Yn*Yn'/T;
%    Sigma_n=eye(K); %
%    for k1=1:K
%       for k2=(1+k1):K
%          Sigma_n(k1,k2)=B{k1}(n,:)*R{k1,k2}*B{k2}(n,:)';
%          Sigma_n(k2,k1)=Sigma_n(k1,k2)';
%       end % k2
%    end %k3
%    detSCV(n)=det(Sigma_n);
% end
% 
% [sortedDetSCV,isort]=sort(detSCV);
% 
% % Sort W for unwhitened data.
% for k=1:K
%    B{k}=B{k}(isort,:);   
% end

W=zeros(N,N,K);
for k=1:K
   W(:,:,k)=B{k}*V{k};
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,V,U]=whiten(x)
% [z,V,U]=whiten(x)
% 
% Whitens the data vector so that E{zz'}=I, where z=V*x.

[N,T]=size(x);

% Step 1. Center the data.
x=x-repmat(mean(x,2),1,T);

% Step 2. Form MLE of data covariance.
covar=x*x'/T;

% Step 3. Eigen decomposition of covariance.
[eigvec, eigval] = eig (covar);

% Step 4. Forming whitening transformation.
V=inv (sqrt (eigval)) * eigvec';
U=eigvec * sqrt(eigval);

% Step 5. Form whitened data
z=V*x;

return
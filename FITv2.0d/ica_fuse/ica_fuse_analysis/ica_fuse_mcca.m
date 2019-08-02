function [C1, C2, A1, A2, Lambda] = ica_fuse_mcca(X1, X2, ncomp1, ncomp2, numComp, featureNames)
%% Do MCCA
%
% Inputs:
% 1. X1 - Data of modality 1
% 2. X2 - Data of modality 2
% 3. ncomp1 - No. of components of modality 1
% 4. ncomp2 - No. of components of modality 2
% 5. numComp - No. of components
%
% Outputs:
% 1. C1 - Reduced data of modality 1
% 2. C1 - Reduced data of modality 2
% 3. A1 - Canonical coefficents of modality 1
% 4. A2 - Canonical coefficients of modality 2
% 5. Lambda - Canonical correlations
%

if (ncomp1 > size(X1, 1)) || (ncomp2 > size(X2, 1))
    error('number of observed signals is less than the component number!!');
end

if (~exist('numComp', 'var'))
    numComp = min([ncomp1, ncomp2]);
end

if (~exist('featureNames', 'var'))
    featureNames = {'1', '2'};
end

%% Dimension reduction of X1
[U1, S1, V1] = svd(X1', 0);
Y1 = (S1(1:ncomp1, 1:ncomp1) \ U1(:, 1:ncomp1)')*X1';
temp = diag(S1).^2;
reminfo = sum(temp(1:ncomp1))/sum(temp);
fprintf('\n');
disp(['Doing SVD on feature ', featureNames{1}]);
disp([num2str(reminfo*100),'% of (non-zero) eigenvalues retained']);
fprintf('\n');

%% Dimension reduction of X1
[U2, S2, V2] = svd(X2',0);
Y2 = (S2(1:ncomp2, 1:ncomp2) \ U2(:, 1:ncomp2)')*X2';
temp = diag(S2).^2;
reminfo = sum(temp(1:ncomp2))/sum(temp);
disp(['Doing SVD on feature ', featureNames{2}]);
disp([num2str(reminfo*100),'% of (non-zero) eigenvalues retained.']);
fprintf('\n');

%% CCA
[C1, C2, Lambda, A1,A2] = ica_fuse_cca(Y1', Y2', numComp);

C1 = pinv(U1(:,1:ncomp1)')*S1(1:ncomp1, 1:ncomp1)*Y1*pinv(A1');
C2 = pinv(U2(:,1:ncomp2)')*S2(1:ncomp2, 1:ncomp2)*Y2*pinv(A2');
C1 = C1';
C2 = C2';


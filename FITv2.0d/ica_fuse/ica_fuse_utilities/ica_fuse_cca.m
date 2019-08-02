function varargout = ica_fuse_cca(X, Y, numComp)
%% Compute canonical coefficients and variates
%
% Inputs:
% 1. X - Data of size N by m1
% 2. Y - Data of size N by m2
% 3. numComp - No. of components
%
% Outputs:
% Variable no. of outputs
% 1. V1 - Canonical coefficients of X
% 2. V2 - Canonical coefficients of Y
% 3. r - Canonical Correlations
% 4. U1 - Canonical variates of X
% 5. U2 - Canonical variates of Y
%

%% Remove mean
X = detrend(X, 0);
Y = detrend(Y, 0);

n = max([1, size(X, 1) - 1]);

if (~exist('numComp', 'var'))
    numComp = min([size(X, 2), size(Y, 2)]);
end

%% Compute covariances
covxx = (X'*X)/n;
covyy = (Y'*Y)/n;
covxy = (X'*Y)/n;
covyx = (Y'*X)/n;

%% Canonical coefficients
[V1, r] = eig(inv(covxx)*covxy*inv(covyy)*covyx);
r = diag(r);
[dd, inds] = sort(r);
inds = inds(end:-1:1);
inds = inds(1:numComp);
r = r(inds);
V1 = V1(:, inds);
r = sqrt(r);

V2 = inv(covyy)*covyx*V1;
V2 = V2.*repmat(1./sqrt(sum(V2.^2)), size(V2, 1), 1);

varargout = {V1, V2, r};

%% Canonical variates
if (nargout > 3)
    U1 = X*V1;
    varargout{end + 1} = U1;
end

if (nargout > 4)
    U2 = Y*V2;
    varargout{end + 1} = U2;
end






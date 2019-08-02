function A = ica_fuse_cov(X, removeMean)
%% Compute covariance matrix. If out of memory error occurs, covariance matrix is
% computed with least memory requirements.
%
% Inputs:
% 1. X - Data of size m by n
% 2. removeMean - Remove mean. Options are 0 and 1.
%
% Outputs:
%
% 1. A - Covariance matrix
%

if (~exist('X', 'var'))
    error('Data variable is not passed');
end

%% Use defaults
if (~exist('removeMean', 'var'))
    removeMean = 1;
end

[m, n] = size(X);

%% Remove mean
if (removeMean)
    X = ica_fuse_remove_mean(X, 1);
end

try
    A = X'*X;
catch
    %% Use less memory for computing covariance matrix
    A = zeros(n, n);

    %% Loop over cols
    for i = 1:n
        currentVector = X(:, i);
        inds = (i:n);
        % Loop over nInd
        for nInd = inds
            temp = currentVector(:)'*X(:, nInd);
            %% Utilise the symmetric nature of the covariance matrix
            A(i, nInd) = temp;
            A(nInd, i) = temp;
        end
        % End loop over nInd
    end
    %% End loop over cols

end

%% Return the final result
if (m > 1)
    A = A ./ (m - 1);
else
    A = A ./ m;
end

function [V, Lambda, V_orig, Lambda_orig] = ica_fuse_pca_reference(data, firstEig, lastEig, reference, saveEig, ...
    doTranspose, removeMean);

% compute covariance matrix by removing the mean and then calculate the
% eigen values and eigen vectors
% THen,projected into the reference direction, recalculate the eigenvalues and
% eigen vectors

% Input: 1. data - 2D array, 1 dimension - voxels, 2 dimension - time points
% 2., 3. and 4.,  firstEig, lastEig and reference, saveEig - scalars
% 5. doTranspose - 'transpose' - computes transpose

%
% firstEig_t2 = firstEig;lastEig_t2 = lastEig;saveEig_t2=saveEig;

if (~exist('saveEig','var')),
    saveEig = 0;
end;

if ~exist('removeMean', 'var')
    removeMean = 'yes';
end

if (saveEig == 3),
    disp('Loading Previous EigenDecomposition 1');edit
    load eig_stuff1;
    firstEig = firstEig_t2;lastEig=lastEig_t2;saveEig=saveEig_t2;
    originalDimension = size(data, 1);
else,
    originalDimension = size(data, 1);

    if exist('doTranspose', 'var')
        if strcmp(doTranspose, 'transpose')
            data = data';
        end
        originalDimension = size(data, 2);
    end

    disp(['Calculating Covariance: ' num2str(originalDimension) '^2']);

    % skip detrending if the removing the mean is already done
    if strcmp(lower(removeMean), 'yes')

        disp('Removing mean from the data ...');

        % check if there is any out of memory error
        try
            tempVar = ones(size(data, 1), 1)*mean(data);
            data = data - tempVar;
            clear tempVar;
        catch
            clear tempVar;
            % code to calculate covariance
            disp('Using slow but less memory constrained for calculating covariance');

            % remove mean by doing loop over all the time points
            for i = 1:size(data, 2)
                data(:, i) = data(:, i) - mean(data(:, i));
            end
        end

    end

    covarianceMatrix = data' * data / size(data, 2);  %covariance for time

    clear data;

    % Calculate the eigenvalues and eigenvectors of covariance matrix.
    disp('Calculating eigendecomposition');

    condCovMat = cond(covarianceMatrix);

    disp(['Condition number of covariance matrix is ', num2str(condCovMat)]);

    [V, Lambda] = eig(covarianceMatrix, 'nobalance');

    % projection into the reference  direction  added by Jingyu
    %-----------&---------------
    NewLambda=abs(V'*reference(:).*diag(Lambda));
    %   NewLambda=((V'*reference').^2)./diag(Lambda);

    %-------------------------
    % Sort the eigenvalues - decending.
    disp('Sorting eigenvalues projected into the reference direction');
    [eigenvalues ind] = sort(NewLambda);
    eigenvalues = flipud(eigenvalues);
    ind = flipud(ind);

    if (saveEig == 4),
        disp('Saving Eigen Decomposition 1');
        save eig_stuff1 Lambda V covarianceMatrix eigenvalues;
        clear covarianceMatrix;
    end;
end;

%estimate the rank (Too long for tica!)...
%maxLastEig = rank(covarianceMatrix);
maxLastEig = lastEig;

disp('Selecting Desired Eigenvalues');
% Remove the smaller eigenvalues
%disp('Removing smaller eigenvalues');
if lastEig < originalDimension
    lowerLimitValue = eigenvalues(lastEig);
else
    lowerLimitValue = eigenvalues(originalDimension);
end
lowerColumns = (NewLambda >= lowerLimitValue);

%disp('Removing larger eigenvalues');
% Remove the larger eigenvalues
higherLimitValue = eigenvalues(firstEig);
higherColumns = NewLambda <= higherLimitValue;

% Combine the results from above
selectedColumns = lowerColumns & higherColumns;

% Select the colums which correspond to the desired range
% of eigenvalues.
V_orig = V(ind,ind);
V = ica_fuse_v_selcol(V,selectedColumns);
Lambda_orig = eigenvalues;
Lambda = diag(NewLambda(selectedColumns));

sumAll=sum(eigenvalues);
sumUsed=sum(diag(Lambda));
retained = (sumUsed/sumAll)*100;
fprintf('%g%% of (non-zero) eigenvalues retained.\n', retained);

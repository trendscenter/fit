function [icasig, betaWeights, betaWeightStr, featureCompStr] = ica_fuse_scaleCompData(meanData, icasig, featureDataLength, featureNames, flipFlag, ConvertToZ)
%% Scale components data

%% Load defaults
ica_fuse_defaults;
global DETREND_NUMBER;

if ~exist('flipFlag', 'var')
    flipFlag = 'flip';
end

if ~exist('ConvertToZ', 'var')
    ConvertToZ = 0;
end

numComp = size(icasig, 1);

%% Initialise variables
startPoint = 1; endPoint = 0;

countN = 0;
betaWeights = repmat(struct('value', []), 1, length(featureDataLength));

if (length(featureDataLength) > 1)
    [rSquare_stat, Betas, ModelIndices] = ica_fuse_multipleRegression(detrend(icasig', 0), ...
        detrend(meanData', 0), numComp, 1, length(meanData), DETREND_NUMBER);
    Betas = Betas(ModelIndices);
end


%% Loop over features
for nn = 1:length(featureDataLength)
    % Update the end point
    endPoint = endPoint + featureDataLength(nn);
    
    refData = meanData(startPoint:endPoint);
    % Component data for that feature
    featureData = (icasig(:, startPoint:endPoint))';
    
    % Do a regression fit where model is component data
    [rSquare_stat, b, ModelIndices] = ica_fuse_multipleRegression(detrend(featureData, 0), ...
        (detrend(refData, 0))', numComp, 1, length(refData), DETREND_NUMBER);
    
    % components beta weights
    b = b(ModelIndices);
    
    if (length(featureDataLength) == 1)
        Betas = b;
    end
    
    % make sure there are no incorrect flippings
    b = abs(b);
    
    % Scale components (flip the sign if b is negative)
    for ii = 1:numComp
        countN = countN + 1;
        icasig(ii, startPoint:endPoint) = abs(b(ii))*featureData(:, ii)';
        if (strcmpi(flipFlag, 'flip') && (sign(Betas(ii)) == -1))
            disp(['Flipping sign for component ', num2str(ii), ' of feature ', featureNames{nn}]);
            icasig(ii, startPoint:endPoint) = sign(Betas(ii))*icasig(ii, startPoint:endPoint);
            b(ii) = sign(Betas(ii))*abs(b(ii));
        end
        % Form structure for printing beta weights information
        compIndex = ica_fuse_returnFileIndex(ii);
        if countN == 1
            featureCompStr = [featureNames{nn}, ' ', compIndex];
        else
            featureCompStr = str2mat(featureCompStr, [featureNames{nn}, ' ', compIndex]);
        end
        
    end
    
    clear featureData;
    
    % Beta coefficients (Includes sign also)
    betaWeights(nn).value = b;
    
    if nn == 1
        betaWeightStr = b;
    else
        betaWeightStr = [betaWeightStr; b];
    end
    
    %% Option for converting to z-scores
    if (ConvertToZ)
        fprintf('\n');
        disp(['Converting components of feature ', featureNames{nn}, ' to z-scores']);
        [tempICAsig] = ica_fuse_convertToZScores(icasig(:, startPoint:endPoint));
        icasig(:, startPoint:endPoint) = tempICAsig;
        clear tempICAsig;
        fprintf('\n');
    end
    
    startPoint = endPoint + 1;
    
    clear refData;
    clear b;
end
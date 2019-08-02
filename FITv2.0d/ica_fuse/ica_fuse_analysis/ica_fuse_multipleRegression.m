function [rSquare_stat, b, ModelIndices, otherIndices, linearRegress, removeTrend, ica, modelX] = ...
    ica_fuse_multipleRegression(model, ica, num_Regress, num_DataSets, diffTimePoints, detrendNumber)

% Performs multiple regression
% Outputs: R-Square Statistic, Coefficients, Model Indices, Other than the
% model indices, line fit to model observed data

% inputs are modelTimecourse and icaTimecourse, number of regressors,
% number of data sets

% By default DETRENDNUMBER is 0

% DETRENDNUMBER - case 0 - Removes the mean
% DETRENDNUMBER - case 1 - Removes the mean and linear trend
% DETRENDNUMBER - case 2 - Uses sine and cosine one cycle, removes mean and
% linear trend
% DETRENDNUMBER - case 3 - Uses sine and cosine two cycles plus sine and cosine one cycle, removes mean and
% linear trend


if ~exist('detrendNumber', 'var')
    % use the defaults here
    ica_fuse_defaults;
    global DETREND_NUMBER;
    detrendNumber = DETREND_NUMBER;
end

if nargin < 2
    error('Need atleast two arguments for the multiple linear regression');
end


% Convert to column vector
if size(ica, 1) == 1
    ica = ica';
end

% number of points taken into consideration
nPoints = size(ica, 1);

% check the dimensions
if size(model, 1) ~= nPoints & size(model, 2) ~= nPoints
    error('Model dimensions doesn''t match with that of observed data');
elseif size(model, 1) ~= nPoints
    model = model';
end


% If the number of regressors doesn't exist
if ~exist('num_Regress', 'var')
    num_Regress = size(model, 2);
    %     disp('-- Number of Regressors is not selected. By default checking the size of the design matrix');
end

% By default number of data sets is 1
if ~exist('num_DataSets', 'var')
    num_DataSets = 1;
    %     disp('-- Number of data sets is not selected. By default selecting data sets equal to 1.');
end

if ~exist('detrendNumber', 'var')
    detrendNumber = 0;
    disp('-- detrendNumber doesn''t exist. By default selecting detrendNumber = 0.');
end

if detrendNumber > 3 | round(detrendNumber) - detrendNumber ~=0
    detrendNumber = 0;
    disp('-- detrendNumber selected is not in options or is not a valid integer. By default selecting detrendNumber = 0.');
end

if ~exist('diffTimePoints', 'var')
    % number of Time Points
    numTimePoints = nPoints / num_DataSets;
    numTimePoints = repmat(numTimePoints, 1, num_DataSets);
else
    numTimePoints = diffTimePoints;
end


% % number of rows and columns of the model timecourse
% [nrows, ncols] = size(model);

% Ramp Function
%rampFun = unitRamp((1:numTimePoints)');
%rampFun = detrend(rampFun, 0); % remove the mean before using ramp Function


% A multivariate model of the data is  Multiple regression solves for unknown coefficients ,, and , by performing a
% least squares fit. Construct and solve the set of simultaneous equations by forming the regression matrix, X,
% and solving for the coefficients

%% Concatenate the design matrix using the number of
%% regressors and the data sets and depending on the detrend number

removeTrend = zeros(size(ica));

linearRegress = zeros(size(ica));

% Start grouping the models along with corresponding detrend terms in a design matrix
switch detrendNumber

    % Detrend with 0 or remove mean
    case 0

        nTerms = 1;

        % Total terms including the regressors along second dimension
        totalTerms = nTerms + num_Regress;

        % Initialise modelX that is to be concatenated
        modelX = zeros(nPoints, totalTerms*num_DataSets);

        % Initialise beta coefficients
        b = zeros(totalTerms*num_DataSets, 1);

        rowStart = 1;
        % Loop over number of data sets
        for nDataSets = 1:num_DataSets

            rowEnd = sum(numTimePoints(1:nDataSets));
            % starting Index
            %rowStart = (nDataSets - 1)*numTimePoints + 1;
            % Ending Index
            %rowEnd = nDataSets*numTimePoints;
            %
            X = [model(rowStart:rowEnd, :), ones(numTimePoints(nDataSets), 1)];

            % modelX which contains the concatenated model
            modelX(rowStart:rowEnd, (nDataSets - 1)*size(X, 2) + 1 :  nDataSets*size(X, 2)) = X;
            % calculate other indices
            otherIndex = totalTerms*(nDataSets - 1) + (num_Regress + 1 : nTerms + num_Regress);
            % calculate model indices
            modelIndex = totalTerms*(nDataSets - 1) + (1:num_Regress);
            % solve beta coeff
            b(totalTerms*(nDataSets - 1) + (1:totalTerms), 1) = ica_fuse_regress(ica(rowStart:rowEnd, 1), X);

            % Calculating how much is to be subtracted from time course
            removeTrend(rowStart:rowEnd, 1) = X(:, num_Regress + 1:totalTerms)*b(otherIndex);

            % Trend is removed in timecourse
            ica(rowStart:rowEnd, 1) = ica(rowStart:rowEnd, 1) - removeTrend(rowStart:rowEnd, 1);

            % use the detrended ica time course and solve the coefficients
            % for model
            b(modelIndex, 1) = ica_fuse_regress(ica(rowStart:rowEnd, 1), X(:, 1:num_Regress));

            % calculate linefit
            linearRegress(rowStart:rowEnd, 1) = X(:, 1:num_Regress)*b(modelIndex, 1);

            % update the row starting
            rowStart = 1 + rowEnd;

            clear X;

        end

        % Removes the mean and linear trend
    case 1

        nTerms = 2;
        totalTerms = nTerms + num_Regress;
        modelX = zeros(nPoints, totalTerms*num_DataSets);

        % Initialise beta coefficients
        b = zeros(totalTerms*num_DataSets, 1);

        rowStart = 1;

        % Loop over number of data sets
        for nDataSets = 1:num_DataSets

            rowEnd = sum(numTimePoints(1:nDataSets));
            % starting Index
            %rowStart = (nDataSets - 1)*numTimePoints + 1;
            % Ending Index
            %rowEnd = nDataSets*numTimePoints;
            %
            rampFun = ica_fuse_unitRamp((1:numTimePoints(nDataSets))');
            rampFun = detrend(rampFun, 0); % remove the mean before using ramp Function

            X = [model(rowStart:rowEnd, :), rampFun, ones(numTimePoints(nDataSets), 1)];

            % modelX which contains the concatenated model
            modelX(rowStart:rowEnd, (nDataSets - 1)*size(X, 2) + 1 :  nDataSets*size(X, 2)) = X;

            % calculate other indices
            otherIndex = totalTerms*(nDataSets - 1) + (num_Regress + 1 : nTerms + num_Regress);
            % calculate model indices
            modelIndex = totalTerms*(nDataSets - 1) + (1:num_Regress);
            % solve beta coeff
            b(totalTerms*(nDataSets - 1) + (1:totalTerms), 1) = ica_fuse_regress(ica(rowStart:rowEnd, 1), X);

            % Calculating how much is to be subtracted from time course
            removeTrend(rowStart:rowEnd, 1) = X(:, num_Regress + 1:totalTerms)*b(otherIndex);

            % Trend is removed in timecourse
            ica(rowStart:rowEnd, 1) = ica(rowStart:rowEnd, 1) - removeTrend(rowStart:rowEnd, 1);

            % use the detrended ica time course and solve the coefficients
            % for model
            b(modelIndex, 1) = ica_fuse_regress(ica(rowStart:rowEnd, 1), X(:, 1:num_Regress));

            % calculate linefit
            linearRegress(rowStart:rowEnd, 1) = X(:, 1:num_Regress)*b(modelIndex, 1);

            % update the row starting
            rowStart = 1 + rowEnd;

            clear X;
            clear rampFun;
        end

        % Uses sine and cosine one cycle, removes mean and
        % linear trend
    case 2

        nTerms = 4;
        totalTerms = nTerms + num_Regress;
        modelX = zeros(nPoints, totalTerms*num_DataSets);

        % Initialise beta coefficients
        b = zeros(totalTerms*num_DataSets, 1);

        rowStart = 1;
        % Loop over number of data sets
        for nDataSets = 1:num_DataSets

            rowEnd = sum(numTimePoints(1:nDataSets));

            % steps in sin and cosine
            unitPoint_sin2phi = (2*pi/(numTimePoints(nDataSets) - 1));

            vec2phi = (0:unitPoint_sin2phi:2*pi)';

            rampFun = ica_fuse_unitRamp((1:numTimePoints(nDataSets))');
            rampFun = detrend(rampFun, 0); % remove the mean before using ramp Function


            % starting Index
            %rowStart = (nDataSets - 1)*numTimePoints + 1;
            % Ending Index
            %rowEnd = nDataSets*numTimePoints;
            %

            X = [model(rowStart:rowEnd, :), detrend(sin(vec2phi), 0), detrend(cos(vec2phi), 0), rampFun, ...
                ones(numTimePoints(nDataSets), 1)];

            % modelX which contains the concatenated model
            modelX(rowStart:rowEnd, (nDataSets - 1)*size(X, 2) + 1 :  nDataSets*size(X, 2)) = X;

            % calculate other indices
            otherIndex = totalTerms*(nDataSets - 1) + (num_Regress + 1 : nTerms + num_Regress);
            % calculate model indices
            modelIndex = totalTerms*(nDataSets - 1) + (1:num_Regress);
            % solve beta coeff
            b(totalTerms*(nDataSets - 1) + (1:totalTerms), 1) = ica_fuse_regress(ica(rowStart:rowEnd, 1), X);

            % Calculating how much is to be subtracted from time course
            removeTrend(rowStart:rowEnd, 1) = X(:, num_Regress + 1:totalTerms)*b(otherIndex);

            % Trend is removed in timecourse
            ica(rowStart:rowEnd, 1) = ica(rowStart:rowEnd, 1) - removeTrend(rowStart:rowEnd, 1);

            % use the detrended ica time course and solve the coefficients
            % for model
            b(modelIndex, 1) = ica_fuse_regress(ica(rowStart:rowEnd, 1), X(:, 1:num_Regress));

            % calculate linefit
            linearRegress(rowStart:rowEnd, 1) = X(:, 1:num_Regress)*b(modelIndex, 1);

            % update the row starting
            rowStart = 1 + rowEnd;

            clear X;
            clear rampFun;
        end

        % Uses sine and cosine two cycles plus sine and cosine one cycle, removes mean and
        % linear trend
    case 3

        nTerms = 6;
        totalTerms = nTerms + num_Regress;
        modelX = zeros(nPoints, totalTerms*num_DataSets);

        % Initialise beta coefficients
        b = zeros(totalTerms*num_DataSets, 1);

        rowStart = 1;
        % Loop over number of data sets
        for nDataSets = 1:num_DataSets

            rowEnd = sum(numTimePoints(1:nDataSets));

            % steps in sin and cosine
            unitPoint_sin2phi = (2*pi/(numTimePoints(nDataSets) - 1));

            unitPoint_sin4phi = (4*pi/(numTimePoints(nDataSets) - 1));

            vec2phi = (0:unitPoint_sin2phi:2*pi)';

            vec4phi = (0:unitPoint_sin4phi:4*pi)';

            rampFun = ica_fuse_unitRamp((1:numTimePoints(nDataSets))');
            rampFun = detrend(rampFun, 0); % remove the mean before using ramp Function

            % starting Index
            %rowStart = (nDataSets - 1)*numTimePoints + 1;
            % Ending Index
            %rowEnd = nDataSets*numTimePoints;
            %
            %             X = [model(rowStart:rowEnd, :), sin(vec4phi) + cos(vec4phi), sin(vec2phi) + cos(vec2phi), ...
            %                     rampFun, ones(numTimePoints(nDataSets),
            %                     1)];

            X = [model(rowStart:rowEnd, :), detrend(sin(vec4phi),0), detrend(cos(vec4phi),0), ...
                detrend(sin(vec2phi),0), detrend(cos(vec2phi),0), ...
                rampFun, ones(numTimePoints(nDataSets), 1)];
            % modelX which contains the concatenated model
            modelX(rowStart:rowEnd, (nDataSets - 1)*size(X, 2) + 1 :  nDataSets*size(X, 2)) = X;

            % calculate other indices
            otherIndex = totalTerms*(nDataSets - 1) + (num_Regress + 1 : nTerms + num_Regress);
            % calculate model indices
            modelIndex = totalTerms*(nDataSets - 1) + (1:num_Regress);
            % solve beta coeff
            b(totalTerms*(nDataSets - 1) + (1:totalTerms), 1) = ica_fuse_regress(ica(rowStart:rowEnd, 1), X);

            % Calculating how much is to be subtracted from time course
            removeTrend(rowStart:rowEnd, 1) = X(:, num_Regress + 1:totalTerms)*b(otherIndex);

            % Trend is removed in timecourse
            ica(rowStart:rowEnd, 1) = ica(rowStart:rowEnd, 1) - removeTrend(rowStart:rowEnd, 1);

            % use the detrended ica time course and solve the coefficients
            % for model
            b(modelIndex, 1) = ica_fuse_regress(ica(rowStart:rowEnd, 1), X(:, 1:num_Regress));

            % calculate linefit
            linearRegress(rowStart:rowEnd, 1) = X(:, 1:num_Regress)*b(modelIndex, 1);

            % update the row starting
            rowStart = 1 + rowEnd;

            clear X;
            clear rampFun;
        end

end
% end switch for detrendNumber

% collect the indices other than the model indices
otherIndices = zeros(1, nTerms*num_DataSets); ModelIndices = zeros(1, num_Regress*num_DataSets);

% Indices other than the model
for ii = 1:num_DataSets
    otherIndices(1, (ii-1)*nTerms + 1 : nTerms*ii) = totalTerms*(ii - 1) + (num_Regress + 1 : nTerms + num_Regress);
    ModelIndices(1, (ii-1)*num_Regress + 1 : num_Regress*ii) = totalTerms*(ii - 1) + (1:num_Regress);
end


% calculate r-square statistic
residual = ica - modelX(:, ModelIndices)*b(ModelIndices, 1);
rSquare_stat = 1 - (norm(residual, 2) / norm(ica - mean(ica), 2))^2;

% % return coefficients and R^2 Statistic
% [b, rSquare_stat] = ica_fuse_regress(ica, modelX);
%
% % Initialise removeTrend
% removeTrend = zeros(size(ica));
%
% % Calculating how much is to be subtracted from time course
% removeTrend = modelX(:, otherIndices)*b(otherIndices);
%
% % Trend is removed in timecourse
% ica = ica - removeTrend;

% Calculating again to get the coefficients after detrending
% [modelB, rSquare_stat] = ica_fuse_regress(ica, modelX(:, ModelIndices));

% calculate linefit
%linearRegress = modelX(:, ModelIndices)*modelB;

% Update the new regression coefficients
%b(ModelIndices) = modelB;
%clear modelB;


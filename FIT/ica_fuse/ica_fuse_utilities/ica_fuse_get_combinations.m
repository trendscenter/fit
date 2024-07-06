function all_comb = ica_fuse_get_combinations(numFeatures, stackAllFeatures, optimalFeatures)
%% Function to get all the combinations of features
%
% Inputs:
% 1. numFeatures - Number of features
% 2. stackAllFeatures - Options are 0 and 1.
% 3. optimalFeatures - Options are 0 and 1.
%
% Outputs:
% all_comb - All combinations cell array
%

if ~exist('stackAllFeatures', 'var')
    stackAllFeatures = 1;
end

if ~exist('optimalFeatures', 'var')
    optimalFeatures = 0;
end


% Initialise all combinations
all_comb{1} = [];

%% Get all combinations possible. First combination is to stack all the
% features
if (numFeatures == 1)
    all_comb{1} = 1;
elseif (numFeatures == 2)
    if (~optimalFeatures)
        all_comb{1} = [1, 2];
    else
        all_comb = {[1, 2]; 1; 2};
    end
else
    if (stackAllFeatures)
        all_comb{1} = (1:numFeatures);
    end
    if optimalFeatures
        %% Pairwise combinations and individual combinations
        pairwise_comb = num2cell(nchoosek(1:numFeatures, 2), 2);
        individual_comb = num2cell(nchoosek(1:numFeatures, 1), 2);
        all_comb = [all_comb; pairwise_comb; individual_comb];
    end
end


if (length(all_comb) == 1) && isempty(all_comb{1})
    error('Error:StackAllFeatures', ['There is no combination of features to run. Please check if you have turned off both\nthe global variables', ...
        ' STACK_ALL_FEATURES and OPTIMIZE_FEATURES in ica_fuse_defaults.m']);
end

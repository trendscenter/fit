function combNumber = ica_fuse_getCombNum(featureNames, currentCombName)
%% Get current combination number using the current combination name and
% feature names
%
% Inputs:
% 1. featureNames - Feature names
% 2. currentCombName - Current combination name
%
% Outputs:
% combNumber - Combination number
% 

matchedStrs = strread(currentCombName, '%s', 'delimiter', '&');

try
    matchedStrs = strtrim(matchedStrs);
catch
    matchedStrs = deblank(matchedStrs);
end

combNumber = zeros(1, length(matchedStrs));

% Loop over matched strings
for nM = 1:length(matchedStrs)
    tempMatch = strmatch(matchedStrs{nM}, featureNames, 'exact');
    combNumber(nM) = tempMatch(1);
end
% End loop over matched strings
function [dat] = ica_fuse_combineFeatures(dataN)
% Combine data over features

ica_fuse_defaults;
global OPTIMIZE_FEATURES;


numFeatures = length(dataN);

numAll = 1;
if strcmpi(OPTIMIZE_FEATURES, 'yes')
    
    if numFeatures == 1   
        numCombinations = 0;
        numIndividual = 0;
    elseif numFeatures == 2
        numCombinations = 0;
        numIndividual = 2;
    else
        numCombinations = (numFeatures)*(numFeatures - 1) / 2;
        numIndividual = numFeatures;
    end
    
else
    
    numCombinations = 0;
    numIndividual = 0;
    
end

% loop over features
for ii = 1:numFeatures
    if ii == 1    
        data = dataN(ii).data;
        combName = dataN(ii).feature_name;
        featureDataLength = size(data, 2);
    else        
        data = [data, dataN(ii).data];
        featureDataLength = [featureDataLength, size(dataN(ii).data, 2)];
        combName = [combName, ' & ', dataN(ii).feature_name];        
    end
    files(ii).name = dataN(ii).files;
    files(ii).feature_name = dataN(ii).feature_name;    
end

% Attach to data structure combination name and data
countN = 1;
dat(countN).data = data;
dat(countN).combName = combName;
dat(countN).files = files;
dat(countN).featureDataLength = featureDataLength;
clear featureDataLength;
clear combName; clear data; clear files;

for nn = 1:numFeatures   
    dat(countN).timeAxis(nn).data = dataN(nn).xAxis;
end

if numCombinations > 0

    for ii = 1 : numFeatures - 1        
        for jj = ii + 1 : numFeatures 
            countN = countN + 1;
            % stack these two together
            dat(countN).data = [dataN(ii).data, dataN(jj).data];          
            dat(countN).files(1).name = dataN(ii).files; dat(countN).files(2).name = dataN(jj).files;            
            % Store features also
            dat(countN).files(1).feature_name = dataN(ii).feature_name; 
            dat(countN).files(2).feature_name = dataN(jj).feature_name;            
            dat(countN).combName = [dataN(ii).feature_name, ' & ', dataN(jj).feature_name];            
            dat(countN).featureDataLength = [size(dataN(ii).data, 2), size(dataN(jj).data, 2)];
        end
    end
    
end

if numIndividual  > 0
    
    for ii = 1:numFeatures
        countN = countN + 1;
        dat(countN).data = dataN(ii).data;          
        dat(countN).files(1).name = dataN(ii).files;
        % store feature name
        dat(countN).files(1).feature_name = dataN(ii).feature_name;
        dat(countN).combName = dataN(ii).feature_name;        
        dat(countN).featureDataLength = size(dataN(ii).data, 2);
    end
    
end




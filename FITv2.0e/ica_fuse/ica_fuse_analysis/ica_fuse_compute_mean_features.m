function [meanData, meanDataGroups] = ica_fuse_compute_mean_features(dataN, dataLength, numSubjects)
%% Compute mean for each feature and mean for each group and feature
%
% Inputs:
% 1. data - Data for features
% 2. dataLength - Data Length
% 3. numSubjects - Number of subjects for the selected groups
%
% Outputs:
% 1. meanData - Mean data for each feature.
% 2. meanDataGroups - Mean data for each feature and group in a cell array
%

meanDataGroups = [];

%% Compute mean for each feature
meanData = zeros(1, sum(dataLength));
startTp = 1;
endTp = 0;
% Loop over features
for nF = 1:length(dataN)
    endTp = endTp + dataLength(nF);
    tmpM = ica_fuse_remove_mean((dataN(nF).data)');
    tmpM = tmpM';
    meanData(startTp:endTp) = mean(tmpM);
    startTp = endTp + 1;
end
% Loop over features

clear startTp endTp tmpM;

%% Compute mean for each feature and group
if (length(numSubjects) > 1)
    meanDataGroups = cell(length(numSubjects), 1);
    startGroup = 1;
    endGroup = 0;
    % Loop over groups
    for nG = 1:length(numSubjects)
        endGroup = endGroup + numSubjects(nG);
        groupInd = (startGroup:endGroup);
        startTp = 1;
        endTp = 0;
        tempData = zeros(1, sum(dataLength));
        % Loop over features
        for nF = 1:length(dataN)
            endTp = endTp + dataLength(nF);
            tmpM = ica_fuse_remove_mean((dataN(nF).data(groupInd, :))');
            tmpM = tmpM';
            %tempData(startTp:endTp) = mean(dataN(nF).data(groupInd, :));
            tempData(startTp:endTp) = mean(tmpM);
            startTp = endTp + 1;
        end
        % End loop over features
        meanDataGroups{nG} = tempData;
        clear tempData;
        startGroup = endGroup + 1;
    end
    % End loop over groups
end





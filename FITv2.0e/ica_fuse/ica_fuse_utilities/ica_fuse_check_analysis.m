function [combnsToRun, stepsToRun] = ica_fuse_check_analysis(fusionInfo)
%% Figure out where the analysis stopped. Return combinations and steps
% that need to be run.
%
% Inputs:
% 1. fusionInfo - Fusion information
%
% Outputs:
% 1. combnsToRun - Combination/Combinations to be run
% 2. stepsToRun - Starting analysis step that needs to be run.
%

%% Initialise output variables
combnsToRun = [];
stepsToRun = [];

%% All combinations
all_comb = fusionInfo.run_analysis.all_comb;

%% Get required variables
pcaFiles = getVar(fusionInfo, 'pcaFiles', all_comb);
icaFiles = getVar(fusionInfo, 'icaFiles', all_comb);
backReconstructFiles = getVar(fusionInfo, 'backReconstructFiles', all_comb);

%% Scaling step error check
if ((~isfield(fusionInfo.run_analysis, 'scaleCompFiles')) || (~isfield(fusionInfo.run_analysis, 'outputFiles')))
    scaleCompFiles = [];
else
    scaleCompFiles = fusionInfo.run_analysis.scaleCompFiles;
end

%% Maximum no. of loops
maxLoops = length(all_comb);
if (maxLoops > 1)
    maxLoops = 2;
end

maxSteps = 5;

files = {pcaFiles, icaFiles, backReconstructFiles, scaleCompFiles};

%% Check where the analysis stopped
for nC = 1:maxLoops

    combNumber = all_comb{nC};

    if isempty(combNumber)
        continue;
    end

    if (nC == 1)

        for nStep = 1:length(files)
            [combnsToRun, stepsToRun] = getCombNo(files{nStep}, nC, nStep + 1, all_comb, maxSteps);
            if ~isempty(combnsToRun)
                return;
            end
        end

    else
        analysisSteps = [length(pcaFiles), length(icaFiles), length(backReconstructFiles)];
        %% Min value and index
        [minVal, index] = min(analysisSteps);

        %% Return the combinations and steps to  be run
        if (minVal < length(all_comb))
            combnsToRun = (minVal + 1:length(all_comb));
            stepsToRun = (index + 1:maxSteps);
            return;
        end
    end

end
%% End for checking where the analysis stopped

function out1 = getVar(fusionInfo, fieldName, all_comb)
%% Get output variable

out1 = [];

if (isfield(fusionInfo.run_analysis, fieldName))
    out1 = getfield(fusionInfo.run_analysis, fieldName);
end

errCheck(all_comb, out1);

function errCheck(all_comb, filesStruct)
%% Error check

if (length(all_comb) < length(filesStruct))
    error('Error:Defaults', ['The default settings for this analysis was changed.\nPlease reset ica_fuse_defaults.m to the original settings ', ...
        ' if you want to resume analysis']);
end



function [combnsToRun, stepsToRun] = getCombNo(fileStruct, combNumber, analysisStep, all_comb, maxSteps)
%% Get combination number where the analysis stopped

combnsToRun = [];
stepsToRun = [];

isInterrupted = isempty(fileStruct);

if (isInterrupted)
    combnsToRun = (combNumber:length(all_comb));
    stepsToRun = (analysisStep:maxSteps);
end
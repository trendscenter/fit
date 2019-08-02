function [optimization_data, legendString, sorted_ind] = ica_fuse_rank_features(sortParameters)
%% Function to rank features.
%
% Inputs:
% 1. sortParameters - Sort parameters data structure
%
% Outputs:
% optimization_data - Optimization data
%

if ~isfield(sortParameters, 'outputDir')
    outputDir = pwd;
else
    outputDir = sortParameters.outputDir;
end

histCriteria = sortParameters.histogramCriteria;
z_threshold = sortParameters.z_threshold;
selGroups = sortParameters.selGroupsVal;
featureDataLength = sortParameters.featureDataLength;
back_reconstruct_files = sortParameters.backReconstructFiles;
div_name = sortParameters.div_name;
div_value = sortParameters.div_value;
groupNames = sortParameters.groupNames;
sortingCriteria = sortParameters.sortingCriteria;
all_comb = sortParameters.all_comb;

group1Name = deblank(groupNames(selGroups(1), :));
group2Name = deblank(groupNames(selGroups(2), :));
sortParameters.selGroupNames = str2mat(group1Name, group2Name);

legendString = cell(length(back_reconstruct_files), 1);

% Initialise optimization_data structure
optimization_data = repmat(struct('divergence', [], 'combinationName', [], 'bestComp', []), 1, length(all_comb));

% good cells
good_inds = ica_fuse_good_cells(all_comb);

countData = 0;
% Loop over number of back reconstruction files
for nn = 1:length(all_comb)

    countData = countData + 1;

    if isempty(all_comb{nn})
        continue;
    end

    % Use Two sample t-test on mixing coefficients
    if strcmpi(sortingCriteria, 'ttest2')

        if strcmpi(histCriteria, 'feature')
            [tempHistData, meanHist] = ica_fuse_calculate_histogram(sortParameters, histCriteria, z_threshold, nn);
            out1 = meanHist(1).data{1};
            out2 = meanHist(1).data{2};
        else
            tempHistData = ica_fuse_calculate_histogram(sortParameters, histCriteria, z_threshold, nn);
            out1 = tempHistData(1).group(1).data;
            out2 = tempHistData(1).group(2).data;
        end

        bestComp = tempHistData(1).bestComp;
        combinationName = tempHistData(1).combinationName;

        % Calculate divergence
        [dName, spatial_divergence] = ica_fuse_divergence(out1, out2, div_name, div_value);

        clear out1 out2;

    else
        % Use spatial divergence

        load(fullfile(outputDir, back_reconstruct_files(nn).name), 'combinationName');
        sortParameters.selectedFeatureVal = 1:length(featureDataLength(nn).Length);
        sortParameters.selectedFeature = str2mat(strread(combinationName, '%s', 'delimiter', '&'));

        % Sort components based on spatial divergence
        sortResults = ica_fuse_sort_components(sortParameters, nn, 0);
        bestComp = sortResults.sorted_comp(1);
        spatial_divergence = sortResults.values(1);
        clear sortResults ;

    end

    optimization_data(countData).divergence = spatial_divergence;
    optimization_data(countData).combinationName = combinationName;
    optimization_data(countData).bestComp = bestComp;
    clear tempHistData;
    if exist('meanHist', 'var')
        optimization_data(countData).meanHist = meanHist;
    end

    legendString{countData} = [num2str(countData), ': ', combinationName];

end
% end loop over number of back reconstruction files

% Use good indices
optimization_data = optimization_data(good_inds);
legendString = legendString(good_inds);

% Sort the KL divergence
[divrgs, sorted_ind] = sort([optimization_data.divergence]);

sorted_ind = sorted_ind(end:-1:1);


% Sort the optimization data
optimization_data = optimization_data(sorted_ind);
legendString = legendString(sorted_ind);

% Fixed in the updates of FITv2.0b (May 14, 2009) when all the features are
% not stacked
good_inds = find(good_inds ~= 0);
sorted_ind = good_inds(sorted_ind);

fprintf('\n');
if strcmpi(sortParameters.sortingCriteria, 'divergence')
    disp('Best component based on spatial divergence are as follows: ');
else
    disp('Best component based on two sample t-test of mixing coefficients are as follows: ');
end
% display best component numbers
for nn = 1:length(sorted_ind)
    disp([optimization_data(nn).combinationName, ': ', ...
        num2str(optimization_data(nn).bestComp)]);
end
fprintf('\n');

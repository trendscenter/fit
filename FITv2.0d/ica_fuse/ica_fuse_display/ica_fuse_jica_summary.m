function ica_fuse_jica_summary(param_file, displayParameters)
% Joint ICA report generator

ica_fuse_defaults;
global ANATOMICAL_FILE;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGE_VALUES;
global ANATOMICAL_PLANE;
global FUSE_PARAM_FILE;


if (~exist('param_file', 'var'))
    param_file = FUSE_PARAM_FILE;
end


if (~isempty(param_file))
    load(param_file);
end


try
    imageValues = displayParameters.image_values;
catch
    imageValues = IMAGE_VALUES;
end

try
    threshold = abs(displayParameters.z_threshold);
catch
    threshold = Z_THRESHOLD;
end

try
    anatomical_plane = lower(displayParameters.anatomical_plane);
catch
    anatomical_plane = ANATOMICAL_PLANE;
end

try
    slices_in_mm = displayParameters.slices_in_mm;
catch
    slices_in_mm = -40:4:72;
end

try
    convertToZ = displayParameters.convert_to_z;
catch
    convertToZ = strcmpi('yes', CONVERT_TO_Z);
end

convertToZStr = 'no';
if (convertToZ)
    convertToZStr = 'yes';
end

disp_opts = {'positive and negative', 'positive', 'absolute value', 'negative'};
if (isnumeric(imageValues))
    imageValueStr = disp_opts{imageValues};
else
    imageValueStr = imageValues;
end
fusionInfo.run_analysis.image_values = imageValueStr;
imageValues = strmatch(lower(imageValueStr), disp_opts, 'exact');


modalities = cellstr(char(fusionInfo.run_analysis.dataInfo(1).feature.modality));
groupNames = cellstr(char(fusionInfo.run_analysis.dataInfo.name));
featureNames = cellstr(char(fusionInfo.run_analysis.dataInfo(1).feature.name));


fusionInfo.run_analysis.modalities = modalities;
fusionInfo.run_analysis.groupNames = groupNames;
fusionInfo.run_analysis.featureNames = featureNames;
fusionInfo.run_analysis.threshold = threshold;
fusionInfo.run_analysis.convert_to_z = convertToZStr;
fusionInfo.run_analysis.anatomical_plane = anatomical_plane;
fusionInfo.run_analysis.anatomical_file = ANATOMICAL_FILE;
fusionInfo.run_analysis.slices_in_mm = slices_in_mm;

cm = [];
text_left_right = [];
outputDir = pwd;

outputDir = fusionInfo.run_analysis.outputDir;
numFeatures = fusionInfo.run_analysis.numFeatures;
selectedFeatureVal = (1:numFeatures);
dataInfo = fusionInfo.run_analysis.dataInfo;
numSubjects = fusionInfo.run_analysis.numSubjects;
outputFiles = fusionInfo.run_analysis.outputFiles;
selectedComp = (1:fusionInfo.run_analysis.numComp);
mask_ind = fusionInfo.run_analysis.mask_ind;
[dims, voxels] = ica_fuse_getFeatureDIM(mask_ind);


%% Fusion Parameters
%
dispParaStr = dispPara(fusionInfo);
fprintf('\n\n');
disp(char(dispParaStr));
fprintf('\n\n');

%%

for nOutputFiles = 1:length(selectedFeatureVal)
    % component files
    componentFiles = ica_fuse_fullFile('directory', outputDir, 'files', ...
        outputFiles(selectedFeatureVal(nOutputFiles)).name);
    if isfield(outputFiles(selectedFeatureVal(nOutputFiles)), 'groupFiles')
        groupCompFiles = outputFiles(selectedFeatureVal(nOutputFiles)).groupFiles;
    else
        groupCompFiles = [];
    end
    feature_name = outputFiles(selectedFeatureVal(nOutputFiles)).feature_name;
    
    featureInfo.feature_name = feature_name;
    featureInfo.groupNames = char(dataInfo.name);
    featureInfo.numSubjects = numSubjects;
    %featureInfo.numSubjects = size(str2mat(dataInfo(1).feature(1).files.name), 1);
    
    % Get the input files
    for nInputFiles = 1:length(dataInfo)
        
        if nInputFiles == 1
            feature_input_files = char(dataInfo(nInputFiles).feature(selectedFeatureVal(nOutputFiles)).files.name);
        else
            feature_input_files = char(feature_input_files, ...
                char(dataInfo(nInputFiles).feature(selectedFeatureVal(nOutputFiles)).files.name));
        end
        
    end
    % end for getting feature files
    
    [compData, anatData, meanData, groupCompData, meanDataLegend, groupCompLegend, text_left_right, ...
        plotType, HInfo] = ica_fuse_loadCompData('component_files', componentFiles, ...
        'group_component_files', groupCompFiles, 'slices_in_mm', slices_in_mm, 'anatomical_file', ...
        ANATOMICAL_FILE, 'anatomical_view', anatomical_plane, 'component_numbers', ...
        selectedComp, 'interp_message', ['Interpolating components ', ...
        'of feature ', feature_name], 'interp_title', ['Interp comp of ', feature_name], ...
        'input_files', feature_input_files, 'voxels', voxels, 'feature_info', featureInfo, ...
        'mask_ind', mask_ind(selectedFeatureVal(nOutputFiles)).ind, 'outputDir', outputDir, ...
        'flip_analyze_images', 1);
    
    
    %     if ~isfield(displayParameters, 'text_left_right')
    %         displayParameters.text_left_right = text_left_right;
    %     end
    
    clear feature_input_files;
    clear maxICAIm; clear minICAIm;
    maxICAIm = zeros(1, length(selectedComp));
    minICAIm = zeros(1, length(selectedComp));
    minInterval = 0;
    maxInterval = 0;
    
    % For images apply display parameters
    if strcmpi(plotType, 'image')
        structDIM = HInfo.DIM;
        % Apply display parameters
        compData = ica_fuse_applyDispParameters(plotType, compData, convertToZ, imageValues, threshold);
        
        % create composite image
        [compData, minICAIm, maxICAIm, minInterval, maxInterval] = ica_fuse_create_composite(anatData, compData, ...
            imageValues, anatomical_plane, structDIM, threshold);
        colorbarLim = [minInterval, maxInterval];
        % Get the associated colormap
        cm = ica_fuse_getColormap(imageValues, 1);
        clear anatData;
        
    end
    
    % Store the component data in outputData structure
    for nComp = 1:length(selectedComp)
        outputData(nOutputFiles).CompData(nComp).data = squeeze(compData(nComp, :, :));
        compIndex = ica_fuse_returnFileIndex(selectedComp(nComp));
        if ~isempty(groupCompData)
            outputData(nOutputFiles).CompData(nComp).groupCompData = squeeze(groupCompData(:, :, nComp, :));
        else
            outputData(nOutputFiles).CompData(nComp).groupCompData = [];
        end
        if ~isempty(groupCompData)
            outputData(nOutputFiles).CompData(nComp).groupCompLegend = groupCompLegend;
        else
            outputData(nOutputFiles).CompData(nComp).groupCompLegend = {};
        end
        outputData(nOutputFiles).CompData(nComp).axesTitle = ['Comp ', compIndex, ' Feature ', feature_name];
        outputData(nOutputFiles).CompData(nComp).plotType = plotType;
        outputData(nOutputFiles).CompData(nComp).maxICAIm = maxICAIm(nComp);
        outputData(nOutputFiles).CompData(nComp).minICAIm = minICAIm(nComp);
        outputData(nOutputFiles).CompData(nComp).minInterval = minInterval;
        outputData(nOutputFiles).CompData(nComp).maxInterval = maxInterval;
        %         if isfield(displayParameters, 'text_left_right')
        %             outputData(nOutputFiles).CompData(nComp).textLeftRight = displayParameters.text_left_right;
        %         else
        outputData(nOutputFiles).CompData(nComp).textLeftRight = [];
        %     end
    end
    outputData(nOutputFiles).meanData = meanData;
    outputData(nOutputFiles).meanDataLegend = meanDataLegend;
    clear meanData;
    clear meanDataLegend;
    clear featureInfo;
    clear compData;
    
end

% Change the data structure such that the components are plotted by feature
countData = 0;

for nComp = 1:length(selectedComp)
    
    if exist('sortCompData', 'var')
        % Add Mixing Coefficients at the top
        countData = countData + 1;
        plotData(countData).data = sortCompData(nComp).data;
        compIndex = ica_fuse_returnFileIndex(selectedComp(nComp));
        plotData(countData).axesTitle = sortCompData(nComp).title;
        plotData(countData).plotType = sortResults.plotType;
        plotData(countData).colorbarMinMaxText = sortCompData(nComp).colorbarText;
        plotData(countData).colorbarLim = sortCompData(nComp).colorbarLim;
        plotData(countData).groupNames = sortResults.selGroupNames;
        plotData(countData).selFeatureNames = sortResults.selFeatureNames;
        % Selected groups, features, component
        plotData(countData).selGroupsVal = sortResults.selGroupsVal;
        plotData(countData).selFeaturesVal = sortResults.selFeaturesVal;
        plotData(countData).compNum = selectedComp(nComp);
        plotData(countData).fusionFile = fusionFile;
        plotData(countData).threshold = threshold;
        %             if isfield(sortCompData, 'histData')
        %                 plotData(countData).histData = sortCompData(nComp).histData;
        %             end
    end
    
    % Plot For a component the features
    for nOutputFiles = 1:length(outputData)
        countData = countData + 1;
        plotData(countData).data = outputData(nOutputFiles).CompData(nComp).data;
        plotData(countData).groupCompData = outputData(nOutputFiles).CompData(nComp).groupCompData;
        plotData(countData).groupCompLegend = outputData(nOutputFiles).CompData(nComp).groupCompLegend;
        plotData(countData).meanData = outputData(nOutputFiles).meanData;
        plotData(countData).meanDataLegend = outputData(nOutputFiles).meanDataLegend;
        plotData(countData).axesTitle = outputData(nOutputFiles).CompData(nComp).axesTitle;
        plotData(countData).plotType = outputData(nOutputFiles).CompData(nComp).plotType;
        minICAIm = outputData(nOutputFiles).CompData(nComp).minICAIm;
        maxICAIm = outputData(nOutputFiles).CompData(nComp).maxICAIm;
        plotData(countData).colorbarMinMaxText = str2mat(num2str(minICAIm), num2str(maxICAIm));
        plotData(countData).colorbarLim = [outputData(nOutputFiles).CompData(nComp).minInterval, ...
            outputData(nOutputFiles).CompData(nComp).maxInterval];
        plotData(countData).groupNames = '';
        plotData(countData).textLeftRight = outputData(nOutputFiles).CompData(nComp).textLeftRight;
    end
    
end
clear outputData;


%% Components
ica_fuse_display_features('plot_data', plotData, 'color_map', cm, 'title_color', 'c', 'time_course_color', '.-c', ...
    'ica_loading_color', {'g'; 'r'; 'c'; 'm'; 'b'}, 'number_per_figure', 1, 'slice_plane', anatomical_plane, 'plot_prev_next', 0);


pause(1);

global figureData;

figureData = fusionInfo.run_analysis;
figureData.output_prefix = fusionInfo.run_analysis.prefix;
figureData.featureNames = char(fusionInfo.run_analysis.dataInfo(1).feature.name);
figureData.groupNames = char(fusionInfo.run_analysis.dataInfo.name);
figureData.all_comb = fusionInfo.run_analysis.all_comb;
figureData.numSubjects = fusionInfo.run_analysis.numSubjects;
[dims, figureData.voxels] = ica_fuse_getFeatureDIM(fusionInfo.run_analysis.mask_ind);


% if (fusionInfo.run_analysis.numGroups > 1)
% T = evalc('get_sort_info');
[T, opt_results, sortR, divR] = evalc('get_sort_info');

%% Optimal features
% ICA is run on different combination of features. Best component for each feature combination is determined using two sample t-test on the mixing coefficients.
% Best component voxels are thresholded using the default Z-threshold and applied as a mask to the original data. ND task histogram is generated for each group and divergence between the
% groups is computed.
if (size(figureData.all_comb, 1) > 1)
    for n = 1:length(opt_results)
        plotBars(opt_results{n}.values, opt_results{n}.xlabelstr, opt_results{n}.ylabelstr, opt_results{n}.titlestr, opt_results{n}.legendstr);
    end
end

%% Two sample t-test results
% Two sample t-test is computed on the mixing coefficents between the selected groups.

% * *a) pval* - p-values
% * *b) tval* - T-values
%

components = (1:figureData.numComp)';
for n = 1:length(sortR)
    plotBars(sortR{n}.tval, sortR{n}.xlabelstr, sortR{n}.ylabelstr, sortR{n}.titlestr, sortR{n}.legendstr);
    disp(sortR{n}.titlestr);
    pval = sortR{n}.pval;
    tval = sortR{n}.tval;
    varNames = {'ComponentNumber', 'Pvalues', 'Tvalues'};
    try
        T = table(components, pval(:), tval(:), 'VariableNames', varNames);
        disp(T);
    catch
        fprintf('%20s\t%20s\t%20s\n', varNames{:});
        fprintf('\n');
        for nComp = 1:length(components)
            fprintf('%20s\t%20s\t%20s\n', num2str(components(nComp), '%d'), num2str(pval(nComp), '%0.3f'), num2str(tval(nComp), '%0.3f'));
        end
        fprintf('\n');
    end
end


%% Divergence results
% For each component, Z-threshold is applied and ND task histogram is
% computed for each group. Divergence is computed between the
% selected groups for each component.
%
for n = 1:length(divR)
    plotBars(divR{n}.values, divR{n}.xlabelstr, divR{n}.ylabelstr, divR{n}.titlestr, divR{n}.legendstr);
    div_values = divR{n}.values;
    varNames = {'ComponentNumber', 'Divergence'};
    try
        T = table(components, div_values(:), 'VariableNames', varNames);
        disp(T);
    catch
        fprintf('%20s\t%20s\t%20s\n', varNames{:});
        fprintf('\n');
        for nComp = 1:length(components)
            fprintf('%20s\t%20s\n', num2str(components(nComp), '%d'), num2str(div_values(nComp), '%0.3f'));
        end
        fprintf('\n');
    end
end


function [opt_results, ttestR, divR] = get_sort_info


opt_results = [];
ttestR = [];
divR = [];

global figureData;

ica_fuse_defaults;

global Z_THRESHOLD_HISTOGRAM;

try
    z_threshold = Z_THRESHOLD_HISTOGRAM;
catch
    z_threshold = 1.5;
end


if (figureData.numGroups == 1)
    return;
end

%figureData = fusionInfo.run_analysis;

[div_name, div_value] = ica_fuse_getDivergencePara;
figureData.div_name = div_name;
figureData.div_value = div_value;
figureData.histogramCriteria = 'feature';
figureData.sortingCriteria = 'ttest2';
%figureData.selGroupsVal = [1, 2];
figureData.z_threshold = z_threshold;
%figureData.selGroupNames = figureData.groupNames;
%figureData.selGroupsVal = (1:size(figureData.selGroupNames, 1));
allGroups = nchoosek(1:figureData.numGroups, 2);
% Selected feature names and indices
figureData.selectedFeature = figureData.featureNames;
figureData.selectedFeatureVal = (1:size(figureData.selectedFeature, 1));

opt_results = cell(1, size(allGroups, 1));
%ttestR = cell(1, size(allGroups, 1));
ttestR = {};
divR = cell(1, size(allGroups, 1));
for nC = 1:size(allGroups, 1)
    figureData.selGroupsVal = allGroups(nC, :);
    figureData.selGroupNames = char(deblank(figureData.groupNames(figureData.selGroupsVal(1), :)), deblank(figureData.groupNames(figureData.selGroupsVal(2), :)));
    [optimization_data, legendString, sorted_ind] = ica_fuse_rank_features(figureData);
    %fH = figure('color', 'w');
    tmp_dvs = [optimization_data.divergence]';
    tmp_combnames = cellstr(char(optimization_data.combinationName));
    tmp.values = tmp_dvs(:)';
    tmp.legendstr = tmp_combnames;
    tmp.xlabelstr = '';
    tmp.ylabelstr = 'Divergence values';
    tmp.titlestr = [deblank(figureData.selGroupNames(1, :)), ' vs ', deblank(figureData.selGroupNames(2, :))];
    opt_results{nC} = tmp;
    
    % Two sample t-test
    figureData.sortingCriteria = 'ttest2';
    sortResults = ica_fuse_sort_components(figureData);
    [dd, sorted_inds] = sort([sortResults.sorted_comp]);
    sortData = sortResults.sortCompData(sorted_inds);
    %     pval = zeros(1, length(sortData));
    %     tval = pval;
    for n = 1:length(sortData)
        ga = sortData(n).data(1).dat;
        gb = sortData(n).data(2).dat;
        if (n == 1)
            pval = zeros(size(ga, 2), length(sortData));
            tval = pval;
        end
        for nP = 1:size(ga, 2)
            [pval(nP, n), tval(nP, n)] = ica_fuse_ttest2(ga(:, nP), gb(:, nP));
        end
    end
    
    for nP = 1:size(pval, 1)
        tmp2.pval = pval(nP, :);
        tmp2.tval = tval(nP, :);
        tmp2.titlestr = [deblank(figureData.selGroupNames(1, :)), ' vs ', deblank(figureData.selGroupNames(2, :))];
        if (size(pval, 1) > 1)
            tmp2.titlestr = [tmp2.titlestr, ' Of Feature ', deblank(figureData.selectedFeature(figureData.selectedFeatureVal(nP), :))];
        end
        tmp2.xlabelstr = 'Components';
        tmp2.ylabelstr = 'T-values';
        tmp2.legendstr = cellstr(strcat('Comp ', num2str((1:length(sortData))')));
        ttestR{end + 1} = tmp2;
    end
    
    figureData.sortingCriteria = 'divergence';
    sortResults = ica_fuse_sort_components(figureData);
    [dd, sorted_inds] = sort([sortResults.sorted_comp]);
    clear tmp2;
    div_vals = sortResults.values(sorted_inds);
    tmp2.values = div_vals;
    tmp2.titlestr = [deblank(figureData.selGroupNames(1, :)), ' vs ', deblank(figureData.selGroupNames(2, :))];
    tmp2.xlabelstr = 'Components';
    tmp2.ylabelstr = 'Divergence';
    tmp2.legendstr = cellstr(strcat('Comp ', num2str((1:length(sortData))')));
    divR{nC} = tmp2;
    % plotBars(tmp_dvs, '', 'Divergence values', [deblank(figureData.selGroupNames(1, :)), 'vs ', deblank(figureData.selGroupNames(2, :))], tmp_combnames);
end




function plotBars(vals, xlabelstr, ylabelstr, titlestr, legendstr)


sz = get(0, 'ScreenSize');
defaultFigPos = get(0, 'DefaultFigurePosition');
defaultFigPos(3) = 0.9*sz(3);
defaultFigPos(2) = 0.5*sz(4) - 0.5*defaultFigPos(4) + 10;
defaultFigPos(1) = 0.5*sz(3) - 0.5*defaultFigPos(3) + 10;

fH = figure('color', 'w', 'position', defaultFigPos);
gH = axes('parent', fH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
cmap = winter(64);

cmap_inds = round(linspace(1, size(cmap, 1), length(vals)));
for nB = 1:length(vals)
    bH = bar(nB, vals(nB), 0.4, 'facecolor', cmap(cmap_inds(nB), :));
    hold on;
end
set(gH, 'Xtick', []);

ylabel(ylabelstr, 'parent', gH);
xlabel(xlabelstr, 'parent', gH);
title(titlestr, 'parent', gH);
legend(legendstr, 'Location', 'bestoutside');


function dispParaStr = dispPara(fusionInfo)


maskType = 'Default';
if (~isempty(fusionInfo.setup_analysis.maskFile))
    maskType = 'User specified';
end

dispParaStr{1} = '....................................................';
dispParaStr{end + 1} = ['Feature names: ', ica_fuse_formatStr(fusionInfo.run_analysis.featureNames, ',')];
dispParaStr{end + 1} = ['Modalities: ', ica_fuse_formatStr(fusionInfo.run_analysis.modalities, ',')];
dispParaStr{end + 1} = ['Number of components: ', num2str(fusionInfo.run_analysis.numComp)];
dispParaStr{end + 1} = ['Group names: ', ica_fuse_formatStr(fusionInfo.run_analysis.groupNames, ',')];
dispParaStr{end + 1} = ['Number of subjects: ', num2str(fusionInfo.run_analysis.numSubjects)];
dispParaStr{end + 1} = ['Type Of PCA: ', fusionInfo.run_analysis.type_pca];
dispParaStr{end + 1} = ['Mask: ', maskType];
dispParaStr{end + 1} = ['Anatomical file: ', fusionInfo.run_analysis.anatomical_file];
dispParaStr{end + 1} = ['Slice Plane: ', upper(fusionInfo.run_analysis.anatomical_plane(1)), fusionInfo.run_analysis.anatomical_plane(2:end)];
dispParaStr{end + 1} = ['Image values: ', upper(fusionInfo.run_analysis.image_values(1)), fusionInfo.run_analysis.image_values(2:end)];
dispParaStr{end + 1} = ['Convert to Z-scores: ', fusionInfo.run_analysis.convert_to_z];
dispParaStr{end + 1} = ['Threshold: ', num2str( fusionInfo.run_analysis.threshold)];
dispParaStr{end + 1} = '....................................................';

function ica_fuse_pica_summary(param_file, displayParameters)
% Para ICA report generator

ica_fuse_defaults;
global ANATOMICAL_FILE;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGE_VALUES;
global ANATOMICAL_PLANE;
global PARA_FUSE_PARAM_FILE;
global FIG_FG_COLOR;


if (~exist('param_file', 'var'))
    param_file = PARA_FUSE_PARAM_FILE;
end


if (~isempty(param_file))
    load(param_file);
end

outputDir = fileparts(param_file);

if (isempty(outputDir))
    outputDir = pwd;
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

sz = get(0, 'ScreenSize');
defaultFigPos = get(0, 'DefaultFigurePosition');
defaultFigPos(3) = 0.9*sz(3);
defaultFigPos(2) = 0.5*sz(4) - 0.5*defaultFigPos(4) + 10;
defaultFigPos(1) = 0.5*sz(3) - 0.5*defaultFigPos(3) + 10;

% Get params from the output files
%slices_in_mm = (-40:4:72);
%threshold = threshold;
convert_to_z = convertToZStr;
imageValues = lower(imageValueStr);
%anatomical_plane = ANATOMICAL_PLANE;
anatomical_file = ANATOMICAL_FILE;
outputFiles = paraICAInfo.run_analysis.outputFiles;
linked_correlation = paraICAInfo.run_analysis.average_correlation;
corrIndices = paraICAInfo.run_analysis.corrIndices;
modalities = cellstr(char(paraICAInfo.run_analysis.dataInfo(1).feature.modality));
groupNames = cellstr(char(paraICAInfo.run_analysis.dataInfo.name));
featureNames = cellstr(char(paraICAInfo.run_analysis.dataInfo(1).feature.name));
numComp = paraICAInfo.run_analysis.numComp;
numSubjects = paraICAInfo.run_analysis.numSubjects;

paraICAInfo.run_analysis.modalities = modalities;
paraICAInfo.run_analysis.groupNames = groupNames;
paraICAInfo.run_analysis.featureNames = featureNames;
paraICAInfo.run_analysis.numComp = numComp;
paraICAInfo.run_analysis.numSubjects = numSubjects;
paraICAInfo.run_analysis.threshold = threshold;
paraICAInfo.run_analysis.convert_to_z = convert_to_z;
paraICAInfo.run_analysis.anatomical_plane = anatomical_plane;
paraICAInfo.run_analysis.anatomical_file = anatomical_file;
paraICAInfo.run_analysis.slices_in_mm =slices_in_mm;
paraICAInfo.run_analysis.image_values = imageValues;
dataInfo = paraICAInfo.run_analysis.dataInfo;
mask_ind = paraICAInfo.run_analysis.mask_ind;
[dims, voxels] = ica_fuse_getFeatureDIM(mask_ind);


numLoops = 1;
if (length(featureNames) == 2)
    [dd, sorted_inds] = sort(abs(linked_correlation), 'descend');
    linked_correlation = linked_correlation(sorted_inds);
    selected_comps{1} = sorted_inds;
    selected_comps{2} = corrIndices(sorted_inds);
    numLoops = length(sorted_inds);
else
    selected_comps = {corrIndices(1), corrIndices(2), corrIndices(3)};
end


feature_files = cell(1, length(featureNames));
loadingCoeff = cell(1, length(featureNames));
allColors = loadingCoeff;
groupLoadingCoeff = cell(length(groupNames), length(featureNames));
for nF = 1:length(featureNames)
    
    feature_files{nF} = ica_fuse_rename_4d_file(outputFiles(nF).name);
    loadingCoeff{nF} = ica_fuse_loadData(fullfile(outputDir, outputFiles(nF).loadingCoeffFile));
    allColors{nF} = colordg(nF);
    
    eG = 0;
    for nG = 1:length(groupNames)
        sG = eG + 1;
        eG = eG + numSubjects(nG);
        groupLoadingCoeff{nG, nF} = loadingCoeff{nF}(sG:eG, :);
    end
    
end


%% Parallel ICA Fusion Parameters
dispParaStr = dispPara(paraICAInfo);

fprintf('\n\n');
disp(char(dispParaStr));
fprintf('\n\n');


%%

%% Correlations between the linked components
% Component number is displayed for each feature and the linked correlation
% is reported.
fNs = featureNames;
fNs = strrep(fNs, ' ', '');
if (length(fNs) == 2)
    varNames = {fNs{1}, fNs{2}, 'Correlation', 'Tvalue', 'Pvalue'};
    [tstat, pval] = r_to_t(linked_correlation(:), sum(numSubjects));
    printToTable({selected_comps{1}(:), selected_comps{2}(:), linked_correlation(:), tstat(:), pval(:)}, varNames, {'%d', '%d', '%0.3f', '%0.4f', '%0.4f'});
else
    % 1 vs 2
    varNames = {fNs{1}, fNs{2}, 'Correlation', 'Tvalue', 'Pvalue'};
    [tstat, pval] = r_to_t(linked_correlation(:), sum(numSubjects));
    printToTable({corrIndices(1), corrIndices(2), linked_correlation(1), tstat(1), pval(1)}, varNames, {'%d', '%d', '%0.3f', '%0.4f', '%0.4f'});
    % 2 vs 3
    varNames = {fNs{2}, fNs{3}, 'Correlation', 'Tvalue', 'Pvalue'};
    printToTable({corrIndices(2), corrIndices(3), linked_correlation(2), tstat(2), pval(2)}, varNames, {'%d', '%d', '%0.3f', '%0.4f', '%0.4f'});
    % 1 vs 3
    varNames = {fNs{1}, fNs{3}, 'Correlation', 'Tvalue', 'Pvalue'};
    printToTable({corrIndices(1), corrIndices(3), linked_correlation(3), tstat(3), pval(3)}, varNames, {'%d', '%d', '%0.3f', '%0.4f', '%0.4f'});
end

clear fNs;

%% Parallel ICA components
% Mixing coefficients and the source maps are displayed. Only the linked
% components are shown. When plotting the component EEG signal, mean of original
% data is also plotted.
%
mixingStruct.numSubjects =  numSubjects;
mixingStruct.groupNames = char(groupNames);


meanEEG = cell(1, length(featureNames));
for nFea = 1:length(featureNames)
    
    if (strcmpi(modalities{nFea}, 'eeg'))
        for nInputFiles = 1:length(dataInfo)
            if (nInputFiles == 1)
                feature_input_files = char(dataInfo(nInputFiles).feature(nFea).files.name);
            else
                feature_input_files = char(feature_input_files, ...
                    char(dataInfo(nInputFiles).feature(nFea).files.name));
            end
        end
        input_data = ica_fuse_loadData(feature_input_files);
        meanData = mean(squeeze(input_data(paraICAInfo.run_analysis.mask_ind(nFea).ind, 2, :)), 2);
        meanEEG{nFea} = meanData;
    end
end


for nComp = 1:numLoops
    
    axesTitle = 'Mixing Coefficients ';
    %mixingStruct.axesTitle = axesTitle;
    compStr = [];
    for nF = 1:length(featureNames)
        tmpL = loadingCoeff{nF}(:, selected_comps{nF}(nComp));
        compStr = [compStr, featureNames{nF}, ' Comp ', ica_fuse_returnFileIndex(selected_comps{nF}(nComp)), ' '];
        mixingStruct.feature(nF).loadingCoeff = tmpL/norm(tmpL, 2);
    end
    mixingStruct.axesTitle = [axesTitle, compStr];
    gH = ica_fuse_getGraphics(axesTitle, 'timecourse', '', 'on');
    set(gH, 'resize', 'on');
    
    ica_fuse_plotParaICALoading('data', mixingStruct, 'parent', gca, 'num_features', length(featureNames), 'colors', allColors, 'legend_string', featureNames, 'titlecolor', FIG_FG_COLOR);
    xlabel('Subjects', 'parent', gca);
    ylabel('Normalized Units', 'parent', gca);
    
    for nFea = 1:length(selected_comps)
        
        comp1 = selected_comps{nFea};
        modality1_comp = fullfile(outputDir, deblank(feature_files{nFea}(comp1(nComp), :)));
        title_str = [featureNames{nFea},  ' Comp ', ica_fuse_returnFileIndex(nComp)];
        if (strcmpi(modalities{nFea}, 'smri') || strcmpi(modalities{nFea}, 'fmri'))
            % fMRI or sMRI
            gH = ica_fuse_getGraphics(title_str, 'graphics', '', 'on');
            set(gH, 'resize', 'on');
            ica_fuse_image_viewer(modality1_comp, 'structfile', anatomical_file, 'threshold', threshold, 'convert_to_zscores', convert_to_z, 'image_values', imageValues, ...
                'slices_in_mm', slices_in_mm, 'anatomical_view', anatomical_plane, 'labels', [featureNames{nFea},  ' Comp ', ica_fuse_returnFileIndex(comp1(nComp))], 'axesh', gca);
        elseif strcmpi(modalities{nFea}, 'eeg')
            % EEG
            %title_str = [featureNames{nFea},  ' Comp ', ica_fuse_returnFileIndex(nComp)];
            tmp_dat = ica_fuse_loadData(modality1_comp);
            tmpMeanData = meanEEG{nFea};
            interpFactor = ceil(voxels / length(tmpMeanData));
            if (interpFactor > 0)
                % resample mean data
                tmpMeanData = ica_fuse_resample(tmpMeanData, voxels, length(tmpMeanData));
                % resample component data
                tmp_dat = ica_fuse_resampleCompData(tmp_dat, length(tmpMeanData), voxels, size(tmp_dat, 1));
            end
            
            gH = ica_fuse_getGraphics(title_str, 'timecourse', '', 'on');
            set(gH, 'resize', 'on');
            ica_fuse_plotTimecourse('parent', gca, 'data', tmp_dat, 'color', 'c', 'titleColor', FIG_FG_COLOR, 'title', title_str, 'YAxisLocation', 'left', 'meanData', ...
                tmpMeanData, 'meanDataLegend', {'Mean'},  'groupCompData',[], 'groupCompLegend', {});
            legend('Mean', 'Component', 'Location', 'best');
        else
            % Gene
            title_str = [featureNames{nFea},  ' Comp ', ica_fuse_returnFileIndex(nComp)];
            tmp_dat = ica_fuse_loadData(modality1_comp);
            tmp_dat = squeeze(tmp_dat(:, 2, :));
            gH = ica_fuse_getGraphics(title_str, 'timecourse', '', 'on');
            set(gH, 'resize', 'on');
            ica_fuse_plotSNP('data',tmp_dat, 'parent', gca, 'xlabel', 'SNPs', 'ylabel', ...
                'Z-scores', 'titleColor', FIG_FG_COLOR, 'color', 'm', 'title', [featureNames{nFea},  ' Comp ', ica_fuse_returnFileIndex(comp1(nComp))]);
        end
        
    end
    
end


%% Two sample t-test results
% Two sample t-test is computed on the mixing coefficents of each feature between the selected groups.
%

if (length(groupNames)  == 1)
    return;
end

allGroups = nchoosek(1:length(groupNames), 2);
sortR = cell(size(allGroups, 1), length(featureNames));
for n = 1:size(allGroups, 1)
    % Current group
    cG = allGroups(n, :);
    
    for nF = 1:length(featureNames)
        
        g1 = groupLoadingCoeff{cG(1), nF};
        g2 = groupLoadingCoeff{cG(2), nF};
        
        tvals = zeros(size(g1, 2), 1);
        pvals = tvals;
        %legendStr = cell(size(g1, 2), 1);
        for nComp = 1:size(g1, 2)
            [pvals(nComp), tvals(nComp)] = ica_fuse_ttest2(g1(:, nComp), g2(:, nComp), 0);
        end
        
        dispStr = [deblank(groupNames{cG(1)}), ' vs ', deblank(groupNames{cG(2)}), ' Of feature ', featureNames{nF}];
        
        legendStr = cellstr(strcat('Comp ', num2str((1:size(g1, 2))')));
        
        components = (1:length(tvals))';
        tmpR.components = components;
        tmpR.pvals = pvals;
        tmpR.tvals = tvals;
        tmpR.title = dispStr;
        tmpR.legend = legendStr;
        
        sortR{n, nF} = tmpR;
        
        fH = figure('color', 'w', 'position', defaultFigPos);
        gH = axes('parent', fH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
        cmap = hsv(64);
        
        cmap_inds = round(linspace(1, size(cmap, 1), length(tvals)));
        for nB = 1:length(tvals)
            bH = bar(nB, tvals(nB), 0.4, 'facecolor', cmap(cmap_inds(nB), :));
            hold on;
        end
        set(gH, 'Xtick', []);
        
        ylabel('T-values', 'parent', gH);
        xlabel('Components', 'parent', gH);
        title(tmpR.title, 'parent', gH);
        legend(legendStr, 'Location', 'bestoutside');
        
        
        
        disp(dispStr);
        % plot_bars(gH, tvals, 'Components', 'T-values', dispStr, legendStr);
        
        set(gcf, 'name', dispStr);
        set(gcf, 'tag', 'tvals');
        
        varNames = {'ComponentNumber', 'Pvalues', 'Tvalues'};
        
        printToTable({components, pvals(:), tvals(:)}, varNames, {'%d', '%d', '%0.3f'});
        
    end
    
end



function plot_bars(gH, vals, xlabelstr, ylabelstr, titlestr, legendstr)


%gH = axes('parent', fH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
cmap = hsv(64);

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



function linecolor = colordg(n)
%COLORDG - Provides a choice of 15 colors for a line plot
%The first seven colors are the same as Matlab's predefined
%values for the PLOT command, i.e.
%
% 'b','g','r','c','m','y','k'
%
%Syntax: linecolor = colordg(n);
%
%Input: N , integer between 1 and 15, giving the following colors
%
% 1 BLUE
% 2 GREEN (pale)
% 3 RED
% 4 CYAN
% 5 MAGENTA (pale)
% 6 YELLOW (pale)
% 7 BLACK
% 8 TURQUOISE
% 9 GREEN (dark)
% 10 YELLOW (dark)
% 11 ORANGE
% 12 MAGENTA (dark)
% 13 GREY
% 14 BROWN (pale)
% 15 BROWN (dark)
%
%Output: LINECOLOR (1 x 3 RGB vector)
%
%Examples:
% 1) h = line(x,y,'Color',colordg(11)); %Picks the orange color
% 2) colordg demo %Creates a figure displaying the 15 colors
% 3) axes; set(gca,'ColorOrder',(colordg(1:15)));
% Overrides the default ColorOrder for the current axes only
% 4) figure; set(gcf,'DefaultAxesColorOrder',(colordg(1:15)));
% Overrides the default ColorOrder for all axes of the current
% figure
% 5) set(0,'DefaultAxesColorOrder',(colordg(1:15)));
% Sets the default ColorOrder for all axes to be created during
% the current matlab session. You may wish to insert this
% command into your startup.m file.
%
%See also: PLOT, LINE, AXES

%Author: Denis Gilbert, Ph.D., physical oceanography
%Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%Web: http://www.qc.dfo-mpo.gc.ca/iml/
%August 2000; Last revision: 26-Sep-2003

if nargin == 0
    error('Must provide an input argument to COLORDG')
end

colorOrder = ...
    [ 0 0 1 % 1 BLUE
    0 1 0 % 2 GREEN (pale)
    1 0 0 % 3 RED
    0 1 1 % 4 CYAN
    1 0 1 % 5 MAGENTA (pale)
    1 1 0 % 6 YELLOW (pale)
    0 0 0 % 7 BLACK
    0 0.75 0.75 % 8 TURQUOISE
    0 0.5 0 % 9 GREEN (dark)
    0.75 0.75 0 % 10 YELLOW (dark)
    1 0.50 0.25 % 11 ORANGE
    0.75 0 0.75 % 12 MAGENTA (dark)
    0.7 0.7 0.7 % 13 GREY
    0.8 0.7 0.6 % 14 BROWN (pale)
    0.6 0.5 0.4 ]; % 15 BROWN (dark)


if isnumeric(n) & n >= 1 & n <= 15
    linecolor = colorOrder(n,:);
elseif strcmp(n,'demo')
    %GENERATE PLOT to display a sample of the line colors
    figure, axes;
    %PLOT N horizontal lines
    for n=1:length(colorOrder)
        h(n) = line([0 1],[n n],'Color',colorOrder(n,:));
    end
    set(h,'LineWidth',5)
    set(gca,'YLim',[0 n+1],'YTick',[1:n],'XTick',[])
    ylabel('Color Number');
else
    error('Invalid input to colordg');
end


function [tstat, pval] = r_to_t(r, n)
% correlation value to t stat

tstat = r./sqrt((1 - r.^2)./(n - 2));
pval = ica_fuse_get_pvalue(tstat, n - 2, 0);


function printToTable(vars, varNames, varStr)

if (~isempty(which('table.m')))
    T = table(vars{:}, 'VariableNames', varNames);
    disp(T);
else
    formatStr = repmat('%20s\t', 1, length(vars));
    formatStr(end-1:end) = '';
    formatStr = [formatStr, '\n'];
    fprintf(formatStr, varNames{:});
    fprintf('\n');
    for n = 1:length(vars{1})
        formVars = cell(1, length(vars));
        for nVar = 1:length(formVars)
            formVars{nVar} = num2str(vars{nVar}(n), varStr{nVar});
        end
        fprintf(formatStr, formVars{:});
    end
    fprintf('\n');
end


function dispParaStr = dispPara(paraICAInfo)


maskType = 'Default';
if (~isempty(paraICAInfo.setup_analysis.maskFile))
    maskType = 'User specified';
end

dispParaStr{1} = '....................................................';
dispParaStr{end + 1} = ['Feature names: ', ica_fuse_formatStr(paraICAInfo.run_analysis.featureNames, ',')];
dispParaStr{end + 1} = ['Modalities: ', ica_fuse_formatStr(paraICAInfo.run_analysis.modalities, ',')];
dispParaStr{end + 1} = ['Number of components: ', num2str(paraICAInfo.run_analysis.numComp)];
dispParaStr{end + 1} = ['Group names: ', ica_fuse_formatStr(paraICAInfo.run_analysis.groupNames, ',')];
dispParaStr{end + 1} = ['Number of subjects: ', num2str(paraICAInfo.run_analysis.numSubjects)];
dispParaStr{end + 1} = ['Type Of Parallel ICA: ', paraICAInfo.run_analysis.type_parallel_ica];
dispParaStr{end + 1} = ['Type Of PCA: ', paraICAInfo.run_analysis.type_pca];
dispParaStr{end + 1} = ['Mask: ', maskType];
dispParaStr{end + 1} = ['Anatomical file: ', paraICAInfo.run_analysis.anatomical_file];
dispParaStr{end + 1} = ['Slice Plane: ', upper(paraICAInfo.run_analysis.anatomical_plane(1)), paraICAInfo.run_analysis.anatomical_plane(2:end)];
dispParaStr{end + 1} = ['Image values: ', upper(paraICAInfo.run_analysis.image_values(1)), paraICAInfo.run_analysis.image_values(2:end)];
dispParaStr{end + 1} = ['Convert to Z-scores: ', paraICAInfo.run_analysis.convert_to_z];
dispParaStr{end + 1} = ['Threshold: ', num2str( paraICAInfo.run_analysis.threshold)];
dispParaStr{end + 1} = '....................................................';

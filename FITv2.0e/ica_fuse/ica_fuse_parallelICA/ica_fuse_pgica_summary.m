function ica_fuse_pgica_summary(param_file, opts)
% Para ICA report generator

ica_fuse_defaults;
global ANATOMICAL_FILE;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGE_VALUES;
global ANATOMICAL_PLANE;
global FIG_FG_COLOR;


load(param_file);

outputDir = fileparts(param_file);

if (isempty(outputDir))
    outputDir = pwd;
end


linked_correlation = pgicaInfo.corr_vals;
corrIndices = pgicaInfo.corr_inds;
loadingFiles = pgicaInfo.loadingFiles;
numSubjects = pgicaInfo.num_subjects;
outputFiles = pgicaInfo.outputFiles;

anatomical_file = ANATOMICAL_FILE;

try
    imageValues = opts.display.image_values;
catch
    imageValues = IMAGE_VALUES;
end

try
    threshold = abs(opts.display.z_threshold);
catch
    threshold = Z_THRESHOLD;
end

try
    anatomical_plane = lower(opts.display.anatomical_plane);
catch
    anatomical_plane = ANATOMICAL_PLANE;
end
slices_in_mm=[];
try
    slices_in_mm = opts.display.slices_in_mm;
catch
    if ~isempty(anatomical_file)
        params = ica_fuse_get_slice_def(ica_fuse_spm_vol(anatomical_file), anatomical_plane);
        slices_in_mm = params.slices;
    end
    %slices_in_mm = -40:4:72;
end

try
    convertToZ = opts.display.convert_to_z;
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


pgicaInfo.anatomical_file = anatomical_file;
pgicaInfo.anatomical_plane = anatomical_plane;
pgicaInfo.threshold = threshold;
pgicaInfo.convert_to_z = convertToZStr;
pgicaInfo.slices_in_mm = slices_in_mm;
pgicaInfo.image_values = imageValueStr;

groupsInfo = [];
try
    groupsInfo = opts.stats.groupsInfo;
catch
end


[dd, sorted_inds] = sort(abs(linked_correlation), 'descend');
linked_correlation = linked_correlation(sorted_inds);
selected_comps{1} = sorted_inds;
selected_comps{2} = corrIndices(sorted_inds);



sz = get(0, 'ScreenSize');
defaultFigPos = get(0, 'DefaultFigurePosition');
defaultFigPos(3) = 0.9*sz(3);
defaultFigPos(2) = 0.5*sz(4) - 0.5*defaultFigPos(4) + 10;
defaultFigPos(1) = 0.5*sz(3) - 0.5*defaultFigPos(3) + 10;


featureNames = {'sMRI', 'fMRI'};
numFeatures = length(loadingFiles);
loadingCoeff{1} = ica_fuse_loadData(fullfile(outputDir, loadingFiles{1}));
loadingCoeff{2} = ica_fuse_loadData(fullfile(outputDir, loadingFiles{2}));
allColors = cell(1, length(loadingCoeff));

groupNames = {''};


if (~isempty(groupsInfo))
    
    groupNames = cellstr(char(groupsInfo.name));
    groupLoadingCoeff = cell(length(groupNames), length(featureNames));
    for nF = 1:numFeatures
        
        allColors{nF} = colordg(nF);
        
        for nG = 1:length(groupsInfo)
            groupLoadingCoeff{nG, nF} = loadingCoeff{nF}(groupsInfo(nG).value, :);
        end
    end
    
end

pgicaInfo.featureNames = featureNames;

%% PGICA-ICA Fusion Parameters
dispParaStr = dispPara(pgicaInfo);

fprintf('\n\n');
disp(char(dispParaStr));
fprintf('\n\n');


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


clear nG nF;

%% SMRI Components
for nFeature = 1:1
    
    axesTitle = ['Mixing Coefficients (', featureNames{nFeature}];
    
    for nComp = 1:length(outputFiles{nFeature})
        
        
        mixingStruct.axesTitle = [axesTitle, ' Comp ', icatb_returnFileIndex(nComp), ')'];
        tmpL = loadingCoeff{nFeature}(:, nComp);
        mixingStruct.feature(1).loadingCoeff = tmpL/norm(tmpL, 2);
        mixingStruct.numSubjects = numSubjects;
        mixingStruct.groupNames = {''};
        
        gH = ica_fuse_getGraphics(axesTitle, 'timecourse', '', 'on');
        set(gH, 'resize', 'on');
        
        ica_fuse_plotParaICALoading('data', mixingStruct, 'parent', gca, 'num_features', 1, 'colors', allColors, 'legend_string', featureNames{nFeature}, 'titlecolor', FIG_FG_COLOR);
        xlabel('Subjects', 'parent', gca);
        ylabel('Normalized Units', 'parent', gca);
        
        title_str = [featureNames{nFeature},  ' Comp ', ica_fuse_returnFileIndex(nComp)];
        % fMRI or sMRI
        gH = ica_fuse_getGraphics(title_str, 'graphics', '', 'on');
        set(gH, 'resize', 'on');
        
        if (isempty(anatomical_file))
            anat_file = deblank(pgicaInfo.files{nFeature}(1, :));
            tmp_slices = [];
        else
            anat_file = anatomical_file;
            tmp_slices = slices_in_mm;
        end
        
        ica_fuse_image_viewer(fullfile(outputDir, outputFiles{nFeature}{nComp}), 'structfile', anat_file, 'threshold', threshold, 'convert_to_zscores', convertToZ, 'image_values', imageValues, ...
            'slices_in_mm', tmp_slices, 'anatomical_view', anatomical_plane, 'labels', [featureNames{nFeature},  ' Comp ', ica_fuse_returnFileIndex(nComp)], 'axesh', gca);
        
    end
    
end


drawnow;


%% fMRI Components
for nFeature = 2
    
    axesTitle = ['Mixing Coefficients (', featureNames{nFeature}];
    
    for nComp = 1:length(outputFiles{nFeature})
        
        
        mixingStruct.axesTitle = [axesTitle, ' Comp ', icatb_returnFileIndex(nComp), ')'];
        tmpL = loadingCoeff{nFeature}(:, nComp);
        mixingStruct.feature(1).loadingCoeff = tmpL/norm(tmpL, 2);
        mixingStruct.numSubjects = numSubjects;
        mixingStruct.groupNames = {''};
        
        gH = ica_fuse_getGraphics(axesTitle, 'timecourse', '', 'on');
        set(gH, 'resize', 'on');
        
        ica_fuse_plotParaICALoading('data', mixingStruct, 'parent', gca, 'num_features', 1, 'colors', allColors, 'legend_string', featureNames{nFeature}, 'titlecolor', FIG_FG_COLOR);
        xlabel('Subjects', 'parent', gca);
        ylabel('Normalized Units', 'parent', gca);
        
        title_str = [featureNames{nFeature},  ' Comp ', ica_fuse_returnFileIndex(nComp)];
        % fMRI or sMRI
        gH = ica_fuse_getGraphics(title_str, 'graphics', '', 'on');
        set(gH, 'resize', 'on');
        
        if (isempty(anatomical_file))
            anat_file = deblank(pgicaInfo.files{nFeature}(1, :));
            tmp_slices = [];
        else
            anat_file = anatomical_file;
            tmp_slices = slices_in_mm;
        end
        
        ica_fuse_image_viewer(fullfile(outputDir, outputFiles{nFeature}{nComp}), 'structfile', anat_file, 'threshold', threshold, 'convert_to_zscores', convertToZ, 'image_values', imageValues, ...
            'slices_in_mm', tmp_slices, 'anatomical_view', anatomical_plane, 'labels', [featureNames{nFeature},  ' Comp ', ica_fuse_returnFileIndex(nComp)], 'axesh', gca);
        
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



function dispParaStr = dispPara(pgicaInfo)

maskFile_modality1 = pgicaInfo.maskFile_modality1;
if (isempty(maskFile_modality1))
    maskFile_modality1 = 'Default';
end


maskFile_modality2 = pgicaInfo.maskFile_modality2;
if (isempty(maskFile_modality2))
    maskFile_modality2 = 'Default';
end

dispParaStr{1} = '....................................................';
dispParaStr{end + 1} = ['Modalities: ', ica_fuse_formatStr(pgicaInfo.featureNames, ',')];
dispParaStr{end + 1} = ['Number of sMRI components: ', num2str(pgicaInfo.numComp1)];
dispParaStr{end + 1} = ['Number of fMRI components (first level, second level): (', num2str(pgicaInfo.numComp2(1)), ', ', num2str(pgicaInfo.numComp2(2)), ')'];
dispParaStr{end + 1} = ['Number of subjects: ', num2str(pgicaInfo.num_subjects)];
dispParaStr{end + 1} = ['Masks (sMRI and fMRI): ', maskFile_modality1, ', ', maskFile_modality2];
dispParaStr{end + 1} = ['Anatomical file: ', pgicaInfo.anatomical_file];
dispParaStr{end + 1} = ['Slice Plane: ', upper(pgicaInfo.anatomical_plane(1)), pgicaInfo.anatomical_plane(2:end)];
dispParaStr{end + 1} = ['Image values: ', upper(pgicaInfo.image_values(1)), pgicaInfo.image_values(2:end)];
dispParaStr{end + 1} = ['Convert to Z-scores: ', pgicaInfo.convert_to_z];
dispParaStr{end + 1} = ['Threshold: ', num2str(pgicaInfo.threshold)];
dispParaStr{end + 1} = '....................................................';
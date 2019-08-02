function ica_fuse_display_histogram(histParameters, selComp, selFeatures, threshold, histogram_criteria)
%% Display histograms
%
%

ica_fuse_defaults;
global UI_FONT_NAME;
global COLORMAPFILE;
global FIG_FG_COLOR;

if ~exist('histParameters', 'var')
    error('Parameters structure required for generating histograms is not provided');
end

load(COLORMAPFILE);
fontSize = 0.05;
groupNames = histParameters.groupNames;
selGroups = histParameters.selGroupsVal;
dataInfo = histParameters.dataInfo;
allFeatureNames = str2mat(dataInfo(1).feature.name);

% Use atmost two features for plotting
if (length(selFeatures) > 2)
    selFeatures = selFeatures(1:2);
end

% Check if the selected features exist or not
if isfield(histParameters, 'selectedFeatures')
    featureNames = histParameters.selectedFeatures;
else
    % Treating this as the first MAT file number which contains all
    % combinations
    featureNames = deblank(allFeatureNames(selFeatures, :));
end


if (isfield(histParameters, 'selectedCombNumber'))
    combNumber = histParameters.selectedCombNumber;
else
    combNumber = 1;
end

histData = ica_fuse_calculate_histogram(histParameters, histogram_criteria, threshold, combNumber, selFeatures, selComp);

graphicsH = repmat(struct('H', []), 1, length(histData(1).group));
%%%%%%%% Display histograms %%%%%%%
% Plot Histogram of groups
for nn = 1:length(histData(1).group)
    hData = histData(1).group(nn).data;
    hData = log(hData + 1);
    size_data = size(hData);
    if nn == 1
        meanHist = zeros(size(hData, 1), size(hData, 2), length(selGroups));
        voxelValues = histData(1).hIndex;
    end
    meanHist(:, :, nn) = mean(hData, length(size_data));
    [graphicsH(nn).H, axesHandle] = plotIm(hData, ['Group ', deblank(groupNames(selGroups(nn), :))], hot);
    axesH = findobj(graphicsH(nn).H, 'tag', 'axes1');
    drawText(axesH, selFeatures, featureNames);
end

% Plot difference histograms
% Plot First group followed by second group
m1 = squeeze(meanHist(:, :, 1));

% Mean for group 1
graphCount = length(graphicsH);
axesTitle = ['Mean of Group ', deblank(groupNames(selGroups(1), :))];
[graphicsH(graphCount + 1).H, axesHandle] = plotIm(m1, axesTitle, hot, axesTitle);
drawText(axesHandle, selFeatures, featureNames);
drawTicks(axesHandle, voxelValues);

graphCount = length(graphicsH);


if (length(selGroups) == 2)
    
    m2 = squeeze(meanHist(:, :, 2));
    
    graphCount = length(graphicsH);
    axesTitle = ['Mean of Group ', deblank(groupNames(selGroups(2), :))];
    % Mean for group 2
    [graphicsH(graphCount + 1).H, axesHandle] = plotIm(m2, axesTitle, hot, axesTitle);
    drawText(axesHandle, selFeatures, featureNames);
    drawTicks(axesHandle, voxelValues);
    
    graphCount = length(graphicsH);
    
    out = m1 - m2;
    group1Name = deblank(groupNames(selGroups(1), :));
    group2Name = deblank(groupNames(selGroups(2), :));
    
    graphCount = length(graphicsH);
    axesTitle = ['Group ', group1Name , ' - Group ', group2Name];
    [graphicsH(graphCount + 1).H, axesHandle] = plotIm(out, ['Group ', group1Name , ' - Group ', group2Name], coldhot, axesTitle);
    CLIM = get(axesHandle, 'CLIM');
    CLIM = [-max(abs(CLIM)), max(abs(CLIM))];
    set(axesHandle, 'CLIM', CLIM);
    drawText(axesHandle, selFeatures, featureNames);
    drawTicks(axesHandle, voxelValues);
    
    % Number of graphs
    graphCount = length(graphicsH);
    
    axesHandle = get(graphicsH(graphCount).H, 'currentAxes');
    axesPos = get(axesHandle, 'position');
    set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
    colorbarPos = [axesPos(1) + axesPos(3) + 0.03, axesPos(2), 0.05, axesPos(4)];
    if ((ica_fuse_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
        colorbarHandle = colorbar('peer', axesHandle);
    else
        colorbarHandle = colorbar;
    end
    
    set(colorbarHandle, 'position', colorbarPos);
    set(colorbarHandle, 'units', 'pixels');
    
    if ((ica_fuse_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
        % Maximum
        xPos = 0.1; yPos = 1.05;
        text(xPos, yPos, group1Name, 'units', 'normalized', 'parent', colorbarHandle, 'color', FIG_FG_COLOR,  ...
            'fontsize', fontSize, 'HorizontalAlignment', 'left', ...
            'FontName', UI_FONT_NAME, 'fontunits', 'normalized');
        
        % Minimum
        yPos = -0.05;
        text(xPos, yPos, group2Name, 'units', 'normalized', 'parent', colorbarHandle, 'color', FIG_FG_COLOR,  ...
            'fontsize', fontSize, 'HorizontalAlignment', 'left', 'FontName', UI_FONT_NAME, 'fontunits', ...
            'normalized');
    else
        labelsH = get(colorbarHandle, 'label');
        set(colorbarHandle, 'color', FIG_FG_COLOR);
        set(labelsH, 'units', 'normalized');
        title(group1Name, 'parent', colorbarHandle, 'color', FIG_FG_COLOR, 'Horizontalalignment', 'center');
        set(labelsH, 'string', group2Name, 'rotation',0);
        set(labelsH, 'position', [0.45,-0.05, 0]);
        set(colorbarHandle, 'YTickLabel', []);
    end
    set(colorbarHandle, 'YTickLabel', []);
    set(colorbarHandle, 'YTick', []);
    
end

drawnow;

graphicsH = plotMarginals(histData, selFeatures, featureNames, selGroups, groupNames, graphicsH);

% Plot Next previous and exit buttons
ica_fuse_plotNextPreviousExitButtons(graphicsH);

disp('Done calculating histograms');


function [figHandle, axesHandle] = plotIm(data, titleStr, cmap, axesTitle)
% Sub function to plot images in a sub plot

ica_fuse_defaults;
global AX_COLOR; % axes color
global LEGEND_COLOR; % LEGEND_COLOR for Matlab 7.0
global FIG_FG_COLOR;
global UI_FONT_NAME;
global UI_FONT_UNITS;
global UI_FONT_SIZE;

figHandle = ica_fuse_getGraphics(titleStr, 'normal', titleStr, 'on');

if length(size(data)) == 3
    numPlots = size(data, 3);
else
    numPlots = 1;
end

data = data(:, end:-1:1, :);
data = permute(data, [2 1 3]);

userdata = zeros(size(data));

for nn = 1:numPlots
    userdata(:, :, nn) = (squeeze(data(:, :, nn)));
end

% Display GUI Options Menu
figOptions = uimenu('parent', figHandle, 'tag', 'extract_data', 'label', 'Extract Data', 'callback', ...
    {@extractDataCallback, figHandle}, 'userdata', userdata);

set(figHandle, 'colormap', cmap);

numRows = ceil(sqrt(numPlots));
numCols = ceil(numPlots/numRows);

fontSize = 0.05;

if numRows*numCols > 1
    fontSize = 0.1;
end

countN = 0;
for nRow = 1:numRows
    
    for nCol = 1:numCols
        countN = countN + 1;
        
        if countN <= numPlots
            if numPlots > 1
                figure(figHandle);
                axesHandle = subplot(numRows, numCols, countN);
            else
                axesHandle = axes('parent', figHandle, 'units', 'normalized', 'position', [0.12 0.12 0.7 0.7]);
            end
            
            set(axesHandle, 'units', 'normalized');
            set(axesHandle, 'color', AX_COLOR);
            set(axesHandle,  'fontunits', 'normalized');
            set(axesHandle, 'FontName', UI_FONT_NAME, 'fontsize', fontSize);
            set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
            set(axesHandle, 'YTickLabel', []);
            set(axesHandle, 'YTick', []);
            set(axesHandle, 'XTickLabel', []);
            set(axesHandle, 'XTick', []);
            subject_str = ['Sub ', num2str(countN)];
            temp = squeeze(data(:, :, countN));
            %temp = temp';
            ImageAxis = image(temp, 'parent', axesHandle, 'CDataMapping', 'scaled');
            
            set(axesHandle, 'clim', [min(data(:)), max(data(:))]);
            clear temp;
            axis(axesHandle, 'off');
            if exist('axesTitle', 'var') && ~isempty(axesTitle)
                title(axesTitle, 'parent', axesHandle);
            else
                title(['Sub ', num2str(countN)], 'parent', axesHandle);
            end
            set(axesHandle, 'tag', ['axes', num2str(countN)]);
        else
            break;
        end
        
    end
    
end


function extractDataCallback(hObject, event_data, handles)
% Extract histogram data

data = get(hObject, 'userdata');

% open input dialog box
prompt = {'Enter a valid file name'};
dlg_title = 'Save Figure Data As';
num_lines = 1;
def = {'histogram_data'};
% save the file with the file name specified
fileName = ica_fuse_inputdlg2(prompt, dlg_title, num_lines, def);

if ~isempty(fileName)
    [dd, fileN, extn] = fileparts(fileName{1});
    if ~isempty(dd)
        outputDir = dd;
    else
        outputDir = pwd;
    end
    clear fileName;
    fileName = [fileN, '.mat'];
    fileName = fullfile(outputDir, fileName);
    ica_fuse_save(fileName, 'data');
    disp(['Figure data is saved in file: ', fileName]);
end


function drawText(axesH, selFeatures, featureNames)

ica_fuse_defaults;
global FIG_FG_COLOR;
global UI_FONT_NAME;

fontSize = 8;
axes(axesH);
axis(axesH, 'on');
set(axesH, 'box', 'off', 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
xlabel([deblank(featureNames(selFeatures(1), :))], 'parent', axesH, 'fontsize', fontSize);
if length(selFeatures) > 1
    ylabel([deblank(featureNames(selFeatures(2), :))], 'parent', axesH, 'fontsize', fontSize);
end
%set(axesH,  'fontunits', 'normalized');
%set(axesH, 'FontName', UI_FONT_NAME, 'fontsize', fontSize);
set(axesH, 'YTick', [], 'XTick', [], 'XTickLabel', [], 'YTickLabel', []);

function drawTicks(axesHandle, voxelValues)
% Draw X Ticks and Y ticks

ica_fuse_defaults;
global FIG_FG_COLOR;
global UI_FONT_NAME;

fontSize = 0.05;

axis(axesHandle, 'on');

set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
set(axesHandle,  'fontunits', 'normalized');
set(axesHandle, 'FontName', UI_FONT_NAME, 'fontsize', fontSize);
set(axesHandle, 'box', 'off', 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);

% Initialise ticks to empty
set(axesHandle, 'XTick', []);
set(axesHandle, 'XTickLabel', []);
set(axesHandle, 'YTick', []);
set(axesHandle, 'YTickLabel', []);

tickStr = ica_fuse_constructString(voxelValues{1});
v1 = eval(tickStr);

v1 = v1(:);

if length(v1) > 10
    ticks = 0:10:length(v1);
    ticks(1) = 1;
else
    ticks = (1:length(v1));
end

set(axesHandle, 'XTick', ticks);
set(axesHandle, 'XTickLabel', num2str(round(v1(ticks))));

if length(voxelValues) > 1
    voxelValues{2} = voxelValues{2}(end:-1:1);
    tickStr = ica_fuse_constructString(voxelValues{2});
    v2 = eval(tickStr);
    v2 = v2(:);
    if length(v2) > 10
        ticks = 0:10:length(v2);
        ticks(1) = 1;
    else
        ticks = (1:length(v2));
    end
    set(axesHandle, 'YTick', ticks);
    set(axesHandle, 'YTickLabel', num2str(round(v2(ticks))));
end



function graphicsH = plotMarginals(histData, selFeatures, featureNames, selGroups, groupNames, graphicsH)
% Plot marginals

if length(selFeatures) == 1
    return;
end

ica_fuse_defaults;
global AX_COLOR;
global UI_FONT_NAME;
global FIG_FG_COLOR;

% Number of groups
numGroups = length(histData(1).group);

groupNames = cellstr(groupNames(selGroups, :));


% Group 1 data
group1Data = histData(1).group(1).data;

data{1} = sum(mean(group1Data, length(size(group1Data)))');

if numGroups > 1
    % Group 2 data
    group2Data = histData(1).group(2).data;
    data{2} = sum(mean(group2Data, length(size(group2Data)))');
end


%%%%%%%%%% Plot Marginals of Feature 1 %%%%%%%%%%%%%%%%%%%%%%
featureName = deblank(featureNames(selFeatures(1), :));

% Edges
xind = histData(1).hIndex{1};

[graphicsH] = plot_marginal_feature(data, xind, featureName, groupNames, graphicsH);

%%%%%%%% End for plotting marginals of feature 1 %%%%%%%%%%%%%


%%%%%%%%%% Plot Marginals of Feature 2 %%%%%%%%%%%%%%%%%%%%%%
if length(histData(1).hIndex) > 1
    clear data;
    featureName = deblank(featureNames(selFeatures(2), :));
    % Edges
    yind = histData(1).hIndex{2};
    data{1} = sum(mean(group1Data, length(size(group1Data))));
    
    if numGroups > 1
        data{2} = sum(mean(group2Data, length(size(group2Data))));
    end
    [graphicsH] = plot_marginal_feature(data, yind, featureName, groupNames, graphicsH);
    
end
%%%%%%%% End for plotting marginals of feature 2 %%%%%%%%%%%%%



function graphicsH = plot_marginal_feature(data, xind, featureName, groupNames, graphicsH)

ica_fuse_defaults;
global AX_COLOR;
global UI_FONT_NAME;
global FIG_FG_COLOR;

fontSize = 0.05;

% Number of figures
graphCount = length(graphicsH);

% Plot marginal histogram of first feature
titleStr = ['Marginal histogram of ', featureName];
graphicsH(graphCount + 1).H = ica_fuse_getGraphics(titleStr, 'normal', titleStr, 'on');
axesH = axes('parent', graphicsH(graphCount + 1).H, 'units', 'normalized', 'position', [0.12 0.12 0.7 0.7]);

% Plot marginal of feature 1 for the first group
plot(ica_fuse_interp(xind, 3), ica_fuse_interp(data{1}, 3), 'c-', 'linewidth', 2);

if length(groupNames) > 1
    hold on;
    plot(ica_fuse_interp(xind, 3), ica_fuse_interp(data{2}, 3), 'b-', 'linewidth', 2);
    hold off;
end

set(axesH, 'fontunits', 'normalized', 'FontName', UI_FONT_NAME, 'fontsize', fontSize);
set(axesH, 'color', AX_COLOR);
set(axesH, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);

title(featureName, 'parent', axesH);
xlabel('Intensity Values', 'parent', axesH);
ylabel('Histogram Count', 'parent', axesH);
ica_fuse_legend(groupNames);

axis(axesH, 'tight');
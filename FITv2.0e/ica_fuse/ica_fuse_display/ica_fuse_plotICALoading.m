function ica_fuse_plotICALoading(varargin)
% Plot ICA Loading parameters

ica_fuse_defaults;
global FIG_FG_COLOR;
global FIG_FG_COLOR;
global UI_FONT_NAME; % font name
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size


titleColor = 'c';
axesTitle = '';
ica_loading_color = {'g'; 'r'; 'c'; 'm'; 'b'};
fontSize = 0.05;

for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'data')
        data = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'parent')
        axesHandle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'color')
        ica_loading_color = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'title')
        axesTitle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'titlecolor')
        titleColor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'xlabel')
        xLabel = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'groupnames')
        groupNames = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'featurenames')
        featureNames = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'yaxislocation')
        yAxisLocation = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'histdata')
        histData = varargin{ii + 1};
    end
    
end

if ~exist('data', 'var')
    error('Data variable must be present');
end

% Get the group data
group1 = data(1).dat;
group2 = data(2).dat;

if ~exist('axesHandle', 'var')
    tag = 'Timecourse';
    GraphicsHandle = ica_fuse_getGraphics(tag, 'timecourse', tag, 'on');
    axesHandle = get(GraphicsHandle, 'CurrentAxes');
end

handles = zeros(1, size(group1, 2)*2);
for nFea = 1:size(group1, 2)
    % Group1 text
    handles(2*nFea - 1) = plot((1:size(group1, 1)), group1(:, nFea), ['-', ica_loading_color{nFea}, 'o'], 'parent', axesHandle);
    hold on;
    % Group2 text
    handles(2*nFea) = plot((size(group1, 1) + 1:size(group1, 1) + size(group2, 1)), group2(:, nFea), ['-', ica_loading_color{nFea}, 'o'], 'parent', axesHandle);
    hold on;
end

% Plot dashed lines
axisYlim = get(axesHandle, 'ylim');
ySpacing = (axisYlim(2) - axisYlim(1))/4;
plotYY = (axisYlim(1):ySpacing:axisYlim(2));
plot(repmat(size(group1, 1), size(plotYY)), plotYY, '--y', 'parent', axesHandle);


% Draw text For Group1
yPos = 0.1;
xPos = 0.025*size(group1, 1) / (size(group1, 1) + size(group2, 1));
text(xPos, yPos, deblank(groupNames(1, :)), 'parent', axesHandle, 'units', 'normalized');

% Draw text For Group2
xPos = (size(group1, 1) + 0.025*size(group2, 1)) / (size(group1, 1) + size(group2, 1));
text(xPos, yPos, deblank(groupNames(2, :)), 'parent', axesHandle, 'units', 'normalized');

set(axesHandle, 'YColor', FIG_FG_COLOR);
set(axesHandle, 'XColor', FIG_FG_COLOR);
set(axesHandle, 'YAxisLocation', yAxisLocation);

if (size(group1, 2) > 1)
    ica_fuse_legend(handles(1:2:end), featureNames(1:size(group1, 2), :));
end
xlabel('Subjects', 'parent', axesHandle);
ylabel('ICA Loading (Arbitrary Units)', 'parent', axesHandle);

title(axesTitle, 'color',  titleColor, 'interpreter', 'tex', 'parent', axesHandle);

axis(axesHandle, 'tight');

% % If histData variable exists
if exist('histData', 'var')
    
    set(axesHandle, 'userdata', histData);
    
    % set context menu and put an option to plot the histograms for subjects
    % for the selected features
    cmenu = uicontextmenu;
    
    % Define the line and associate it with the context menu
    set(axesHandle, 'UIContextMenu', cmenu);
    % Define the context menu items
    histogramItem = uimenu(cmenu, 'Label', 'Plot histogram of features', 'Callback', ...
        {@histogramCallback, axesHandle});
    
end


%%%%%%%%%% Define Function callbacks %%%%%%%%%%%%%%%%%%%%%%%

function histogramCallback(hObject, event_data, handles)
% Histogram Plot Callback

ica_fuse_defaults;

global Z_THRESHOLD_HISTOGRAM;

histData = get(handles, 'userdata');

parentFigHandle = get(handles, 'parent');

% Get the related information from the histData structure
compNum = histData.compNum;
fusionFile = histData.fusionFile;
selGroupsVal = histData.selGroupsVal;
selFeaturesVal = histData.selFeaturesVal;
threshold = Z_THRESHOLD_HISTOGRAM;

outputDir = fileparts(fusionFile);
load(fusionFile);

% Form histParameters structure
histParameters.numGroups = fusionInfo.run_analysis.numGroups;
histParameters.numSubjects = fusionInfo.run_analysis.numSubjects;
histParameters.selGroupsVal = selGroupsVal;
histParameters.dataInfo = fusionInfo.run_analysis.dataInfo;
histParameters.groupNames = str2mat(histParameters.dataInfo.name);
histParameters.backReconstructFiles = fusionInfo.run_analysis.backReconstructFiles;
histParameters.mask_ind = fusionInfo.run_analysis.mask_ind;
histParameters.featureDataLength = fusionInfo.run_analysis.featureDataLength;
histParameters.normalize = fusionInfo.run_analysis.normalize;
histParameters.featureNormPara = fusionInfo.run_analysis.featureNormPara;
histParameters.outputDir = outputDir;

%% Check this for backward compatibility
if (isfield(fusionInfo.run_analysis, 'all_comb'))
    histParameters.all_comb = fusionInfo.run_analysis.all_comb;
else
    % Assuming all the features are stacked
    stackAllFeatures = 1;
    optimalFeatures = 0;
    if (length(fusionInfo.run_analysis.pcaFiles) > 1)
        optimalFeatures = 1;
    end
    % Get all the combinations
    histParameters.all_comb = ica_fuse_get_combinations(fusionInfo.run_analysis.numFeatures, stackAllFeatures, optimalFeatures);
end

if (~isfield(fusionInfo.run_analysis, 'dims'))
    [dims, voxels] = ica_fuse_getFeatureDIM(fusionInfo.run_analysis.mask_ind);
else
    dims = fusionInfo.run_analysis.dims;
    voxels = fusionInfo.run_analysis.voxels;
end

histParameters.dims = dims;
histParameters.voxels = voxels;

set(parentFigHandle, 'pointer', 'watch');

% Call histogram plot function
ica_fuse_display_histogram(histParameters, compNum, selFeaturesVal, threshold, 'feature');

set(parentFigHandle, 'pointer', 'arrow');

fprintf('\n');

function [figHandle] = plotIm(data, titleStr, cmap, axesTitle)
% Sub function to plot images in a sub plot

ica_fuse_defaults;
global AX_COLOR; % axes color
global LEGEND_COLOR; % LEGEND_COLOR for Matlab 7.0
global FIG_FG_COLOR;
global UI_FONT_NAME;
global UI_FONT_UNITS;
global UI_FONT_SIZE;

figHandle = ica_fuse_getGraphics(titleStr, 'normal', titleStr, 'on');

set(figHandle, 'colormap', cmap);

if length(size(data)) == 3
    numPlots = size(data, 3);
else
    numPlots = 1;
end

numRows = ceil(sqrt(numPlots));
numCols = ceil(numPlots/numRows);

fontSize = (0.05);

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
                axesHandle = axes('parent', figHandle, 'units', 'normalized', 'position', [0.1 0.1 0.7 0.7]);
            end
            
            set(axesHandle, 'units', 'normalized', 'color', AX_COLOR, 'FontName', UI_FONT_NAME, 'fontunits', ...
                'normalized', 'fontsize', fontSize);
            set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
            set(axesHandle, 'YTickLabel', []);
            set(axesHandle, 'YTick', []);
            set(axesHandle, 'XTickLabel', []);
            set(axesHandle, 'XTick', []);
            subject_str = ['Sub ', num2str(countN)];
            temp = squeeze(data(:, :, countN));
            temp = temp';
            ImageAxis = image(temp, 'parent', axesHandle, 'CDataMapping', 'scaled');
            
            set(axesHandle, 'clim', [min(data(:)), max(data(:))]);
            clear temp;
            axis(axesHandle, 'off');
            if exist('axesTitle', 'var')
                title(axesTitle, 'parent', axesHandle);
            else
                title(['Sub ', num2str(countN)], 'parent', axesHandle);
            end
        else
            break;
        end
        
    end
    
end

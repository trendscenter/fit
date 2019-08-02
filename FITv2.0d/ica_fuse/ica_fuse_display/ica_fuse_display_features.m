function ica_fuse_display_features(varargin)


plotType = 'timecourse';
time_course_color = 'm';
%cmap = colormap;
titleColor = 'c';
number_per_figure = 4;
plot_prev_next = 1;


for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'plot_data')
        plotData = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'title_color')
        titleColor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'color_map')
        cmap = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'time_course_color')
        time_course_color = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'ica_loading_color')
        ica_loading_color = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'number_per_figure')
        number_per_figure = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'slice_plane')
        slice_plane = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'plot_prev_next')
        plot_prev_next = varargin{ii + 1};
    end
end


ica_fuse_defaults;
global AX_COLOR; % axes color
global LEGEND_COLOR; % LEGEND_COLOR for Matlab 7.0
global FIG_FG_COLOR;
global UI_FONT_NAME;
global UI_FONT_UNITS;
global UI_FONT_SIZE;

%fontSize = 0.05;

% Get the axes position based on the number of images per figure
% Possible options include 4, 6 (Three per row)

% Number of plots
numPlots = length(plotData);

numCols = ceil(sqrt(number_per_figure));

numRows = round(number_per_figure / numCols);

% Number of rows and columns
%numRows = 2; numCols = ceil(number_per_figure / numRows);

% Number of figures possible
numFigures = ceil(numPlots / number_per_figure);

% Number of last axes
numLastAxes = mod(numPlots, number_per_figure);

matchIndex1 = strmatch('ica loading', lower(str2mat(plotData.plotType)), 'exact');
matchIndex2 = strmatch('image', lower(str2mat(plotData.plotType)), 'exact');

% Offsets
yOffset = 0.065; xOffset = 0.052;

% Axes width and height
axesWidth = ((1 - xOffset - numCols*xOffset) / numCols) - 0.01;

if ~isempty(matchIndex1) & ~isempty(matchIndex2)
    xOffset = 0.072;
    if numCols ~= 1
        axesWidth = ((1 - (numCols - 1)*1.5*xOffset - 2*xOffset) / numCols) - 0.01;
    else
        xOffset = 0.1;
        axesWidth = ((1 - 2*xOffset) / numCols) - 0.01;
    end
end

if isempty(matchIndex1) & ~isempty(matchIndex2)
    if numCols ~= 1
        axesWidth = ((1 - (numCols - 1)*1.75*xOffset - 2*xOffset) / numCols) - 0.01;
    else
        xOffset = 0.1;
        axesWidth = ((1 - 2*xOffset) / numCols) - 0.01;
    end
end

leftXOffset = xOffset;
axesHeight = ((1 - yOffset*numRows) / numRows) - 0.01;

fontSize = (0.025)*sqrt(numRows*numCols);

axesWidth = min([axesHeight, axesWidth]);
axesHeight = axesWidth;

if numLastAxes == 0
    numLastAxes = number_per_figure;
end

count = 0;
% Loop over  figures
for nFig = 1:numFigures
    
    % Get the graphics handle
    graphicsH(nFig).H = ica_fuse_getGraphics(['Fusion Analysis ', num2str(nFig)], 'graphics', ['JointICA ', num2str(nFig)], 'on');
    
    
    % get the number of axes
    if nFig ~= numFigures
        numAxes = number_per_figure;
    else
        numAxes = numLastAxes;
    end
    
    countData = 0;
    % Loop over number of axes
    for getRow = 1:numRows
        
        axesXOrigin = leftXOffset;
        
        for getCol = 1:numCols
            
            count = count + 1;
            countData = countData + 1;
            
            if count <= length(plotData)
                
                graphicsH(nFig).userdata(countData).plotType = plotData(count).plotType;
                
                axesYOrigin = (1 - getRow*axesHeight - getRow*yOffset);
                
                %axesXOrigin = leftXOffset + (getCol - 1)*axesWidth + (getCol - 1)*xOffset;
                
                axesH = axes('Parent', graphicsH(nFig).H, 'units', 'normalized', 'position', ...
                    [axesXOrigin, axesYOrigin, axesWidth, axesHeight], 'color', AX_COLOR, 'FontName', UI_FONT_NAME, ...
                    'fontunits', 'normalized', 'fontsize', fontSize);
                
                graphicsH(nFig).userdata(countData).axesPosition = get(axesH, 'position');
                
                % Figure data
                graphicsH(nFig).userdata(countData).data = plotData(count).data;
                graphicsH(nFig).userdata(countData).time_course_color = time_course_color;
                graphicsH(nFig).userdata(countData).titleColor = titleColor;
                graphicsH(nFig).userdata(countData).axesTitle = plotData(count).axesTitle;
                graphicsH(nFig).userdata(countData).yAxisLocation = 'left';
                graphicsH(nFig).userdata(countData).meanData = plotData(count).meanData;
                graphicsH(nFig).userdata(countData).meanDataLegend = plotData(count).meanDataLegend;
                graphicsH(nFig).userdata(countData).groupCompData = plotData(count).groupCompData;
                graphicsH(nFig).userdata(countData).groupCompLegend = plotData(count).groupCompLegend;
                graphicsH(nFig).userdata(countData).colorbarLIM = [];
                graphicsH(nFig).userdata(countData).axesCLIM = [];
                graphicsH(nFig).userdata(countData).cmap = cmap;
                graphicsH(nFig).userdata(countData).colorbarPos = [];
                graphicsH(nFig).userdata(countData).colorbarMinMaxText = [];
                graphicsH(nFig).userdata(countData).textLeftRight = [];
                graphicsH(nFig).userdata(countData).slice_plane = [];
                graphicsH(nFig).userdata(countData).ica_loading_color = [];
                graphicsH(nFig).userdata(countData).groupNames = [];
                graphicsH(nFig).userdata(countData).offset = [xOffset, yOffset];
                
                if strcmpi(plotData(count).plotType, 'image')
                    % Get the function to plot the image
                    
                    % Colorbar position
                    colorbarPos = get(axesH, 'position');
                    colorbarPos(1) = colorbarPos(1) + colorbarPos(3) + 0.4*xOffset;
                    colorbarPos(3) = 0.35*xOffset;
                    colorbarPos(2) = colorbarPos(2) + 0.2*colorbarPos(4);
                    colorbarPos(4) = 0.4*colorbarPos(4);
                    
                    
                    minClim = min(plotData(count).colorbarLim);
                    maxClim = 2*max(plotData(count).colorbarLim); %max(plotData(count).data(:));
                    graphicsH(nFig).userdata(countData).axesCLIM = [minClim, maxClim];
                    graphicsH(nFig).userdata(countData).colorbarLIM = plotData(count).colorbarLim;
                    graphicsH(nFig).userdata(countData).colorbarPos = colorbarPos;
                    
                    graphicsH(nFig).userdata(countData).colorbarMinMaxText = plotData(count).colorbarMinMaxText;
                    graphicsH(nFig).userdata(countData).textLeftRight = plotData(count).textLeftRight;
                    graphicsH(nFig).userdata(countData).slice_plane = slice_plane;
                    
                    if strcmpi(slice_plane, 'axial')
                        % Plot the image
                        ica_fuse_plotImage('parent', axesH, 'data', plotData(count).data, 'CDataMapping', 'scaled', 'axesCLIM', [minClim, maxClim], ...
                            'colormap', cmap, 'title', plotData(count).axesTitle, 'titlecolor', titleColor, 'colorbarLIM', plotData(count).colorbarLim, ...
                            'colorbarPosition', colorbarPos, 'colorbarMinMaxText', plotData(count).colorbarMinMaxText, ...
                            'text_left_right',  plotData(count).textLeftRight);
                    else
                        % Plot the image
                        ica_fuse_plotImage('parent', axesH, 'data', plotData(count).data, 'CDataMapping', 'scaled', 'axesCLIM', [minClim, maxClim], ...
                            'colormap', cmap, 'title', plotData(count).axesTitle, 'titlecolor', titleColor, 'colorbarLIM', plotData(count).colorbarLim, ...
                            'colorbarPosition', colorbarPos, 'colorbarMinMaxText', plotData(count).colorbarMinMaxText);
                    end
                    
                    axesXOrigin = colorbarPos(1) + colorbarPos(3) + xOffset;
                    
                elseif strcmpi(plotData(count).plotType, 'timecourse')
                    % get the function to plot the time course
                    
                    ica_fuse_plotTimecourse('parent', axesH, 'data', plotData(count).data, 'color', time_course_color, 'titleColor', ...
                        titleColor, 'title', plotData(count).axesTitle, 'YAxisLocation', 'left', ...
                        'meanData', plotData(count).meanData, 'meanDataLegend', plotData(count).meanDataLegend, ...
                        'groupCompData', plotData(count).groupCompData, 'groupCompLegend', ...
                        plotData(count).groupCompLegend);
                    
                    axesXOrigin = axesXOrigin + axesWidth + xOffset;
                    
                elseif strcmpi(plotData(count).plotType, 'snp')
                    
                    ica_fuse_plotSNP('parent', axesH, 'data', plotData(count).data, 'color', time_course_color, 'titleColor', ...
                        titleColor, 'title', plotData(count).axesTitle, 'xlabel', 'Subjects', 'ylabel', '');
                    
                    axesXOrigin = axesXOrigin + axesWidth + xOffset;
                    
                elseif strcmpi(plotData(count).plotType, 'ica loading')
                    
                    if isfield(plotData, 'selFeaturesVal') & isfield(plotData, 'selGroupsVal')
                        histData.selGroupsVal = plotData(count).selGroupsVal;
                        histData.selFeaturesVal = plotData(count).selFeaturesVal;
                        histData.compNum = plotData(count).compNum;
                        histData.fusionFile = plotData(count).fusionFile;
                        %                         if isfield(plotData(count), 'histData')
                        %                             histData.data = plotData(count).histData;
                        %                         end
                    end
                    
                    graphicsH(nFig).userdata(countData).ica_loading_color = ica_loading_color;
                    graphicsH(nFig).userdata(countData).groupNames = plotData(count).groupNames;
                    graphicsH(nFig).userdata(countData).selFeatureNames = plotData(count).selFeatureNames;
                    
                    % Plot ICA loading
                    ica_fuse_plotICALoading('parent', axesH, 'data', plotData(count).data, 'color', ica_loading_color, ...
                        'titleColor', titleColor, 'title', plotData(count).axesTitle, 'groupNames', ...
                        plotData(count).groupNames, 'YAxisLocation', 'left', 'histData', histData, 'featurenames', plotData(count).selFeatureNames);
                    
                    axesXOrigin = axesXOrigin + axesWidth + xOffset;
                    
                end
            end
            
        end
        
    end
    % end loop over number of axes
    
end
% end loop over number of figures


% Plot previous, next and exit push buttons
if (plot_prev_next)
    ica_fuse_plotNextPreviousExitButtons(graphicsH);
end

for nFig = 1:numFigures
    set(graphicsH(nFig).H, 'WindowButtonDownFcn', @mouseClickCallback);
end


function mouseClickCallback(hObject, event_data, handles)

ica_fuse_defaults;

% axes color
global AX_COLOR;

% FONT DEFAULTS
global UI_FONT_NAME; % font name
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size

% Identify selectionType
selectionType = get(hObject, 'SelectionType');

objTag = get(gco, 'Tag');

if strcmpi(objTag, 'legend')
    return;
end

% Object Type
objType = get(gco, 'Type');

% Identify double click
if (strcmpi(selectionType, 'normal') | strcmpi(selectionType, 'open')) & ...
        (strcmpi(objType, 'axes') | strcmpi(objType, 'image') | strcmpi(objType, 'line'))
    
    objData = get(hObject, 'userdata');
    
    figIndex = objData.index;
    
    userdata = objData.GraphicsHandle(figIndex).userdata;
    
    set(hObject, 'units', 'normalized');
    
    currentPoint = get(hObject, 'currentPoint');
    
    for nn = 1:length(userdata)
        
        if strcmpi(userdata(nn).plotType, 'image')
            colorbarPos = userdata(nn).colorbarPos;
            axesLocation = userdata(nn).axesPosition;
            % Include the axes until colorbar
            axesLocation(3) = colorbarPos(1) - axesLocation(1) + colorbarPos(3);
        else
            axesLocation = userdata(nn).axesPosition;
        end
        
        
        pointLiesInside = ica_fuse_check_point_inside(currentPoint, axesLocation);
        
        if pointLiesInside
            data = userdata(nn).data; % data
            titleColor = userdata(nn).titleColor; % title color
            time_course_color = userdata(nn).time_course_color; % timecourse color
            axesTitle = userdata(nn).axesTitle; % axes title
            yAxisLocation = userdata(nn).yAxisLocation; % y axis location
            meanData = userdata(nn).meanData; % mean data
            meanDataLegend = userdata(nn).meanDataLegend; % mean data legend
            groupCompData = userdata(nn).groupCompData; % group component data
            groupCompLegend = userdata(nn).groupCompLegend; % group component legend
            colorbarLim = userdata(nn).colorbarLIM;
            axesCLIM = userdata(nn).axesCLIM;
            cmap = userdata(nn).cmap;
            colorbarPos = userdata(nn).colorbarPos;
            colorbarMinMaxText = userdata(nn).colorbarMinMaxText;
            textLeftRight = userdata(nn).textLeftRight;
            slice_plane = userdata(nn).slice_plane;
            ica_loading_color = userdata(nn).ica_loading_color;
            groupNames = userdata(nn).groupNames;
            if (isfield(userdata, 'selFeatureNames'))
                featureNames = userdata(nn).selFeatureNames;
            end
            offset = userdata(nn).offset;
            
            
            [graphicsHandle] = ica_fuse_getGraphics(axesTitle, 'displaygui', axesTitle, 'on');
            
            axesH = axes('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'color', AX_COLOR, 'FontName', UI_FONT_NAME, ...
                'fontunits', UI_FONT_UNITS, 'fontsize', UI_FONT_SIZE);
            
            try
                set(hObject, 'pointer', 'watch');
                
                % Plot time course
                if strcmpi(userdata(nn).plotType, 'timecourse')
                    
                    ica_fuse_plotTimecourse('parent', axesH, 'data', data, 'color', time_course_color, 'titleColor', ...
                        titleColor, 'title', axesTitle, 'YAxisLocation', yAxisLocation, 'meanData', meanData, ...
                        'meanDataLegend', meanDataLegend, 'groupCompData', groupCompData, 'groupCompLegend', groupCompLegend);
                elseif strcmpi(userdata(nn).plotType, 'image')
                    % Colorbar position
                    colorbarPos = get(axesH, 'position');
                    colorbarPos(1) = colorbarPos(1) + colorbarPos(3) + 0.4*offset(1);
                    colorbarPos(3) = 0.35*offset(1);
                    % Plot the image
                    if strcmpi(slice_plane, 'axial')
                        ica_fuse_plotImage('parent', axesH, 'data', data, 'CDataMapping', 'scaled', 'axesCLIM', axesCLIM, ...
                            'colormap', cmap, 'title', axesTitle, 'titlecolor', titleColor, 'colorbarLIM', colorbarLim, ...
                            'colorbarPosition', colorbarPos, 'colorbarMinMaxText', colorbarMinMaxText, ...
                            'text_left_right',  textLeftRight);
                    else
                        ica_fuse_plotImage('parent', axesH, 'data', data, 'CDataMapping', 'scaled', 'axesCLIM', axesCLIM, ...
                            'colormap', cmap, 'title', axesTitle, 'titlecolor', titleColor, 'colorbarLIM', colorbarLim, ...
                            'colorbarPosition', colorbarPos, 'colorbarMinMaxText', colorbarMinMaxText);
                        
                    end
                elseif strcmpi(userdata(nn).plotType, 'ica loading')
                    ica_fuse_plotICALoading('parent', axesH, 'data', data, 'color', ica_loading_color, ...
                        'titleColor', titleColor, 'title', axesTitle, 'groupNames', ...
                        groupNames, 'YAxisLocation', yAxisLocation, 'featurenames', featureNames);
                    
                end
                set(hObject, 'pointer', 'arrow');
            catch
                set(hObject, 'pointer', 'arrow');
                ica_fuse_displayErrorMsg;
            end
            
            return;
        end
    end
end
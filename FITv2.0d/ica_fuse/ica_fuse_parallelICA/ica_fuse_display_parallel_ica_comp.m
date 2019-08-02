function ica_fuse_display_parallel_ica_comp(varargin)
% Display parallel ICA components

time_course_color = 'm';
titleColor = 'c';
number_per_figure = 4;

for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'plot_data')
        plotData = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'title_color')
        titleColor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'color_map')
        cmap = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'time_course_color')
        time_course_color = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'number_per_figure')
        number_per_figure = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'slice_plane')
        slice_plane = varargin{ii + 1};
    end
end


% Number of components
numComponents = length(plotData);
% Number of modalities
numModalities = length(plotData(1).feature);
% Total number of plots
numPlots = numComponents*(numModalities + 1);


% if (numModalities == 3)
%     numCols = numModalities;
% else
numCols = ceil(sqrt(number_per_figure));
%end

numRows = round(number_per_figure / numCols);
% Number of figures possible
numFigures = ceil(numComponents / number_per_figure);


% Offsets
yOffset = 0.01; xOffset = 0.065;

if number_per_figure == 1
    yOffset = 0.025;
    xOffset = 0.08;
end

% Font size
fontSize = (0.05)*sqrt(numRows*numCols);

if fontSize > 1
    fontSize = 0.85;
end

% Box widths
boxWidth = ((1 - xOffset - numCols*xOffset) / numCols);

boxHeight = ((1 - 2*yOffset*numRows) / numRows);


startXOffset = 0.6*xOffset;
count = 0;
% Loop over  figures
for nFig = 1:numFigures
    
    % Get the graphics handle
    figHandle = ica_fuse_getGraphics(['Parallel ICA ', num2str(nFig)], 'graphics', ['ParallelICA ', num2str(nFig)], 'on');
    
    % Loop over number of rows
    for getRow = 1:numRows
        
        boxXOrigin = startXOffset;
        % Loop over number of cols
        for getCol = 1:numCols
            
            count = count + 1;
            
            if count <= length(plotData)
                
                boxYOrigin = (1 - getRow*boxHeight - getRow*yOffset);
                
                % Box Position
                boxPos = [boxXOrigin, boxYOrigin, boxWidth, boxHeight];
                
                % Draw Loading coefficients
                plotObjects(figHandle, boxPos, plotData(count), cmap, titleColor, time_course_color, slice_plane, fontSize, number_per_figure);
            end
            
            boxXOrigin = boxXOrigin + boxWidth + xOffset;
            
        end
        % End loop over number of rows
    end
    % End loop over number of rows
    
    graphicsH(nFig).H = figHandle;
    graphicsH(nFig).userdata = get(figHandle, 'userdata');
    clear userData;
    
end
% End loop over figures


% Plot previous, next and exit push buttons
ica_fuse_plotNextPreviousExitButtons(graphicsH);

for nFig = 1:numFigures
    set(graphicsH(nFig).H, 'WindowButtonDownFcn', {@mouseClickCallback, graphicsH(nFig).H});
end

function plotObjects(figHandle, boxPos, dataStruct, cmap, titleColor, time_course_color, slice_plane, fontSize, number_per_figure)
% Use boxposition and draw objects

ica_fuse_defaults;

global AX_COLOR;
global UI_FONT_NAME;
global FIG_FG_COLOR;
global FIG_FG_COLOR;

userdata = get(figHandle, 'userdata');

xOffset = 0.04;
yOffset = 0.06;

if number_per_figure == 1
    xOffset = 0.06;
elseif number_per_figure == 4
    xOffset = 0.045;
end

startYOffset = 0.03;

% numFeatures
numFeatures = length(dataStruct.feature);

numAxes = numFeatures + 1;

axes1Height = 0.3*(boxPos(4) - 2*yOffset - startYOffset);
axes2Height = 0.7*(boxPos(4) - 2*yOffset - startYOffset);

fontSize1 = axes2Height*fontSize;
fontSize2 = axes1Height*fontSize;

axesWidth = boxPos(3) - xOffset;
% Loading coefficient axes

axesOb(1).pos = [boxPos(1) + xOffset, boxPos(2) + boxPos(4) - startYOffset - axes1Height, axesWidth, axes1Height];
if numAxes == 2
    axesWidth = axesOb(1).pos(3);
    axesOb(2).pos = [axesOb(1).pos(1) + 0.5*axesOb(1).pos(3) - 0.5*axesWidth, boxPos(2) + yOffset, axesWidth, 0.9*axes2Height];
    %legendString{length(legendString) + 1} = dataStruct.feature(1).feature_name;
    colors = {'.-g'};
    %legendString{length(legendString) + 1} = dataStruct.feature(1).feature_name;
elseif numAxes == 3
    %axesOb(1).pos = [boxPos(1) + xOffset, boxPos(2) + boxPos(4) - startYOffset - axes1Height, axesWidth, axes1Height];
    axesWidth = 0.425*axesOb(1).pos(3);
    axesWidth = min([axesWidth, axes2Height]); axes2Height = axesWidth;
    axesOb(2).pos = [axesOb(1).pos(1), boxPos(2) + yOffset, axesWidth, axes2Height];
    if (strcmpi(dataStruct.feature(1).modality, 'fmri') || strcmpi(dataStruct.feature(1).modality, 'smri')) && ...
            strcmpi(dataStruct.feature(2).modality, 'gene')
        axesOb(3).pos = [axesOb(2).pos(1) + axesOb(2).pos(3) + 1.85*xOffset, boxPos(2) + yOffset, axesWidth, axes2Height];
    else
        axesOb(3).pos = [axesOb(2).pos(1) + axesOb(2).pos(3) + 1.4*xOffset, boxPos(2) + yOffset, axesWidth, axes2Height];
    end
    
    colors = {'.-g', '.-r'};
    
else
    
    %         for nA = 1:4
    %             if (nA == 1)
    %                 sh = subplot(2, 3, 1:3);
    %             else
    %                 sh = subplot(2, 3, 3 + nA - 1);
    %             end
    %             axesOb(nA).pos = get(sh, 'position');
    %         end
    
    
    axesWidth = 0.29*axesOb(1).pos(3);
    axesWidth = min([axesWidth, axes2Height]); axes2Height = axesWidth;
    axesOb(2).pos = [axesOb(1).pos(1), boxPos(2) + yOffset, axesWidth, axes2Height];
    axesOb(3).pos = [axesOb(2).pos(1) + axesOb(2).pos(3) + 1.4*xOffset, boxPos(2) + yOffset, axesWidth, axes2Height];
    axesOb(4).pos = [axesOb(3).pos(1) + axesOb(3).pos(3) + 1.4*xOffset, boxPos(2) + yOffset, axesWidth, axes2Height];
    
    colors = {'.-g', '.-r', '.-y'};
    
end

legendString = cellstr(char(dataStruct.feature.feature_name));

countData = length(userdata);
for nA = 1:numAxes
    countData = countData + 1;
    axesH = axes('Parent', figHandle, 'units', 'normalized', 'position', axesOb(nA).pos, 'color', AX_COLOR, ...
        'Fontname', UI_FONT_NAME, 'fontunits', 'normalized', 'fontsize', fontSize);
    
    % Store user data
    userdata(countData).axesPosition = axesOb(nA).pos;
    userdata(countData).time_course_color = time_course_color;
    userdata(countData).titleColor = titleColor;
    userdata(countData).colorbarLIM = [];
    userdata(countData).axesCLIM = [];
    userdata(countData).cmap = cmap;
    userdata(countData).colorbarPos = [];
    userdata(countData).colorbarMinMaxText = [];
    userdata(countData).textLeftRight = [];
    userdata(countData).slice_plane = [];
    userdata(countData).offset = [xOffset, yOffset];
    userdata(countData).numFeatures = numFeatures;
    userdata(countData).ica_loading_colors = colors;
    userdata(countData).ica_loading_legend = legendString;
    
    try
        userdata(countData).meanData = dataStruct.feature(nA - 1).meanData;
    catch
        userdata(countData).meanData  = [];
    end
    
    if nA == 1
        
        userdata(countData).data = dataStruct;
        userdata(countData).axesTitle = '';
        userdata(countData).plotType = 'ica loading';
        % Plot parallel ICA loading coefficients
        ica_fuse_plotParaICALoading('data', dataStruct, 'parent', axesH, 'num_features', numFeatures, 'colors', colors, ...
            'legend_string', legendString, 'titleColor', titleColor);
    else
        
        if strcmpi(dataStruct.feature(nA - 1).modality, 'eeg')
            userdata(countData).xlabel = 'Time';
            userdata(countData).ylabel = 'Data Units';
        else
            userdata(countData).xlabel = 'SNPs';
            userdata(countData).ylabel = 'Z-score';
        end
        % Figure data
        userdata(countData).data = dataStruct.feature(nA - 1).data;
        userdata(countData).axesTitle = dataStruct.feature(nA - 1).titleStr;
        
        if strcmpi(dataStruct.feature(nA - 1).modality, 'fmri') || strcmpi(dataStruct.feature(nA - 1).modality, 'smri')
            userdata(countData).plotType = 'image';
            % Plot Image
            axesPos = get(axesH, 'position');
            axesWidth = min([axesPos(3), axesPos(4)]);
            axesPos(3:4) = axesWidth;
            set(axesH, 'position', axesPos);
            
            % Colorbar position
            colorbarPos = get(axesH, 'position');
            colorbarPos(1) = colorbarPos(1) + colorbarPos(3) + 0.4*xOffset;
            colorbarPos(3) = 0.35*xOffset;
            
            colorbarLim = dataStruct.feature(nA - 1).colorbarLim;
            
            % Colorbar min max text
            minClim = min(colorbarLim);
            maxClim = 2*max(colorbarLim);
            colorbarMinMaxText = str2mat(num2str(dataStruct.feature(nA - 1).minICAIm), num2str(dataStruct.feature(nA - 1).maxICAIm));
            
            axesCLIM = [minClim, maxClim];
            
            % Update user data
            userdata(countData).colorbarLIM = colorbarLim;
            userdata(countData).axesCLIM = axesCLIM;
            userdata(countData).cmap = cmap;
            userdata(countData).colorbarPos = colorbarPos;
            userdata(countData).colorbarMinMaxText = colorbarMinMaxText;
            userdata(countData).textLeftRight = dataStruct.feature(nA - 1).text_left_right;
            userdata(countData).slice_plane = slice_plane;
            
            if strcmpi(slice_plane, 'axial')
                % Plot the image
                ica_fuse_plotImage('parent', axesH, 'data', dataStruct.feature(nA - 1).data, 'CDataMapping', 'scaled', 'axesCLIM', ...
                    axesCLIM, 'colormap', cmap, 'title', dataStruct.feature(nA - 1).titleStr, 'titlecolor', titleColor, ...
                    'colorbarLIM', colorbarLim, 'colorbarPosition', colorbarPos, 'colorbarMinMaxText', ...
                    colorbarMinMaxText, 'text_left_right', dataStruct.feature(nA - 1).text_left_right);
            else
                % Plot the image
                ica_fuse_plotImage('parent', axesH, 'data', dataStruct.feature(nA - 1).data, 'CDataMapping', 'scaled', 'axesCLIM', ...
                    axesCLIM, 'colormap', cmap, 'title', dataStruct.feature(nA - 1).titleStr, 'titlecolor', titleColor, ...
                    'colorbarLIM', colorbarLim, 'colorbarPosition', colorbarPos, 'colorbarMinMaxText', ...
                    colorbarMinMaxText);
            end
            
        elseif strcmpi(dataStruct.feature(nA - 1).modality, 'eeg')
            
            userdata(countData).plotType = 'eeg';
            ica_fuse_plotTimecourse('parent', axesH, 'data', dataStruct.feature(nA - 1).data, 'color', time_course_color, 'titleColor', ...
                titleColor, 'title', dataStruct.feature(nA - 1).titleStr, 'YAxisLocation', 'left', ...
                'meanData', userdata(countData).meanData, 'groupcompdata', []);
            ica_fuse_legend('Mean', 'Comp');
            xlabel(userdata(countData).xlabel);
            yh = ylabel(userdata(countData).ylabel, 'units', 'normalized');
            yhpos = get(yh, 'position');
            set(yh, 'position', [1-yhpos(1), yhpos(2), yhpos(3)]);
            
        else
            
            userdata(countData).plotType = 'snp';
            % Plot SNPs
            ica_fuse_plotSNP('data', dataStruct.feature(nA - 1).data, 'parent', axesH, 'xlabel', userdata(countData).xlabel, 'ylabel', ...
                userdata(countData).ylabel, 'titleColor', titleColor, 'color', time_course_color, 'title', dataStruct.feature(nA - 1).titleStr);
            
        end
        
    end
    
end


% Set userdata to figure
set(figHandle, 'userdata', userdata);

%%%%%%%%%%%%% Function Callbacks %%%%%%%%%%%%%%%

function mouseClickCallback(hObject, event_data, handles, meanData)

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
            colorbarLim = userdata(nn).colorbarLIM;
            axesCLIM = userdata(nn).axesCLIM;
            cmap = userdata(nn).cmap;
            colorbarPos = userdata(nn).colorbarPos;
            colorbarMinMaxText = userdata(nn).colorbarMinMaxText;
            textLeftRight = userdata(nn).textLeftRight;
            slice_plane = userdata(nn).slice_plane;
            offset = userdata(nn).offset;
            numFeatures = userdata(nn).numFeatures;
            colors = userdata(nn).ica_loading_colors;
            legendString = userdata(nn).ica_loading_legend;
            % SNPS x and y labels
            x_label = userdata(nn).xlabel;
            y_label = userdata(nn).ylabel;
            
            % Open figure
            [graphicsHandle] = ica_fuse_getGraphics(axesTitle, 'displaygui', axesTitle, 'on');
            
            axesH = axes('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'color', AX_COLOR, 'FontName', UI_FONT_NAME, ...
                'fontunits', UI_FONT_UNITS, 'fontsize', UI_FONT_SIZE);
            
            try
                set(hObject, 'pointer', 'watch');
                
                % Plot time course
                if strcmpi(userdata(nn).plotType, 'snp')
                    % Plot SNPs
                    ica_fuse_plotSNP('data', data, 'parent', axesH, 'xlabel', x_label, 'ylabel', ...
                        y_label, 'titleColor', titleColor, 'color', time_course_color, 'title', axesTitle);
                elseif strcmpi(userdata(nn).plotType, 'eeg')
                    ica_fuse_plotTimecourse('parent', axesH, 'data', data, 'color', time_course_color, 'titleColor', ...
                        titleColor, 'title', axesTitle, 'YAxisLocation', 'left', 'meanData', userdata(nn).meanData, 'groupcompdata', []);
                    ica_fuse_legend('Mean', 'Comp');
                    xlabel(x_label);
                    yh = ylabel(y_label, 'units', 'normalized');
                    yhpos = get(yh, 'position');
                    set(yh, 'position', [1-yhpos(1), yhpos(2), yhpos(3)]);
                elseif strcmpi(userdata(nn).plotType, 'image')
                    % Colorbar position
                    colorbarPos = get(axesH, 'position');
                    colorbarPos(1) = colorbarPos(1) + colorbarPos(3) + 0.5*offset(1);
                    colorbarPos(3) = 0.5*offset(1);
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
                    % Plot parallel ICA loading coefficients
                    ica_fuse_plotParaICALoading('data', data, 'parent', axesH, 'num_features', numFeatures, 'colors', colors, ...
                        'legend_string', legendString, 'titleColor', titleColor);
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
function ica_fuse_plotTimecourse(varargin)
% Plot the timecourse

ica_fuse_defaults;
global FIG_FG_COLOR;
global UI_FONT_NAME; % font name
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size


titleColor = 'c';
axesTitle = '';
time_course_color = 'm';
yAxisLocation = 'left';
fontSize = 0.05;
legendString = {};
line_width = 1.5;

for ii = 1:2:nargin
    
    if strcmpi(varargin{ii}, 'data')
        data = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'parent')
        axesHandle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'color')
        time_course_color = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'title')
        axesTitle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'titlecolor')
        titleColor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'yaxislocation')
        yAxisLocation = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'meandata')
        meanData = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'meandatalegend')
        legendString = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'groupcompdata')
        groupCompData = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'groupcomplegend')
        groupCompLegend = varargin{ii + 1};
    end
    
end

if ~exist('data', 'var')
    error('Data variable must be present');
end

if size(data, 1) == 1 & size(data, 2) > 1
    data = data';
end

if size(data, 2) == 1
    xAxis = (1:length(data))';
    yAxis = data;
else
    xAxis = data(:, 1);
    yAxis = data(:, 2);
end

if exist('meanData', 'var')
    meanDataColors = {'.-y', '.-m', '.-r', '.-g', '.-b'};
    for nn = 1:size(meanData, 2)
        colorCount = mod(nn, length(meanDataColors));
        if colorCount == 0
            colorCount = length(meanDataColors);
        end
        plot(xAxis, meanData(:, nn), meanDataColors{colorCount}, 'linewidth', line_width);
        hold on;
    end
    hold on;
end

if exist('axesHandle', 'var')
    if isempty(groupCompData)
        plot(xAxis, yAxis, time_course_color, 'parent', axesHandle, 'linewidth', line_width);
        legendString{length(legendString) + 1} = 'Component Signal';
    else
        
        % Component data colors
        compDataColors = cell(1, length(meanDataColors));
        for nn = 1:length(meanDataColors)
            compDataColors{nn} = meanDataColors{end - nn + 1};
        end
        
        for nn = 1:size(meanData, 2)
            colorCount = mod(nn, length(compDataColors));
            if colorCount == 0
                colorCount = length(compDataColors);
            end
            plot(xAxis, squeeze(groupCompData(:, 2, nn)), compDataColors{colorCount}, 'parent', axesHandle, 'linewidth', line_width);
            legendString{length(legendString) + 1} = groupCompLegend{nn};
            hold on;
        end
    end
else
    tag = 'Timecourse';
    [GraphicsHandle] = ica_fuse_getGraphics(tag, 'timecourse', tag, 'on');
    plot_handle = plot(xAxis, yAxis, time_course_color, 'linewidth', line_width);
    axesHandle = get(plot_handle, 'parent');
end

hold off;

% Prepare legend string
if isempty(legendString)
    legendString{1} = 'Component Signal';
    %else
    %legendString{end + 1} = ['Component Signal'];
end

ica_fuse_legend(legendString);

grid(axesHandle, 'on');
axis(axesHandle, 'tight');

title(axesTitle, 'color',  titleColor, 'parent', axesHandle);

set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR, 'YAxisLocation', yAxisLocation);
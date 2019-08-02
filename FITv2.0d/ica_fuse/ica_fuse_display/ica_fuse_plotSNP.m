function ica_fuse_plotSNP(varargin)
% Plot SNPS

% Load defaults
ica_fuse_defaults;
global FIG_FG_COLOR;
global FIG_FG_COLOR;


titleColor = 'c';
axesTitle = '';
time_course_color = 'm';
yAxisLocation = 'left';
xlabelStr = 'SNPs';
ylabelStr = 'Z-score';
line_width = 1.5;

% Loop over number of arguments
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
    elseif strcmpi(varargin{ii}, 'xlabel')
        xlabelStr = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'ylabel')
        ylabelStr = varargin{ii + 1};
    end
end
% End loop over number of arguments


% Plot SNPs
plot(detrend(data, 0), time_course_color, 'parent', axesHandle, 'linewidth', line_width);
xlabel(xlabelStr, 'parent', axesHandle);
ylabel(ylabelStr, 'parent', axesHandle);
title(axesTitle, 'parent', axesHandle, 'color', titleColor);
axis(axesHandle, 'tight');
set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
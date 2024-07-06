function ica_fuse_legend(varargin)

ica_fuse_defaults;
global LEGEND_COLOR;

try
    legendH = legend(varargin{:});

    set(legendH, 'textColor', LEGEND_COLOR);
    set(legendH, 'EdgeColor', LEGEND_COLOR);
    set(legendH, 'Location', 'best');
catch

    legend(varargin{:});
end
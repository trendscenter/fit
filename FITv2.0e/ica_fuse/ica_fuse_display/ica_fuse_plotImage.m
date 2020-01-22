function ica_fuse_plotImage(varargin)

ica_fuse_defaults;
global FIG_FG_COLOR;
global UI_FONT_NAME; % font name
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size

titleText = '';
cDataMapping = 'scaled';
titleColor = 'c';

fontSize = 0.05;

% get input variables
for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'parent')
        axesHandle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'data')
        data = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'cdatamapping')
        cDataMapping = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'axesclim')
        axesClim = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'colormap')
        cmap = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'titlecolor')
        titleColor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'title')
        titleText = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'colorbarlim')
        colorbarLim = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'colorbarposition')
        colorbarPosition = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'colorbarminmaxtext')
        colorbarMinMaxText = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'text_left_right')
        text_left_right = varargin{ii + 1};
    end
    
end
% end for getting input variables

% Check the necessary vars
if ~exist('axesHandle', 'var')
    error('axes handle must be passed');
end


if ~exist('data', 'var')
    error('data must be present to plot the image');
end
% End for checking necessary vars

%%%%%%%%%%%%%%%%%%%%
% Plot image
%%%%%%%%%%%%%%%%%%%
ImageAxis = image(data, 'parent', axesHandle, 'CDataMapping', cDataMapping);

if ~exist('axesClim', 'var')
    minClim = min(colorbarLim);
    maxClim = 2*max(colorbarLim);
    %     minClim = min(data(:));
    %     maxClim = max(data(:));
    axesClim = [minClim, maxClim];
end

if exist('cmap', 'var')
    handleFig = get(axesHandle, 'parent');
    set(handleFig, 'Colormap', cmap);
end

set(axesHandle, 'CLIM', axesClim);

axis(axesHandle , 'off');
axis(axesHandle, 'image');
title(titleText, 'color',  titleColor, 'parent', axesHandle);


if exist('text_left_right', 'var')
    if ~isempty(text_left_right)
        [axesPos] =   get(axesHandle, 'position');
        diff_xPos = colorbarPosition(1) + colorbarPosition(3) - axesPos(1) - axesPos(3);
        xPos = 1 + 0.1*diff_xPos;
        yPos = 0.5;
        text(xPos, yPos, text_left_right(1), 'units', 'normalized', 'parent', axesHandle, 'color', FIG_FG_COLOR,  ...
            'fontsize', fontSize, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
            'FontName', UI_FONT_NAME, 'fontunits', 'normalized',  'fontweight', 'bold');
        
        xPos = -0.05;
        text(xPos, yPos, text_left_right(2), 'units', 'normalized', 'parent', axesHandle, 'color', FIG_FG_COLOR,  ...
            'fontsize', fontSize, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
            'FontName', UI_FONT_NAME, 'fontunits', 'normalized', 'fontweight', 'bold');
    end
    
end

if exist('colorbarPosition', 'var')
    %%%%%%%%%%%%%%%%%%%%%
    % Plot colorbar
    %%%%%%%%%%%%%%%%%%%%%%
    
    fields_props = {'color', 'Fontname', 'fontunits', 'fontsize'};
    
    vals = cell(1, 2*length(fields_props));
    for n = 1:length(fields_props)
        vals{2*n - 1} = fields_props{n};
        vals{2*n} = get(axesHandle, fields_props{n});
    end
    
    % cbAxesHandle = axes('parent', get(axesHandle, 'parent'), vals{:});
    % set(cbAxesHandle, 'position', colorbarPosition);   
    
    if ((ica_fuse_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
        cbAxesHandle = axes('parent', get(axesHandle, 'parent'), vals{:}, 'units', 'normalized');
        set(cbAxesHandle, 'position', colorbarPosition);
        axis(cbAxesHandle, 'off');
        colorbarHandle = colorbar('peer', cbAxesHandle);
    else
        colorbarHandle = colorbar;
    end
    
    set(colorbarHandle, 'position', colorbarPosition);
    set(colorbarHandle, 'units', 'pixels');    
    
    ChildH = get(colorbarHandle, 'Children');
    imInd = strmatch('image', lower(get(ChildH, 'Type')), 'exact');
    set(ChildH(imInd), 'YData', axesClim);
    
    if exist('colorbarLim', 'var')        
%         ChildH = get(colorbarHandle,'Children');
%         imInd = strmatch('image', lower(get(ChildH, 'Type')), 'exact');
%         set(ChildH(imInd), 'YData', colorbarLim);
        set(colorbarHandle, 'YLim', colorbarLim);
    end
    
    set(colorbarHandle, 'YTick', []);
    set(colorbarHandle, 'YTickLabel', []);
     set(colorbarHandle, 'XTickLabel', []);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Plot Text %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist('colorbarMinMaxText', 'var')
        colorbarMin = deblank(colorbarMinMaxText(1, :));
        colorbarMax = deblank(colorbarMinMaxText(2, :));
        % Plot colorbar text
        if ((ica_fuse_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
            % Maximum
            xPos = 0.01; yPos = 1.1;
            text(xPos, yPos, colorbarMax, 'units', 'normalized', 'parent', colorbarHandle, 'color', FIG_FG_COLOR,  ...
                'fontsize', fontSize, 'HorizontalAlignment', 'left', ...
                'FontName', UI_FONT_NAME, 'fontunits', 'normalized');
            
            % Minimum
            yPos = -0.1;
            text(xPos, yPos, colorbarMin, 'units', 'normalized', 'parent', colorbarHandle, 'color', FIG_FG_COLOR,  ...
                'fontsize', fontSize, 'HorizontalAlignment', 'left', ...
                'FontName', UI_FONT_NAME, 'fontunits', 'normalized');
        else
            labelsH = get(colorbarHandle, 'label');
            set(colorbarHandle, 'color', FIG_FG_COLOR);
            set(labelsH, 'units', 'normalized');
            title(colorbarMax, 'parent', colorbarHandle, 'color', FIG_FG_COLOR, 'Horizontalalignment', 'center');
            set(labelsH, 'string', colorbarMin, 'rotation',0);
            set(labelsH, 'position', [0.45,-0.02, 0]);
            set(colorbarHandle, 'YTickLabel', []);            
        end
    end
    
end

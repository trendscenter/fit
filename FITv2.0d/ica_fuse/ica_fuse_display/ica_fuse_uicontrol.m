function [handles] = ica_fuse_uicontrol(varargin)
% Applies default background and foreground colors for uicontrol using
% ica_fuse_defaults.m file
%
% Input:
% Same inputs as for uicontrol.
%
% Output:
% 1. handles vector
% On Matlab 7 and higher for popup box panels handle followed by popup
% handle will be returned.
% 


% load defaults
ica_fuse_defaults;

global UI_BG_COLOR;
global UI_FG_COLOR;
global BUTTON_BG_COLOR;
global BUTTON_FG_COLOR;
global FIG_BG_COLOR;
global FIG_FG_COLOR;
global UI_FONT_NAME;
global UI_FONT_UNITS;
global UI_FONT_SIZE;

% slider step
tag = '';
horizontal_alignment = 'center';
labels = '';
callback_function = '';
handle_visibility = 'on';
handle_enable = 'on';
fontSize = UI_FONT_SIZE;
fontName = UI_FONT_NAME;
fontUnits = UI_FONT_UNITS;
fontWeight = 'normal';
userData = [];
sliderStep = [0.01, 0.1];
maxValue = 1;
minValue = 0;
toolTipString = '';
value = 0;
panel_tag = 'panel';

% loop over number of arguments
for ii = 1:2:nargin    
    if strcmpi(varargin{ii}, 'tag')
        % tag for controls
        tag = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'backgroundcolor')
        % background color
        backgroundcolor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'foregroundcolor')
        % foreground color
        foregroundcolor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'horizontalalignment')
        % horizontal alignment
        horizontal_alignment = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'string')
        % string
        labels = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'parent')
        % parent
        parentHandle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'callback')
        % callback function
        callback_function = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'visible')
        % handle visibility
        handle_visibility = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'enable')
        % handle enable
        handle_enable = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'fontsize')
        % font size
        fontSize = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'fontunits')
        % font units
        fontUnits = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'fontname')
        % font name
        fontName = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'fontweight')
        % font weight
        fontWeight = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'style')
        % style
        styleType = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'units')
        % units
        handleUnits = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'position')
        % position
        position = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'userdata')
        % userdata
        userData = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'max')
        % maximum value
        maxValue = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'min')
        % minimum value
        minValue = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'sliderstep')
        % slider step
        sliderStep = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'tooltipstring')
        % tool tip string
        toolTipString = varargin{ii + 1};        
    elseif strcmpi(varargin{ii}, 'value')
        % value
        value = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'panel_tag')
        % panel tag
        panel_tag = varargin{ii + 1};
    end           
end

% set the default uicontrol colors
if strcmpi(styleType, 'pushbutton')
    if ~exist('backgroundcolor', 'var')
        backgroundcolor = BUTTON_BG_COLOR;
    end    
    if ~exist('foregroundcolor', 'var')
        foregroundcolor = BUTTON_FG_COLOR;
    end
    
else
    if ~exist('backgroundcolor', 'var')
        backgroundcolor = UI_BG_COLOR;
    end
    
    if ~exist('foregroundcolor', 'var')
        foregroundcolor = UI_FG_COLOR;
    end
    
end
% end for naming control colors

% set the default value
if strcmpi(styleType, 'popup') | strcmpi(styleType, 'popupmenu')
    if value == 0
        value = 1;
    end
end

if strcmpi(styleType, 'listbox')
    if ~isempty(value)
        if value == 0
            value = 1;
        end
    end
end
% end for setting the default value

% get the version of Matlab
whichVersion = str2num(version('-release'));

if isempty(whichVersion)
    whichVersion = 14;
end
     
% check the version of Matlab
if whichVersion <= 13        
    controlHandle = uicontrol('parent', parentHandle, 'BackgroundColor', backgroundcolor, 'ForegroundColor', foregroundcolor, ...
        'string', labels, 'units', handleUnits, 'position', position, 'fontunits', fontUnits, 'fontname', fontName, 'FontSize', ...
        fontSize, 'style', styleType, 'callback', callback_function, 'Visible', handle_visibility, 'enable', handle_enable, ...
        'tag', tag, 'min', minValue, 'max', maxValue, 'sliderstep', sliderStep, 'userdata', userData, 'value', value, ...
        'HorizontalAlignment', horizontal_alignment, 'toolTipString', toolTipString);
    handles(1) = controlHandle;
else
    % For versions greater than 13 plot popup on a panel
    
    if strcmpi(styleType, 'popup') | strcmpi(styleType, 'popupmenu')
        % plot popup box on panel for Matlab version 14 and later to fix
        % flashing of popup box on Mac and Linux
        panelHandle = uipanel('Parent', parentHandle, 'units', 'normalized', 'position', position, 'BackgroundColor', backgroundcolor, ...
            'bordertype', 'none', 'tag', panel_tag);        
        controlHandle = uicontrol('parent', panelHandle, 'BackgroundColor', backgroundcolor, 'ForegroundColor', foregroundcolor, ...
            'string', labels, 'units', handleUnits, 'position', [0 0 1 1], 'fontunits', fontUnits, 'fontname', fontName, 'FontSize', ...
            fontSize, 'style', styleType, 'callback', callback_function, 'Visible', handle_visibility, 'enable', handle_enable, ...
            'tag', tag, 'min', minValue, 'max', maxValue, 'sliderstep', sliderStep, 'userdata', userData, 'value', value, ...
            'HorizontalAlignment', horizontal_alignment, 'toolTipString', toolTipString);                
        handles(1) = panelHandle;
        handles(2) = controlHandle;
    else
        % others use the same
        controlHandle = uicontrol('parent', parentHandle, 'style', styleType, 'units', handleUnits, 'BackgroundColor', backgroundcolor, ...
            'ForegroundColor', foregroundcolor, 'string', labels, 'position', position, 'fontunits', fontUnits, ...
            'fontname', fontName, 'FontSize', fontSize, 'callback', callback_function, ...
            'Visible', handle_visibility, 'enable', handle_enable, 'tag', tag, 'min', minValue, 'max', maxValue, ...
            'sliderstep', sliderStep, 'userdata', userData, 'value', value, 'HorizontalAlignment', ...
            horizontal_alignment, 'toolTipString', toolTipString);
        handles(1) = controlHandle;        
    end
end
% end for plotting controls

set(controlHandle, 'HorizontalAlignment', horizontal_alignment);
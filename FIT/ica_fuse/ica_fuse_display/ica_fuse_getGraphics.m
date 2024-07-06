function [GraphicsHandle] = ica_fuse_getGraphics(TitleFig, GraphicsType, Tag, visible)
% Get the graphics window

% load defaults
ica_fuse_defaults;


global FIG_BG_COLOR; % Figure background color
global FIG_FG_COLOR; % Figure foreground color
global UI_BG_COLOR; % background color for user interface controls except pushbutton
global BUTTON_BG_COLOR; % background color for pushbutton
global UI_FG_COLOR; % Foreground color for controls except pushbutton
global BUTTON_FG_COLOR; % Foreground color for pushbutton

% axes color
global AX_COLOR;
global LEGEND_COLOR;

% FONT DEFAULTS
global UI_FONT_NAME; % font name
global UI_FONT_UNITS; % font units
global UI_FONT_SIZE; % font size

global SCREENSIZE;

if ~exist('TitleFig','var')
    TitleFig = 'Graphics';
end

if ~exist('GraphicsType','var')
    GraphicsType = 'normal';
end

if ~exist('Tag','var')
    Tag = 'Graphics';
end

if ~exist('visible', 'var')
    visible = 'on';
end

screenSize = SCREENSIZE;
% get the minimum screen size pixels
S = min([screenSize(3), screenSize(4)]);
%[R] = ica_fuse_aspectRatio; 
R = [1 1];
paperOrient = 'portrait';

switch lower(GraphicsType)
   
    case 'normal'
        extendLeft = round(S*R(1)*.55);
        extendUp = round(S*R(2)*.55);
        x0= (screenSize(3)/2)-(extendLeft/2);
        y0 = (screenSize(4)/2)-(extendUp/2);
        Rect = [x0 y0 extendLeft extendUp];
    case 'graphics'
        extendLeft = round(S*R(1)*.85);
        extendUp = round(S*R(2)*.85);
        x0= (screenSize(3)/2)-(extendLeft/2);
        y0 = (screenSize(4)/2)-(extendUp/2);
        Rect = [x0 y0 extendLeft extendUp];
    case 'displaygui'
        extendLeft = round(S*R(1)*.6);
        extendUp = round(S*R(2)*.6);
        x0= (screenSize(3)/2)-(extendLeft/2);
        y0 = (screenSize(4)/2)-(extendUp/2);
        Rect = [x0 y0 extendLeft extendUp];
    case 'statusbar'
        extendLeft = round(S*R(1)*.2);
        extendUp = round(S*R(2)*.3);
        x0= screenSize(3)*.1;
        y0 = screenSize(4)*.1;
        Rect = [x0 y0 extendLeft extendUp];
    case 'timecourse'
        extendLeft = round(screenSize(3)*R(1)*.9);
        extendUp = round(S*R(2)*.3);
        x0= screenSize(3)*.05;
        y0 = screenSize(4)*.5;
        Rect = [x0 y0 extendLeft extendUp];  
        paperOrient = 'landscape';
    otherwise        
        extendLeft = round(S*R(1)*.5);
        extendUp = round(S*R(2)*.5);
        x0= (screenSize(3)/2)-(extendLeft/2);
        y0 = (screenSize(4)/2)-(extendUp/2);
        Rect = [x0 y0 extendLeft extendUp];        
end


GraphicsHandle = figure('units', 'pixels', 'Tag', Tag, 'Position', Rect, 'MenuBar', 'figure', 'color', FIG_BG_COLOR, ...
                        'DefaultTextColor', FIG_FG_COLOR, 'DefaultAxesColor', AX_COLOR, 'DefaultTextInterpreter', 'none', ...
                        'DefaultAxesYColor', 'k', 'DefaultAxesZColor', 'k', 'DefaultPatchFaceColor', 'k', 'DefaultPatchEdgeColor', ...
                        'k', 'DefaultSurfaceEdgeColor', 'k', 'DefaultLineColor', 'k', 'DefaultUicontrolInterruptible', 'on', ...                                                     
                        'Renderer', 'zbuffer', 'RendererMode' , 'manual', 'Visible', visible, 'Name', sprintf(TitleFig), ...
                        'paperorientation', paperOrient, 'InvertHardcopy', 'off', 'resize', 'off', 'PaperPositionMode', 'auto');       
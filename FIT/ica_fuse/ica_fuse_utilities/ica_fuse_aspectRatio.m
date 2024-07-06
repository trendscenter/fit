function [aspectRatio] = ica_fuse_aspectRatio 
% get the aspect ratio of the figures
%
% Input:
% 
% Output:
% aspectRatio - aspect ratio calculated using the defaults from
% ica_fuse_defaults

% load defaults
ica_fuse_defaults;

global SCREENSIZE;
global WSCREEN;

dimForSquareFigure = [50, 50, 853, 800];
S0   = SCREENSIZE;
WS = WSCREEN;

xDiff = (dimForSquareFigure(3)*WS(3));
yDiff = (dimForSquareFigure(4)*WS(4));
xRatio = xDiff/yDiff;
yRatio = yDiff/xDiff;

if xRatio > 1
    xRatio = 1;
else
    yRatio = 1;
end

% aspect ratio
aspectRatio = [xRatio, yRatio];
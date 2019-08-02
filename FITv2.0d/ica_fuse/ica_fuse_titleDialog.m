function ica_fuse_titleDialog
% Title dialog displays the authors, organization and the collaborators

ica_fuse_defaults;

% Font colors
global FIG_BG_COLOR;
global FIG_FG_COLOR;
global BUTTON_BG_COLOR;
global BUTTON_FG_COLOR;

global AX_COLOR;

% Font name, units and size
global UI_FONT_NAME;
global UI_FONT_UNITS;
global UI_FONT_SIZE;

keystr1 = sprintf(['keyVal = get(gcbf, ''currentcharacter''); \n']);

keystr2 = sprintf(['if keyVal == 13 \n delete(gcbf); clear keyVal; \n end\n']);

keypressCallback = [keystr1 keystr2];

% set the defaults for the figure window
figProperty = {'Resize','off', 'windowstyle', 'modal', ...
    'MenuBar', 'none',...
    'DefaultTextColor', FIG_FG_COLOR,...
    'DefaultTextInterpreter', 'none',...
    'DefaultAxesColor', AX_COLOR,...
    'DefaultAxesXColor', 'k',...
    'DefaultAxesYColor', 'k',...
    'DefaultAxesZColor', 'k',...
    'DefaultPatchFaceColor', 'k',...
    'DefaultPatchEdgeColor', 'k',...
    'DefaultSurfaceEdgeColor', 'k',...
    'DefaultLineColor', 'k',...
    'DefaultUicontrolInterruptible', 'on',...
    'PaperType', 'usletter', ...
    'PaperUnits', 'normalized', ...
    'PaperPositionmode', 'auto', ...
    'InvertHardcopy', 'off',...
    'Renderer', 'zbuffer',...
    'color', FIG_BG_COLOR, 'resize', 'off', ...
    'keypressfcn', keypressCallback,  'NumberTitle', 'off'};

h = figure('name', 'About FIT', figProperty{1:length(figProperty)});

pos = get(h, 'position');

screenSize = get(0, 'screensize');
figurePos(1) = 0.5*(screenSize(1) + screenSize(3)) - 0.5*pos(3);
figurePos(2) = 0.5*(screenSize(2) + screenSize(4)) - 0.5*pos(4);
figurePos(3) = 0.95*pos(3);
figurePos(4) =  0.7*pos(4);

set(h, 'position', figurePos);

% positions of variables
axisPos = [0 0 1 1];
titlePos = [0.5 0.9];
normalPos = [0.04 0.6];
okPos(3) = 0.12; okPos(4) = 0.08;
okPos(1) = 0.75 - 0.5*okPos(3); okPos(2) = 0.1;
moreInfoPos(3) = 0.15; moreInfoPos(4) = 0.08;
moreInfoPos(1) = 0.25 - 0.5*moreInfoPos(3); moreInfoPos(2) = okPos(2);

% title color
titleColor = [0 0.9 0.9];

% fonts
titleFont = 16;
subtitleFont = 13;
normalFont = 11;

% set axis handle to off
axisHandle = axes('Parent', h, 'Position', axisPos, 'Visible', 'off');


% Name of the toolbox
text('units', 'normalized', 'string', 'Fusion ICA Toolbox', 'position', titlePos, 'fontsize', titleFont, 'HorizontalAlignment', 'center', ...
    'fontweight', 'bold', 'FontName', UI_FONT_NAME, 'color', titleColor, 'FontAngle', 'italic', 'parent', axisHandle);

% change the position of the text
titlePos(2) = titlePos(2) - 0.1;

% Caption of the toolbox
text('units', 'normalized', 'string', 'FITv2.0d', 'position', titlePos, 'fontsize', subtitleFont, 'HorizontalAlignment', 'center', ...
    'fontweight', 'normal', 'FontName', UI_FONT_NAME, 'parent', axisHandle);

% Display the remaining things like the version number, organization,
% authors
fusionInfo = {'', '', ['\bfRelease Date: \rm', num2str('29-Aug-2017')], '', char('\bfAuthors: \rm\bf The FIT Development Team'), ...
    '', '\bfOrganization: \rmThe MIND Research Network', '', '\bfWebsite: \rmhttp://icatb.sourceforge.net'};

% change the position of the text
normalPos(2) = titlePos(2) - 0.25;


text('units', 'normalized', 'string', fusionInfo, 'position', normalPos, 'fontsize', normalFont, 'HorizontalAlignment', 'left', ...
    'fontweight', 'normal', 'FontName', UI_FONT_NAME, 'interpreter', 'tex', 'parent', axisHandle);

% Callback for the ok button
old_pushString = 'OK';
okCallback = 'close(gcbf)';
ok_pushbutton = ica_fuse_uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', old_pushString, 'callback', okCallback, ...
    'fontunits', UI_FONT_UNITS, 'fontname', UI_FONT_NAME, 'FontSize', UI_FONT_SIZE, 'position', okPos, ...
    'backgroundcolor', BUTTON_BG_COLOR, 'ForegroundColor', BUTTON_FG_COLOR);

clear old_pushString;

old_pushString = 'More Info';

% More information give the info about the collaborators
moreInfo_pushbutton = ica_fuse_uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', ...
    old_pushString, 'callback', {@fusionMoreInfo, h, figProperty}, 'position', ...
    moreInfoPos, 'fontunits', UI_FONT_UNITS, 'fontname', UI_FONT_NAME, 'FontSize', UI_FONT_SIZE, 'ForegroundColor', ...
    BUTTON_FG_COLOR, 'backgroundcolor', BUTTON_BG_COLOR);


% click on push button to display more information about fusion
function fusionMoreInfo(hObject, eventdata, handles, figProperty)

ica_fuse_defaults;

% Font colors
global FIG_BG_COLOR;
global FIG_FG_COLOR;
global BUTTON_BG_COLOR;
global BUTTON_FG_COLOR;

global AX_COLOR;

% Font name, units and size
global UI_FONT_NAME;
global UI_FONT_UNITS;
global UI_FONT_SIZE;


% set the defaults
figHandle = figure('name', 'More Info about FIT', figProperty{1:length(figProperty)});

pos = get(figHandle, 'position');

screenSize = get(0, 'screensize');
figurePos(1) = 0.5*(screenSize(1) + screenSize(3)) - 0.5*pos(3);
figurePos(2) = 0.5*(screenSize(2) + screenSize(4)) - 0.5*pos(4);
figurePos(3) = 0.95*pos(3);
figurePos(4) =  0.7*pos(4);

set(figHandle, 'position', figurePos);

% set axis off
axisHandle = axes('Parent', figHandle, 'Position', [0 0 1 1], 'Visible', 'off');

normalFont = 10;

% Display the credits

% Add changes here
D(1).string = 'Fusion ICA toolbox includes joint ICA, CCA + joint ICA and parallel ICA methods for fusion analysis.';
D(size(D, 2) + 1).string = '';
D(size(D, 2) + 1).string = 'Additional contributions were provided by:';
D(size(D, 2) + 1).string = '1. Kent Kiehl at The MIND Research Network, Albuquerque, NM.';

%% end for changes
D(size(D,2)+1).string = '';
D(size(D,2)+1).string = 'Many additional ICA algorithms were generously contributed by Andrzej Cichocki.';
D(size(D,2)+1).string = '';
D(size(D,2)+1).string = 'ICA algorithms are also available in the ICALab toolbox at http://www.bsp.brain.riken.jp/ICALAB/ICALABImageProc/.';
D(size(D,2)+1).string = '';
D(size(D,2)+1).string = 'More information on ICA can be found in the book "Adaptive Blind Signal and Image Processing" by Andrzej Cichocki and Shun-ichi Amari.';
D(size(D,2)+1).string = '';
D(size(D,2)+1).string = 'FIT is released under the GNU General Public License. FIT uses some functions from the SPM and GIFT.';
D(size(D,2)+1).string = '';
D(size(D,2)+1).string = 'SPM: FIT uses spm_vol family of functions. Please visit SPM at http://www.fil.ion.ucl.ac.uk/spm/.';
D(size(D,2)+1).string = 'GIFT: FIT uses file selection windows and dialog boxes from the GIFT.';
D(size(D,2)+1).string = '';
D(size(D,2)+1).string = 'For more information please contact Vince Calhoun at vcalhoun@unm.edu or visit FIT project website http://icatb.sourceforge.net/.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = {D.string};

position = get(figHandle, 'position');

xOffSet = 0.02; yOffSet = 0.04;

okPos(3) = 0.12; okPos(4) = 0.08;
okPos(1) = 0.5 - 0.5*okPos(3); okPos(2) = yOffSet;


% define the position of the listbox
position(1) = xOffSet; position(2) = okPos(2) + okPos(4) + yOffSet; position(3) = 1 - 2*xOffSet; position(4) = 1 - position(2) - yOffSet;


handle_scroll = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style','listbox', ...
    'position', position, 'string', str, 'foregroundcolor', FIG_FG_COLOR, ...
    'horizontalalignment','left', 'backgroundcolor', BUTTON_BG_COLOR, 'FontSize', normalFont, 'fontunits', ...
    UI_FONT_UNITS, 'fontname', UI_FONT_NAME);

% Apply conditions for dialog box differently for different platforms
if ispc
    set(handle_scroll, 'enable', 'inactive');
else
    set(handle_scroll, 'enable', 'on');
end


%[newString, newPos] = textwrap(handle_scroll, str);
maxChars = 75;

% wrap the string inside the uicontrol
[newString, newPos] = textwrap(handle_scroll, str, maxChars);

set(handle_scroll, 'String', newString);

set(handle_scroll, 'min', 0, 'max', 2);

% make no selection
set(handle_scroll, 'value', []);

cancelCallback = 'delete(gcbf)';

old_pushString = 'Return';

ok_pushbutton = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', 'string', ...
    old_pushString, 'callback', cancelCallback, 'position', okPos, 'fontunits', UI_FONT_UNITS, ...
    'fontname', UI_FONT_NAME, 'FontSize', UI_FONT_SIZE, 'backgroundcolor', BUTTON_BG_COLOR, 'ForegroundColor', BUTTON_FG_COLOR);


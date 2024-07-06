function ica_fuse_parallel_ICAToolbox
% Parallel ICA fusion toolbox

ica_fuse_defaults;
global FIG_BG_COLOR;
global UI_FONT_NAME;

displayTag = 'parallel_ica_gui';

%%%%%% Delete any previous figures of display GUI %%%%%
displayH = findobj('tag', displayTag);

for ii = 1:length(displayH)
    delete(displayH(ii));
end
%%%%%% end for deleting previous figures of display GUI %%%%

titleText = 'ParallelICA_v1.0_alpha';

% Parallel ICA fusion figure
[displayHandle] = ica_fuse_getGraphics(titleText, 'normal', displayTag, 'on');
set(displayHandle, 'menubar', 'none');
displayPos = get(displayHandle, 'position');
displayPos(3:4) = [400 300];
set(displayHandle, 'position', displayPos);
ica_fuse_movegui(displayHandle, 'center');

set(displayHandle, 'NumberTitle', 'off');

%%%%%%%%% Menus %%%%%%%%%

% File -> New -> Open -> Close
file_menu = uimenu('parent', displayHandle, 'label', 'File');
new_menu = uimenu('parent', file_menu, 'label', 'New', 'callback', {@setupAnalysisCallback, displayHandle});
open_menu = uimenu('parent', file_menu, 'label', 'Open', 'callback', {@openCallback, displayHandle});
close_menu = uimenu('parent', file_menu, 'label', 'Close', 'callback', {@closeCallback, displayHandle});

% Tools -> Display
tools_menu = uimenu('parent', displayHandle, 'label', 'Tools');
run_menu = uimenu('parent', tools_menu, 'label', 'Run', 'callback', {@runAnalysisCallback, displayHandle});
display_menu = uimenu('parent', tools_menu, 'label', 'Display', 'callback', {@displayCallback, displayHandle});



titleColor = [0 0.9 0.9];
titleFont = 13;
xOffset = 0.02; yOffset = 0.04;
frameWidth = 1 - 2*xOffset; frameHeight = (1 - 2*yOffset);

buttonWidth = 0.3; buttonHeight = 0.06;

textHeight = 0.05; textWidth = 0.4;
frame1Pos = [xOffset, 1 - yOffset - frameHeight, frameWidth, frameHeight];
% User Interface Control
frame1H = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'frame', 'position', ...
    frame1Pos, 'backgroundcolor', FIG_BG_COLOR);

% title color
analysisStr = {titleText};
text1Pos = [0.5 - 0.5*textWidth, 1 - yOffset - textHeight, textWidth, textHeight];
text1H = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', 'position', ...
    text1Pos, 'backgroundcolor', FIG_BG_COLOR, 'foregroundcolor', titleColor, 'string', analysisStr, ...
    'horizontalalignment', 'center');
[newStr, newPos] = textwrap(text1H, analysisStr);

text1Pos(4) = newPos(4);
set(text1H, 'string', newStr, 'position', newPos);


buttonSpacing = (frame1Pos(2) + frame1Pos(4)) / 3;

% Draw Setup Analysis Display buttons
setupPos = [frame1Pos(1) + 2*xOffset, text1Pos(2) - buttonSpacing, buttonWidth, ...
        buttonHeight];
runPos = [frame1Pos(1) + frame1Pos(3) - 2*xOffset - buttonWidth, text1Pos(2) - buttonSpacing, ...
        buttonWidth, buttonHeight];
displayPos = [0.5 - 0.5*buttonWidth, frame1Pos(2) + buttonSpacing, buttonWidth, buttonHeight];

setupAnalysisH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', setupPos, 'tag', 'setup_analysis', 'string', 'Setup Analysis', 'horizontalalignment', 'center', ...
    'callback', {@setupAnalysisCallback, displayHandle});

addValue = 0.02;

% Set position for setup analysis button
setPosition(setupAnalysisH, addValue);

runAnalysisH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', runPos, 'tag', 'run_analysis', 'string', 'Run Analysis', 'horizontalalignment', 'center', 'callback', ...
    {@runAnalysisCallback, displayHandle});

% Set position for run analysis button
setPosition(runAnalysisH, addValue);

displayH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', displayPos, 'tag', 'display', 'string', 'Display', 'horizontalalignment', 'center', 'callback', ...
    {@displayCallback, displayHandle});

% Set position for display button
setPosition(displayH, addValue);
displayPos = get(displayH, 'position');
displayPos(1) = 0.5 - 0.5*displayPos(3);
set(displayH, 'position', displayPos);



%%%%%%%%%%% Define Function Callbacks %%%%%%%%

function setupAnalysisCallback(handleObj, event_data, handles)
% Setup Analysis

ica_fuse_setup_parallelICA;


function runAnalysisCallback(handleObj, event_data, handles)
% Run Analysis

ica_fuse_run_parallelICA;


function closeCallback(handleObj, event_data, handles);
% Close figure

delete(handles);

function openCallback(handleObj, event_data, handles)
% Open fMRI gene fusion file

ica_fuse_defaults;
global PARALLEL_ICA_INFO_MAT_FILE ;


fusionFile = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
    ['*', PARALLEL_ICA_INFO_MAT_FILE , '*.mat'], 'title', 'Select parallel ICA information file');      

load(fusionFile);

if ~exist('paraICAInfo', 'var')
    error(['Selected file: ', fusionFile, ' is not a valid fmri-gene fusion info file']);
end

% open fusion file
ica_fuse_setup_parallelICA([], fusionFile);


function displayCallback(handleObj, event_data, handles)
% Display Callback

ica_fuse_parallelICA_displayGUI;


function setPosition(hObject, addValue)
% Set object position

objPos = get(hObject, 'position');
% Set position
extentPos = get(hObject, 'extent');
objPos(3) = extentPos(3) + addValue; objPos(4) = extentPos(4) + addValue;
set(hObject, 'position', objPos);
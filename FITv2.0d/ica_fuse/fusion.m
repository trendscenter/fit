function fusion(fusionType)
% Fusion ICA Toolbox

ica_fuse_defaults;

% Colors
global FIG_BG_COLOR;
global FIG_FG_COLOR;

ica_fuse_delete_gui({'jointICA_fusion', 'parallel_ica_toolbox', 'CCA_ICA_fusion'});

try
    % add this statement to fix the button color
    feature('JavaFigures', 0);
catch
end

if isempty(which('ica_fuse_run_analysis.m'))
    addpath(genpath(fileparts(which('fusion.m'))), '-end');
end

if ~isappdata(0, 'FileRootsData')
    % Use MATLAB Server to get drives on Windows when matlab -nojvm option
    % is used
    ica_fuse_getdrives;
    % Close MATLAB Server
    ica_fuse_closeMatlabServer;
end

fitVerNum = 'FITv2.0d';

checkFusionPath;

figTag = 'fusion_ica';

% delete a previous figure of group ica
checkGUI = findobj('tag', figTag);

if ~isempty(checkGUI)
    for ii = 1:length(checkGUI)
        delete(checkGUI(ii));
    end
end

figHandle = figure('units', 'pixels', 'name', ['Fusion ICA Toolbox (', fitVerNum, ')'], 'color', FIG_BG_COLOR, 'position', [100 100 560 420], 'menubar', 'none', ...
    'resize', 'off', 'numbertitle', 'off', 'visible', 'off', 'tag', figTag);

movegui(figHandle, 'center');

set(figHandle, 'units', 'normalized');

xOffset = 0.04; yOffset = 0.04;

% Frame 1
frame1Pos = [0.02 0.02 0.96 0.96];
frame1H = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'frame', 'position', frame1Pos, ...
    'backgroundcolor', FIG_BG_COLOR);

% Frame 2
frameWidth = 0.8;
frameHeight = 0.2;
frameXPos = (1 - frameWidth) / 2;
frameYPos = (1 - frameHeight - 1.5*yOffset);
frame2Pos = [frameXPos frameYPos frameWidth frameHeight];
frame2H = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'frame', 'position', frame2Pos, ...
    'backgroundcolor', FIG_BG_COLOR);


textPos = frame2Pos;
textPos(4) = 0.05;
textPos(2) = textPos(2) + 1.5*textPos(4);
textPos(3) = 0.6;
textPos(1) = 0.5*(1 - textPos(3));

% Text Position
titleTextH = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'text', 'position', textPos, ...
    'backgroundcolor', FIG_BG_COLOR, 'foregroundcolor', [0 1 1], 'string', cellstr(str2mat('Fusion ICA Toolbox (FIT)', fitVerNum)), ...
    'fontsize', 13, 'max', 2);

ica_fuse_wrapStaticText(titleTextH);
newPosition = get(titleTextH, 'position');
textPos(2) = frame2Pos(2) + 0.5*newPosition(4);
textPos(4) = newPosition(4);
set(titleTextH, 'position', textPos);

% Frame 3
frameHeight = 0.4;
frameYPos = frameYPos - 1.5*yOffset - frameHeight;
frame3Pos = [frameXPos frameYPos frameWidth frameHeight];

frame3H = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'frame', 'position', frame3Pos, ...
    'backgroundcolor', FIG_BG_COLOR);

% Methods Text
textPos = frame3Pos;
textWidth = 0.2; textHeight = 0.05;
textPos(1) = 0.5 - 0.5*textWidth;
textPos(2) = textPos(2) + textPos(4) - 0.5*textHeight;
textPos(3) = textWidth;
textPos(4) = textHeight;

ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'text', 'position', textPos, ...
    'backgroundcolor', FIG_BG_COLOR, 'foregroundcolor', [0 1 1], 'string', 'Methods');

% Button width and height
buttonWidth = 0.16; buttonHeight = 0.06;
buttonYOrigin = frame3Pos(2) + 0.5*frame3Pos(4) - 0.5*buttonHeight;
buttonXOrigin = frame3Pos(1) + 2*xOffset;

%% Parallel ICA
paraICAButtonPos = [buttonXOrigin, buttonYOrigin, buttonWidth, buttonHeight];
paraICAButtonH = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', paraICAButtonPos, 'string', 'Parallel ICA', 'tag', 'paraICA_button', 'callback', {@paraICACallback, figHandle}, ...
    'tooltipstring', 'Open parallel ICA toolbox ...');

extentPos = get(paraICAButtonH, 'extent');
paraICAButtonPos(3) = extentPos(3) + 0.01;
paraICAButtonPos(4) = extentPos(4) + 0.01;
set(paraICAButtonH, 'position', paraICAButtonPos);

%% Joint ICA
buttonXOrigin = paraICAButtonPos(1) + paraICAButtonPos(3) + 2*xOffset; %frame3Pos(1) + frame3Pos(3) - 2*xOffset - buttonWidth;
jointICAButtonPos = [buttonXOrigin, buttonYOrigin, buttonWidth, buttonHeight];

jointICAButtonH = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', jointICAButtonPos, 'string', 'Joint ICA', 'tag', 'jointICA_button', 'callback', {@jointICACallback, figHandle}, ...
    'tooltipstring', 'Open joint ICA toolbox ...');

extentPos = get(jointICAButtonH, 'extent');
jointICAButtonPos(3) = extentPos(3) + 0.01;
jointICAButtonPos(4) = extentPos(4) + 0.01;

set(jointICAButtonH, 'position', jointICAButtonPos);


%% tIVA
buttonXOrigin = jointICAButtonPos(1) + jointICAButtonPos(3) + 2*xOffset; %frame3Pos(1) + frame3Pos(3) - 2*xOffset - buttonWidth;
tIVAButtonPos = [buttonXOrigin, buttonYOrigin, buttonWidth, buttonHeight];

tIVAButtonH = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', tIVAButtonPos, 'string', 'tIVA', 'tag', 'tiva_button', 'callback', {@tIVACallback, figHandle}, ...
    'tooltipstring', 'Open tIVA ...');

extentPos = get(tIVAButtonH, 'extent');
tIVAButtonPos(3) = extentPos(3) + 0.01;
tIVAButtonPos(4) = extentPos(4) + 0.01;

set(tIVAButtonH, 'position', tIVAButtonPos);


%% CCA + Joint ICA
%buttonXOrigin = paraICAButtonPos(1) + paraICAButtonPos(3) + 0.5*(jointICAButtonPos(1) -paraICAButtonPos(1) - paraICAButtonPos(3)) - 0.5*buttonWidth;
buttonXOrigin = frame3Pos(1) + 2*xOffset;
CCAJointICAButtonPos = [buttonXOrigin, buttonYOrigin - 1.25*yOffset - 0.5*buttonHeight, buttonWidth, buttonHeight];
CCAICAButtonH = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', CCAJointICAButtonPos, 'string', 'CCA + Joint ICA', 'tag', 'CCA_ICA_button', 'callback', {@CCAICACallback, figHandle}, ...
    'tooltipstring', 'Open CCA + Joint ICA toolbox ...');
extentPos = get(CCAICAButtonH, 'extent');
CCAJointICAButtonPos(3) = extentPos(3) + 0.01;
CCAJointICAButtonPos(4) = extentPos(4) + 0.01;

%CCAJointICAButtonPos(1) = buttonXOrigin;
%CCAJointICAButtonPos(1) = frame3Pos(1) + 0.5*frame3Pos(3) - 0.5*CCAJointICAButtonPos(3);
set(CCAICAButtonH, 'position', CCAJointICAButtonPos);

%% MCCA
buttonXOrigin = CCAJointICAButtonPos(1) + CCAJointICAButtonPos(3) + 2*xOffset; %frame3Pos(1) + frame3Pos(3) - 2*xOffset - buttonWidth;
mccaButtonPos = [buttonXOrigin, CCAJointICAButtonPos(2), buttonWidth, buttonHeight];

mccaButtonH = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', mccaButtonPos, 'string', 'MCCA', 'tag', 'mcca_button', 'callback', {@mccaCallback, figHandle}, ...
    'tooltipstring', 'Open MCCA ...');

extentPos = get(mccaButtonH, 'extent');
mccaButtonPos(3) = extentPos(3) + 0.01;
mccaButtonPos(4) = extentPos(4) + 0.01;

set(mccaButtonH, 'position', mccaButtonPos);

align([mccaButtonH, jointICAButtonH], 'center', 'None')

% Frame 4
frameHeight = 0.2;
frameYPos = yOffset;
frame4Pos = [frame3Pos(1) frameYPos frame3Pos(3) frameHeight];

frame4H = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'frame', 'position', frame4Pos, ...
    'backgroundcolor', FIG_BG_COLOR);

% About button
buttonWidth = 0.12; buttonHeight = 0.06;
aboutButtonPos = [frame4Pos(1) + xOffset, frame4Pos(2) + 0.5*frame4Pos(4) - 0.5*buttonHeight, buttonWidth, buttonHeight];
aboutButtonH = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', aboutButtonPos, 'string', 'About', 'tag', 'about_button', 'callback', 'ica_fuse_titleDialog', ...
    'tooltipstring', 'Information about authors and collaborators');


% Help button
helpButtonPos = [0.5 - 0.5*buttonWidth, frame4Pos(2) + 0.5*frame4Pos(4) - 0.5*buttonHeight, buttonWidth, buttonHeight];
helpButtonH = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', helpButtonPos, 'string', 'Help', 'tag', 'help_button', 'callback', 'ica_fuse_openHelp', ...
    'tooltipstring', 'Open FIT help manual');


exitButtonPos = [0.75 - 0.5*buttonWidth, frame4Pos(2) + 0.5*frame4Pos(4) - 0.5*buttonHeight, buttonWidth, buttonHeight];
exitButtonH = ica_fuse_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', exitButtonPos, 'string', 'Exit', 'tag', 'exit_button', 'callback', {@exitCallback, figHandle});

if exist('fusionType', 'var')
    if strcmpi(fusionType, 'jointica')
        % Open joint ica
        jointICACallback(jointICAButtonH, [], figHandle);
    elseif strcmpi(fusionType, 'paraica')
        % Open parallel ica
        paraICACallback(paraICAButtonH, [], figHandle);
    elseif strcmpi(fusionType, 'ccaica')
        % Open CCA + Joint ICA
        CCAICACallback(CCAICAButtonH, [], figHandle);
    elseif strcmpi(fusionType, 'mcca')
        % Open mcca
        mccaCallback(mccaButtonH, [], figHandle);
    elseif strcmpi(fusionType, 'tiva')
        % Open tiva
        tIVACallback(mccaButtonH, [], figHandle);
    elseif strcmpi(fusionType, 'exit') || strcmpi(fusionType, 'quit')
        % quit
        exitCallback(exitButtonH, [], figHandle);
    elseif strcmpi(fusionType, 'help')
        ica_fuse_openHelp;
    elseif strcmpi(fusionType, 'ver')
        disp('.........................................................................');
        disp('............... Fusion ICA Toolbox version v2.0d ............................');
        disp('.........................................................................');
        disp('Fusion ICA Toolbox: v2.0d ');
        disp('Joint ICA Toolbox: v2.0a');
        disp('Parallel ICA Toolbox: v1.0b');
        disp('CCA + Joint ICA Toolbox: v1.0b');
        disp('MCCA: v1.0a');
        disp('tIVA: v1.0a');
        fprintf('\n\n');
    else
        delete(figHandle);
        error('Unknown modality specified');
    end
    
else
    set(figHandle, 'visible', 'on');
end

%%%%%%%%%%%%%% Function Callbacks %%%%%%%%%%%%%%%

function paraICACallback(hObject, event_data, handles)
% Parallel ICA callback
% Purpose: Opens parallel ICA toolbox

set(handles, 'pointer', 'watch');

disp('Opening parallel ICA fusion ...');

fprintf('\n');

parallelICA_fusion;

set(handles, 'pointer', 'arrow');

delete(handles);


function jointICACallback(hObject, event_data, handles)
% jointICA callback
% Purpose: Opens jointICA toolbox

set(handles, 'pointer', 'watch');

disp('Opening joint ICA fusion ...');

fprintf('\n');

jointICA_fusion;

set(handles, 'pointer', 'arrow');

delete(handles);


function CCAICACallback(hObject, event_data, handles)
% CCA + Joint ICA callback
% Purpose: Opens CCA + jointICA toolbox

set(handles, 'pointer', 'watch');

disp('Opening CCA + Joint ICA fusion ...');

fprintf('\n');

CCA_ICA_fusion;

set(handles, 'pointer', 'arrow');

delete(handles);


function mccaCallback(hObject, event_data, handles)
% MCCA callback
global DATA_REDUCTION_TYPE;
global ICA_ALGORITHM_NAME;

set(handles, 'pointer', 'watch');

disp('Opening MCCA ...');

fprintf('\n');

DATA_REDUCTION_TYPE = 'MCCA';
ICA_ALGORITHM_NAME = 'None';

mCCA_fusion;

set(handles, 'pointer', 'arrow');

delete(handles);


function tIVACallback(hObject, event_data, handles)
% tIVA callback
global ICA_ALGORITHM_NAME;

set(handles, 'pointer', 'watch');

disp('Opening tIVA ...');

fprintf('\n');

ICA_ALGORITHM_NAME = 'IVA-G';

tIVA_fusion;

set(handles, 'pointer', 'arrow');

delete(handles);

function exitCallback(hObject, event_data, handles)
% Exit callback

disp('Quitting fusion');
try
    if (exist('fusion_finish.m', 'file') == 2) | (exist('fusion_finish.m', 'file') == 6)
        fprintf( 'Executing fusion_finish ...\n' );
        % execute fusion finish file
        fusion_finish;
    else
        delete(get(0, 'children'));
    end
    % end for checking finish file
    
catch
    ica_fuse_displayErrorMsg;
end



function checkFusionPath
% Check fusion path

% check fusion start up file
if (exist('fusion_startup.m', 'file') == 2) | (exist('fusion_startup.m', 'file') == 6)
    fprintf( 'Executing fusion_startup ...\n' );
    % execute fusion start up file
    fusion_startup;
    
else
    
    % Get the full file path of the fusion
    fusionPath = which('fusion.m');
    
    % Folder location of the fusion
    fusionPath = fileparts(fusionPath);
    
    % all directories on path
    pathstr = path;
    
    % Get the directories on path
    if ispc
        allDirs = strread(pathstr, '%s', 'delimiter', ';');
    else
        allDirs = strread(pathstr, '%s', 'delimiter', ':');
    end
    
    if ~isempty(allDirs)
        
        [indices] = regexp(allDirs, 'ica_fuse$');
        indices = good_cells(indices);
        matchedDirs = allDirs(indices);
        
        if length(matchedDirs) > 1
            error('Fix MATLAB path such that it has only one version of fusion at a time');
        elseif length(matchedDirs) == 1
            if strcmpi(fusionPath, matchedDirs{1})
                return;
            else
                error('Fix MATLAB path such that it has only one version of fusion at a time');
            end
        end
        
    end
    
end
% end for adding path

function ind = good_cells( mycell )
% Find good cells

if ~iscell(mycell)
    mycell = {mycell};
end

ind = cellfun('isempty', mycell);

% Good cells
ind = (ind == 0)';



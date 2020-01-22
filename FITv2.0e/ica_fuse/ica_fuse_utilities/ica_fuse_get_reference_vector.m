function [reference] = ica_fuse_get_reference_vector(groupNames, numSubjects)
% Get reference vector based on group information

% Load defaults
ica_fuse_defaults;

global HELP_FG_COLOR;
global UI_FONT_NAME;
global UI_FONT_SIZE;

groupTag = 'Group_Info';

% Close all figures with the tag Group_Info
checkHandles = findobj(0, 'tag', groupTag);
if ~isempty(checkHandles)
    delete(checkHandles);
end

titleFig = 'Enter reference vector for groups';
% Plot Figure
[graphicsHandle] = ica_fuse_getGraphics(titleFig, 'normal', groupTag, 'off');

if iscell(groupNames)
    groupNames = str2mat(groupNames);
end

handles_data.groupNames = groupNames;
handles_data.numSubjects = numSubjects;

set(graphicsHandle, 'menubar', 'none', 'userdata', handles_data);

if ispc
    set(graphicsHandle, 'windowstyle', 'modal');
end


%%%%%%%%%%%%% Draw Title here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title color
titleColor = [0 0.9 0.9];

% fonts
titleFont = 13;

axisH = axes('Parent', graphicsHandle, 'position', [0 0 1 1], 'visible', 'off');

xPos = 0.5; yPos = 0.97;

text(xPos, yPos, titleFig, 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', 'fontsize', titleFont, 'HorizontalAlignment', ...
    'center', 'FontName', UI_FONT_NAME, 'parent', axisH);

%%%%%%%%%%%%% end for drawing title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Offsets
xOffset = 0.05; yOffset = 0.05;

% List width and height
listWidth = 0.42; listHeight = 0.42;
listTextHeight = 0.05;
listOrigin = 0.3;

listTextPos = [xOffset, listOrigin + listHeight + yOffset - 0.5*listTextHeight, listWidth, listTextHeight];

% Plot listbox text
listTextH = ica_fuse_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', listTextPos, 'String', 'Groups', 'fontsize', UI_FONT_SIZE - 1, 'tag', 'groupListTxt');
ica_fuse_wrapStaticText(listTextH);
listTextPos = get(listTextH, 'position');
listTextPos(2) = listOrigin + listHeight + yOffset;
set(listTextH, 'position', listTextPos);

listPos = [xOffset, listOrigin, listWidth, listHeight];
% Plot listbox
listH = ica_fuse_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listPos, 'String', groupNames, 'fontsize', UI_FONT_SIZE - 1, 'min', 0, 'max', 1, ...
    'value', 1, 'callback', {@listCallback, graphicsHandle}, 'tag', 'groupList', 'userdata', pwd);

% Plot editbox text
editTextPos = [listTextPos(1) + listTextPos(3) + xOffset, listOrigin + listHeight + yOffset - 0.5*listTextHeight, listWidth, listTextHeight];
editTextH = ica_fuse_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', editTextPos, 'String', ['Enter vector for ', deblank(groupNames(1, :)), ' of length ', num2str(numSubjects(1))], ...
    'fontsize', UI_FONT_SIZE - 1, 'tag', 'groupText');
ica_fuse_wrapStaticText(editTextH);
editTextPos = get(editTextH, 'position');
editTextPos(2) = listOrigin + listHeight + yOffset;
set(editTextH, 'position', editTextPos);

editPos = editTextPos;
editPos(2) = listPos(2);
editPos(4) = listPos(4);

% Define the context menu
cmenu = uicontextmenu;

drawnow;

editH = ica_fuse_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', '', 'fontsize', UI_FONT_SIZE - 1, 'max', 2, 'min', 0, 'tag', 'groupVector', ...
    'callback', {@editCallback, graphicsHandle}, 'horizontalalignment', 'left', 'tooltipstring', 'You can also enter vector using ascii file when you right click on edit box.');

set(editH, 'uicontextmenu', cmenu);

% Define the context menu items
item1 = uimenu(cmenu, 'Label', 'Load Ascii File', 'callback', {@editContextMenuCallback, editH, listH});

% Plot Ok and cancel buttons
buttonWidth = 0.12; buttonHeight = 0.05;
okPos = [0.5 - 0.5*buttonWidth, 0.1, buttonWidth, buttonHeight];
cancelPos = okPos;
cancelPos(1) = 0.25 - 0.5*buttonWidth;

helpPos = cancelPos;
helpPos(1) = 0.75 - 0.5*buttonWidth;

cancelH = ica_fuse_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', cancelPos, 'string', 'Cancel', 'callback', 'delete(gcbf)');

referenceData = repmat(struct('data', ''), length(numSubjects), 1);
okH = ica_fuse_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', okPos, 'string', 'OK', 'userdata', referenceData, 'tag', 'groupInfoOk', 'callback', ...
    {@okCallback, graphicsHandle});

helpH = ica_fuse_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', helpPos, 'string', '?', 'foregroundcolor', HELP_FG_COLOR, 'callback', {@groupInfoCallback});

reference = [];

try
    set(graphicsHandle, 'visible', 'on');
    waitfor(graphicsHandle);
catch
end

if isappdata(0, 'groupInfoData')
    reference = getappdata(0, 'groupInfoData');
    rmappdata(0, 'groupInfoData');
end

function listCallback(hObject, event_data, handles)
% Listbox callback

% Get group info
handles_data = get(handles, 'userdata');
groupNames = handles_data.groupNames;
numSubjects = handles_data.numSubjects;

% Current entry in listbox
getVal = get(hObject, 'value');
getStr = get(hObject, 'string');
selectedStr = deblank(getStr(getVal, :));

okH = findobj(handles, 'tag', 'groupInfoOk');
okData = get(okH, 'userdata');

% Update group info text
groupTextH = findobj(handles, 'tag', 'groupText');
set(groupTextH, 'string', ['Enter vector for ', deblank(groupNames(getVal, :)), ' of length ', num2str(numSubjects(getVal))]);
ica_fuse_wrapStaticText(groupTextH);

groupVectorH = findobj(handles, 'tag', 'groupVector');
set(groupVectorH , 'string', okData(getVal).data);


function editCallback(hObject, event_data, handles)
% Editbox callback

groupListH = findobj(handles, 'tag', 'groupList');
% Group list value
groupListVal = get(groupListH, 'value');

okH = findobj(handles, 'tag', 'groupInfoOk');
okData = get(okH, 'userdata');


getStr = get(hObject, 'string');

try
    getVal = str2num(getStr);
    getVal = getVal(:)';
catch
    error(lasterr);
end

% Set data
if groupListVal
    okData(groupListVal).data = num2str(getVal);
    set(okH, 'userdata', okData);
end


function groupInfoCallback(hObject, event_data, handles)
% Open help dialog box about groups

msgStr = ['You have the option to enter reference vector for each group. This information will be used when calculating PCA.', ...
        ' Click on each entry in left listbox to enter the corresponding reference vector for that group. Option is provided to enter group information using ascii file when you use right click on edit box.'];

ica_fuse_dialogBox('title', 'Reference Vector', 'textBody', msgStr, 'textType', 'large');


function okCallback(hObject, event_data, handles)
% Ok callback

okData = get(hObject, 'userdata');

% Get group info
handles_data = get(handles, 'userdata');
groupNames = handles_data.groupNames;
numSubjects = handles_data.numSubjects;

% Do error check
for nS = 1:length(numSubjects)
    if isempty(okData(nS).data)
        error(['Reference vector for group ', deblank(groupNames(nS, :)), ' is not entered']);
    end
    
    temp = str2num(okData(nS).data);
    
    if (length(temp) ~= numSubjects(nS))
        error('Error:GroupInfo', 'Length of reference vector (%s) for group %s \n doesn''t match the number of subjects (%s) in that group', ...
            num2str(length(temp)), deblank(groupNames(nS, :)), num2str(numSubjects(nS)));
    end
    
    okData(nS).data = temp(:)';
end
% End for error checking

setappdata(0, 'groupInfoData', [okData.data]');

delete(handles);


function editContextMenuCallback(hObject, event_data, editH, listH)
% Context menu callback

getStr = get(listH, 'string');
getVal = get(listH, 'value');
startPath = get(listH, 'userdata');

selectedStr = deblank(getStr(getVal, :));

asciiFile = ica_fuse_selectEntry('title', ['Select ascii file for group ', selectedStr], 'filter', '*.asc', 'typeEntity', 'file', ...
    'typeSelection', 'single', 'startPath', startPath);

if ~isempty(asciiFile)
    set(listH, 'userdata', fileparts(asciiFile));
    temp = load(asciiFile, '-ascii');
    
    if ~isempty(temp)
        if isnumeric(temp)
            temp = num2str(temp);
        end
        set(editH, 'string', temp);
        ica_fuse_executeCallback(editH);
    end
    
end
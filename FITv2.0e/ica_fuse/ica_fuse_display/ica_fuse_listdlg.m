function [selection, value] = ica_fuse_listdlg(varargin)
% inputs: prompt string, List string, name of the ok button, name of the
% cancel button, window style, name of the figure window, selection mode,
% list size in pixels, help button data, movegui option

% output: selection and value associated with the click of ok or cancel
% button


% get defaults
titleFig = '';
smode = 2;
promptString = {'Select Entries'};
liststring = [];


%listsize = [300 400];
okString = 'Ok';
cancelString = 'Cancel';
moveguiOption = 'center';
maxSelection = [];
windowStyle = 'normal';

if mod(length(varargin), 2) ~= 0
    % input arguments must be in pairs
    error('Arguments to ica_fuse_listdlg must in pairs.')
end

% loop over number of parameters
for ii = 1:2:length(varargin)
    switch lower(varargin{ii})
        case 'promptstring'
            promptString = varargin{ii+1};
        case 'liststring'
            liststring = varargin{ii+1};
        case 'okstring'
            okString = varargin{ii+1};
        case 'cancelstring'
            cancelString = varargin{ii+1};
        case 'windowstyle'
            windowStyle = varargin{ii+1};
        case 'title_fig'
            titleFig = varargin{ii+1};
        case 'selectionmode'
            switch lower(varargin{ii+1})
                case 'single'
                    smode = 1;
                case 'multiple'
                    smode = 2;
            end
        case 'listsize'
            listsize = varargin{ii+1};
        case 'maxselection'
            maxSelection = varargin{ii+1};
        case 'movegui'
            moveguiOption = varargin{ii+1};
        case 'help'
            helpDetails = varargin{ii+1};
        otherwise
            error(['Unknown parameter name passed to ica_fuse_listdlg.  Name was ' varargin{ii}])
    end
end

% Setup figure for GUI
[InputHandle] = ica_fuse_getGraphics( titleFig, 'normal',  titleFig, 'on');
set(InputHandle, 'menubar', 'none');

figPos = get(InputHandle, 'position');

if ~exist('listsize', 'var')    
    listsize(1) = 0.75*figPos(3); listsize(2) = 0.85*figPos(4);
end

figPos = [figPos(1) figPos(2) listsize];

set(InputHandle, 'position', figPos);

% move the figure to the specified position
ica_fuse_movegui(InputHandle, moveguiOption);

% set up positions for listbox, ok , cancel, help and select all buttons

% offsets
xOffset = 0.05; yOffset = 0.03;

% Ok button position
okWidth = 0.2; okHeight = 0.06;

if exist('helpDetails', 'var')
    okPos = [0.5 - 0.5*okWidth yOffset + 0.5*okHeight okWidth okHeight];
else
    okPos = [0.75 - 0.5*okWidth yOffset + 0.5*okHeight okWidth okHeight];
end

% cancel button position
cancelWidth = 0.2; cancelHeight = 0.06;
cancelPos = [0.25 - 0.5*cancelWidth yOffset + 0.5*cancelHeight cancelWidth cancelHeight];

% Help button position
if exist('helpDetails', 'var')
    helpWidth = 0.2; helpHeight = 0.06;
    helpPos = [0.75 - 0.5*helpWidth yOffset + 0.5*helpHeight helpWidth helpHeight];
end

% select all button position
if smode ~= 1
    selectAllWidth = 0.4; selectAllHeight = 0.06;
    selectAllPos = [0.5 - 0.5*selectAllWidth okPos(2) + 0.5*okPos(4) + yOffset + 0.5*selectAllHeight selectAllWidth selectAllHeight];
end

% text position for prompt String
textWidth = 0.9; textHeight = 0.1;
textPos = [0.5 - 0.5*textWidth 0.9 - 0.5*textHeight textWidth textHeight];

%%%%%%%%%%%%%%%%%%%%%%%%% Plot UI Controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prompt string text
promptTextH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', textPos, 'horizontalalignment', 'center', 'string', promptString);

if ~iscell( promptString)
    promptString =  {promptString};
end

[promptString, newPos] = textwrap(promptTextH, promptString);
textPos(4) = newPos(4);
set(promptTextH, 'string', promptString, 'position', textPos);

% listbox Position
if smode ~= 1
    listOrigin = selectAllPos(2) + selectAllPos(4) + yOffset;
else
    listOrigin = okPos(2) + okPos(4) + yOffset;
end

listHeight = (textPos(2) - listOrigin - yOffset);
listWidth = 0.9;
listPos = [0.5 - 0.5*listWidth listOrigin listWidth listHeight];

% prompt string text
listH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listPos, ...
                                       'horizontalalignment', 'center', 'string', liststring, 'callback', {@doListboxClick}, ...
                                       'max', smode, 'value', 1);

% Draw Ok and cancel buttons here
% Draw push button
okData.maxSelection = maxSelection;

okH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', okString, ...
                                       'tooltipstring', 'Select Entries...', 'tag', 'ok', 'userdata', okData);

cancelH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', cancelPos, 'string', cancelString,  ...
                                            'tooltipstring', 'Cancel Entries...', 'tag', 'cancel');

% plot help button if necessary
if exist('helpDetails', 'var')
    % Help push button
    helpH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', helpPos, 'string', '?', ...
                                             'tooltipstring', 'help', 'tag', 'help', 'callback', @helpCallback, 'userdata', helpDetails);
end

% plot help button if necessary
if smode == 2
    % Help push button
    selectAllH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', selectAllPos, 'string', 'Select-All', ...
                                                    'tooltipstring', 'Select all', 'tag', 'select-all', 'callback', {@selectAllCallback, listH});
    set(listH, 'callback', {@doListboxClick, selectAllH});
end


% set the callbacks here
set(okH, 'callback', {@okCallback, listH});
set(cancelH, 'callback', {@cancelCallback, listH});

try
    set(InputHandle, 'visible','on');
    %set(fig, 'color', BG_COLOR);
    uiwait(InputHandle);
catch
    if ishandle(InputHandle)
        delete(InputHandle);
    end
end

if isappdata(0, 'ListDialogAppData')
    ad = getappdata(0, 'ListDialogAppData');
    selection = ad.selection;
    value = ad.value;
    rmappdata(0,'ListDialogAppData')
else
    % figure was deleted
    selection = [];
    value = 0;
end


% define function callbacks here

% help callback
function helpCallback(hObject, evd, handles)

% help button callback
helpDetails = get(hObject, 'userdata');

if isfield(helpDetails, 'str')
    str = helpDetails.str;
else
    str = '';
end

if isfield(helpDetails, 'title')
    titleString = helpDetails.title;
else
    titleString = '';
end

figHandle = ica_fuse_dialogBox('title', titleString, 'textBody', str, 'textType', 'large');

waitfor(figHandle);

% select all callback
function selectAllCallback(selectall_btn, evd, listbox)

set(selectall_btn, 'enable', 'off');
getString = get(listbox,'string');
if iscell(getString)
    sizeString = length(getString);
else
    sizeString = size(getString, 1);
end
set(listbox, 'value', 1:sizeString);

% Ok callback
function okCallback(ok_btn, evd, listbox)

ad.value = 1;
ad.selection = get(listbox,'value');
okData = get(ok_btn, 'userdata');

maxSelection = [];

if isstruct(okData)
    if isfield(okData, 'maxSelection')
        maxSelection = okData.maxSelection;
    end
end

% check if the entries are selected than the maximum specified
try
    % check if max selection is not empty
    if ~isempty(maxSelection)        
        % check if selection is more than the specified maximum selection
        if length(ad.selection) > maxSelection
            set(listbox, 'value', []);
            errorMsg = ['Maximum selection allowed is ', num2str(maxSelection)];
            ica_fuse_errorDialog(errorMsg);        
        else            
            ad.selection = get(listbox, 'value');
            setappdata(0, 'ListDialogAppData', ad)
            delete(gcbf);
        end 
    else
        ad.selection = get(listbox, 'value');
        setappdata(0, 'ListDialogAppData', ad)
        delete(gcbf);
    end
catch
    ica_fuse_errorDialog(lasterr); 
end


% cancel callback
function cancelCallback(cancel_btn, evd, listbox)

ad.value = 0;
ad.selection = [];
setappdata(0,'ListDialogAppData', ad)
delete(gcbf);

function doListboxClick(listbox, evd, selectall_btn)

% if this is a doubleclick, doOK
if strcmp(get(gcbf, 'SelectionType'),'open')
    okCallback([],[],listbox);
elseif nargin == 3
    if length(get(listbox, 'string')) == length(get(listbox, 'value'))
        set(selectall_btn, 'enable', 'off')
    else
        set(selectall_btn, 'enable', 'on')
    end
end
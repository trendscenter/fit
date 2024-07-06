function Answer = ica_fuse_inputDialog(varargin)
% input dialog box - Dialog box for getting the user input.
% The uicontrols are pop up and edit box. 
% 
% Input:
% 1. inputText - structure containg the type of uicontrol,
%    prompt string, optional arguments (default answers).
% 2. Title - title of the input dialog
%
% Output:
%
% Answer - cell array containing the answers selected by the user,
% returns {} when the cancel button is pressed
% 
% Note:
% units for the uicontrols are normalized
% returns the values for the variables can be string or number (By default
% answer is string)

handle_visibility = 'on';

if nargin < 2
    error('Atleast two arguments are required');
end

if mod(nargin, 2) ~= 0
    error('Arguments must be in pairs');
end

for i = 1:2:nargin
    if strcmp('inputtext', lower(varargin{i}))
        inputText = varargin{i + 1};
    elseif strcmp('title', lower(varargin{i}))
        titleFig = varargin{i + 1};
    elseif strcmp('handle_visibility', lower(varargin{i}))
        handle_visibility = varargin{i + 1};
    end
end

handle_visibility = lower(handle_visibility);

if ~exist('titleFig', 'var')
    titleFig = 'Input Dialog';
end

% number of UIcontrols excluding the push buttons
numUIcontrols = length(inputText);

if ~isfield(inputText, 'enable')
    % loop over input texts
    for ii = 1:length(inputText)
        inputText(ii).enable = 'on';
    end
end
% end for checking enable property

okTag = 'OK';
cancelTag = 'Cancel';

[figHandle] = ica_fuse_plot_controls_fig(inputText, titleFig, handle_visibility, okTag, cancelTag);

okHandle = findobj(figHandle, 'tag', okTag);
cancelHandle = findobj(figHandle, 'tag', cancelTag);

% done callback
set(okHandle, 'callback', {@okCallback, figHandle});

% cancel callback
set(cancelHandle, 'callback', {@cancelCallback, figHandle});


% set the figure data
set(figHandle, 'UserData', inputText, 'keypressfcn', {@keypressCallback, figHandle, okHandle});

if strcmpi(handle_visibility, 'on')
    
    % wait for figure handle
    try
        waitfor(figHandle);
    catch
        delete(figHandle);
    end
    
else    
    % execute the done callback
    okCallback(okHandle, [], figHandle);
end

Answer = {};
% check application data
if isappdata(0, 'inputAppData')
    % get the application data
    Answer = getappdata(0, 'inputAppData');
    rmappdata(0, 'inputAppData'); % remove the application data
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Define function callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ok callback
function okCallback(handleObj, event_data, handles)

% OK push button callback for the input dialog
% get the figure data
inputParameters = get(handles, 'userdata');
answer = {};
% loop over number of parameters
for ii = 1:length(inputParameters)
    tag = [inputParameters(ii).tag];
    answerH = findobj(handles, 'tag', tag); % find tag
    if isfield(inputParameters(ii), 'answerType')
        answerType = inputParameters(ii).answerType; % answer type
    else
        answerType = 'string';
    end
    getStyle = inputParameters(ii).uiType;
    getString = get(answerH, 'string');  % get the answer string
    if iscell(getString)
        getString = str2mat(getString);
    end
    getValue = get(answerH, 'value'); % get the answer value
    % check the uicontrol type
    if strcmp(lower(getStyle), 'edit')
        answerString = deblank(getString);
    else
        answerString = deblank(getString(getValue, :));
    end
    % check the answer type
    if strcmp(lower(answerType), 'numeric')
        answer{ii} = str2num(answerString);
    else
        answer{ii} = answerString;
    end
end
% set the application data
setappdata(0, 'inputAppData', answer);
delete(handles);

% cancel callback
function cancelCallback(handleObj, evd, handles)

% close figure window
delete(handles);

function keypressCallback(hObject, event_data, handles, okHandle);

keyVal = get(handles, 'currentcharacter');

% if enter key is detected
if keyVal == 13
    okCallback(okHandle, [], handles);
end
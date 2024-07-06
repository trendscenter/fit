function [InputHandle] = ica_fuse_plot_controls_fig(inputText, figureTag, handle_visibility, okTag, cancelTag, plotHelp)
% Plot on a figure the user interface controls with Done and Cancel buttons
% Function callbacks can be defined outside of the function

if ~exist('handle_visibility', 'var')
    handle_visibility = 'on';
end

if ~exist('figureTag', 'var')
    figureTag = 'input_dialog';
end

if ~exist('okTag', 'var')
    okTag = 'done';
end

if ~exist('cancelTag', 'var')
    cancelTag = 'cancel';
end

% Plot Help If Necessary
if ~exist('plotHelp', 'var')
    plotHelp = 0;
end

ica_fuse_defaults;
global HELP_FG_COLOR;
global UI_FONT_SIZE;

promptPrefix = 'prompt';

% Setup figure for GUI
[InputHandle] = ica_fuse_getGraphics(figureTag, 'normal', figureTag, handle_visibility);

% Set no menu bar for the figure
set(InputHandle, 'Menubar', 'none');

% plot all the user interface controls and their options to the right
% plot a help button to the right

% offsets for x and y
xOffset = 0.02; yOffset = 0.05;

% UI control width and height
uicontrol_Width = 0.3; uicontrol_Height = 0.05;

% number of UIcontrols excluding the push buttons
numUIcontrols = length(inputText);

% Define the push button positions

% Ok Button push button position
okPos(3) = 0.15; okPos(4) = 0.05;
okPos = [0.75 - 0.5*okPos(3) 0.94 - 0.5*okPos(4) okPos(3) okPos(4)];

% Cancel button position
cancelPos = [0.25 - 0.5*okPos(3) okPos(2) okPos(3) okPos(4)];

% set the yPos
yPos = okPos(2) - 1.5*yOffset;

% Position of prompt string
promptPos = [0.02 yPos  0.5 0.1];

% Position of answer string
answerPos = promptPos;
answerPos(1) = promptPos(1) + promptPos(3) + xOffset; answerPos(3) = uicontrol_Width; answerPos(4) = uicontrol_Height;
helpPos = answerPos;
helpPos(1) = helpPos(1) + helpPos(3) + xOffset;
helpPos(3) = 0.06;


%%%%%%%%%%% plot all the uicontrols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countTag = 0;
for ii = 1:numUIcontrols
    if isfield(inputText, 'width')
        answerPos(3) = inputText(ii).width;
    end
    % increment count
    countTag = countTag + 1;
    tempH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', ...
        promptPos, 'string', inputText(ii).promptString, 'tag', [promptPrefix, inputText(ii).tag], 'panel_tag', ...
        ['panel_text', inputText(ii).tag], 'fontsize', UI_FONT_SIZE - 1);
    
    % get the last handle of tempH
    promptHandle = tempH(end);
    
    if ~iscell(inputText(ii).promptString)
        newTextString = {inputText(ii).promptString};
    end
    
    % wrap the prompt string and get the new position
    [newTextString, newPos] = textwrap(promptHandle, newTextString);
    
    % select the same height as the text wrap
    promptPos(4) = newPos(4);
    
    promptPos(2) = newPos(2) - 0.5*promptPos(4);
    answerPos(2) = newPos(2) - 0.5*answerPos(4);
    
    
    % include panel only for scrolling
    storeTag{countTag} = get(tempH(1), 'tag');
    
    % store all the initial positions
    initialYPositions(countTag) = promptPos(2);
    
    set(promptHandle, 'string', newTextString);
    
    set(tempH(1), 'position', promptPos);    
    
    clear tempH;
    
    % store all the initial positions
    countTag = countTag + 1;
    
    tempH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', inputText(ii).uiType, 'position', ...
        answerPos, 'string', inputText(ii).answerString, 'tag', inputText(ii).tag, ...
        'enable', inputText(ii).enable, 'panel_tag', ['panel', inputText(ii).tag], 'value', inputText(ii).value, ...
        'fontsize', UI_FONT_SIZE - 1, 'HorizontalAlignment', 'center');        
    
    storeTag{countTag} = get(tempH(1), 'tag');
    
    initialYPositions(countTag) = answerPos(2);
    
    if plotHelp
        
        % Store all the initial positions
        helpPos(2) = newPos(2) - 0.5*helpPos(4);
        countTag = countTag + 1;    
        helpH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
            helpPos, 'string', '?', 'tag', ['help_', inputText(ii).tag], 'enable', 'on', 'fontsize', UI_FONT_SIZE - 1, ...
            'HorizontalAlignment', 'center', 'foregroundcolor', HELP_FG_COLOR);            
        storeTag{countTag} = get(helpH(1), 'tag');
        
        initialYPositions(countTag) = helpPos(2);
        
        % Update help button position
        helpPos(2) = helpPos(2) - 0.5*helpPos(4) - yOffset;
        
    end
    
    promptPos(2) = promptPos(2) - 0.5*promptPos(4) - yOffset;
    answerPos(2) = answerPos(2) - 0.5*answerPos(4) - yOffset;    
    
end
%%%%%%%%%%% End for plotting all the uicontrols %%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% Plot push buttons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define push buttons
cancelHandle = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
    cancelPos, 'string', cancelTag, 'Tag', cancelTag);

% define push buttons
okHandle = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
    okPos, 'string', okTag, 'Tag', okTag);

checkPosition = promptPos(2);

% define slider data
sliderData.tag = storeTag;
sliderData.maxHeight = okPos(2); % - 0.5*okPos(4);
sliderData.minHeight = checkPosition - yOffset;
sliderData.initialYPositions = initialYPositions;

% Plot slider if necessary
if initialYPositions(end) < 0.01
    
    sliderPos = [0.965 0 0.04 1];
    
    % Plot slider
    maxVal = 0;
    minVal = -abs(sliderData.minHeight);
    % control the slider step size
    sliderStep = [0.05 0.2];
    
    sliderH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'slider', 'position', ...
        sliderPos, 'min', minVal, 'max', maxVal, 'value', maxVal, 'callback', {@verSliderCallback, InputHandle}, ...
        'userdata', sliderData, 'sliderstep', sliderStep);       
    
end

%%%%%%%%%%%%%%% Define callbacks for the controls %%%%%%%%%%%%%

function verSliderCallback(handleObj, evd, handles)
% vertical slider callback

% execute the slider callback
ica_fuse_verSliderCallback(handleObj, handles);
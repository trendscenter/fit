function [maskFiles] = ica_fuse_selectMask(featureNames, modalities, dataInfo)
% Select mask for each modaility

ica_fuse_defaults;
global BUTTON_BG_COLOR;
global UI_FONT_NAME;
global UI_FONT_UNITS;
global UI_FONT_SIZE;

oldDir = pwd;

titleFig = 'Select mask for features';

% Setup figure for GUI
[InputHandle] = ica_fuse_getGraphics(titleFig, 'normal', titleFig);

set(InputHandle, 'menubar', 'none');

%%%%%%%%%%%%% Draw Title here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title color
titleColor = [0 0.9 0.9];

% fonts
titleFont = 13;

axisH = axes('Parent', InputHandle, 'position', [0 0 1 1], 'visible', 'off');

xPos = 0.5; yPos = 0.97;

text(xPos, yPos, titleFig, 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', 'fontsize', ...
    titleFont, 'HorizontalAlignment', 'center', 'FontName', UI_FONT_NAME, 'parent', axisH);

%%%%%%%%%%%%% end for drawing title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% draw listbox on the left handside
xOffset = 0.04; yOffset = 0.05;

% draw a frame here
frameOrigin = yOffset;
frameWidth = 0.65; frameHeight = 0.06;
framePos = [0.5 - 0.5*frameWidth frameOrigin frameWidth frameHeight];

% plot frame
%frameH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'frame', 'position', framePos);

% selected Text Position here
selectedTextPos = framePos;
selectedTextPos(2) = framePos(2) + 0.005; selectedTextPos(4) = framePos(4) - 0.01;
selectedTextPos(1) = framePos(1) + 0.005; selectedTextPos(3) = framePos(3) - 0.01;

% plot selected text here
% selectedTextH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', ...
%     selectedTextPos, 'string', 'Files Selected: ', 'tag', 'files_selected_text_box');

% Listbox positions
listWidth = 0.65;
listPos = framePos(2) + framePos(4) + yOffset;
listHeight = (yPos - listPos - 2*yOffset);
listPos = [xOffset listPos listWidth listHeight];

listData.modalities = modalities;
listData.indices = repmat(struct('ind', ''), 1, length(modalities));

% plot the listbox here
listH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listPos, ...
    'string', featureNames, 'tag', 'data_listbox', 'userdata', listData, 'callback', ...
    {@listCallback, InputHandle});
clear subject_string;

% calculate the positioning of the pushbuttons
pushButtonWidth = 0.2; pushButtonHeight = 0.05;
pushButtonXPos = listPos(1) + listPos(3) + 1.5*xOffset;
pushButtonYPos = listPos(2) + listPos(4) - yOffset;
% push button position
pushPos1 = [pushButtonXPos pushButtonYPos pushButtonWidth pushButtonHeight];
pushPos2 = pushPos1; pushPos2(2) = pushPos2(2) - pushPos2(4) - yOffset;
pushPos3 = pushPos2; pushPos3(2) = pushPos3(2) - pushPos3(4) - yOffset;
pushPos4 = pushPos3; pushPos4(2) = pushPos4(2) - pushPos4(4) - yOffset;
pushPos5 = pushPos4; pushPos5(2) = pushPos5(2) - pushPos5(4) - yOffset;

viewH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'position', pushPos1, 'style', 'pushbutton', ...
    'string', 'View', 'tag', 'view', 'callback', {@viewFilesCallback, listH});

changeH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'position', pushPos2, 'style', 'pushbutton', ...
    'string', 'Change', 'tag', 'change', 'callback', {@changeCallback, listH});

helpString = sprintf(['Click on the left listbox and select the files for each dataset.', ...
    ' After the selection is done, press Ok button.', ' The selected files for that ' ...
    'data set can be viewed by clicking View button. The selected files can be changed by clicking on change button.']);
helpData.string = helpString; helpData.title =  'Data Selection';

helpH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'position', pushPos3, 'style', 'pushbutton', ...
    'BackgroundColor', BUTTON_BG_COLOR, 'ForegroundColor', [1 1 0], 'fontunits', UI_FONT_UNITS, 'fontname', ...
    UI_FONT_NAME, 'FontSize', UI_FONT_SIZE, 'string', '?', 'tag', 'help', 'fontweight', 'bold', 'userdata', helpData, ...
    'callback', {@helpCallback, InputHandle});

okH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'position', pushPos4, 'style', 'pushbutton', ...
    'string', 'OK', 'tag', 'OK', 'userdata', cell(1, length(featureNames)), 'callback', ...
    {@okCallback, InputHandle});

figureData = struct('changeFiles', 'no', 'dataInfo', dataInfo, 'startPath', pwd);

set(InputHandle, 'userdata', figureData);


try
    set(InputHandle, 'visible','on');
    waitfor(InputHandle);
catch
    if ishandle(InputHandle)
        delete(InputHandle);
    end
end

% get the application data
if isappdata(0, 'okAppData')
    maskFiles = getappdata(0, 'okAppData');
    % remove application data
    rmappdata(0, 'okAppData');
    cd(oldDir);
else
    cd(oldDir);
    maskFiles = {};
    error('Figure Window was quit');
end



%%%%%%%%%%%% Object Callbacks %%%%%%%%%%%%%%%%

function listCallback(handleObj, event_data, handles)
% Listbox callback

ica_fuse_defaults;
global MRI_DATA_FILTER;


try

    figureData = get(handles, 'userdata');

    dataInfo = figureData.dataInfo;

    startPath = figureData.startPath;

    % get the value
    getValue = get(handleObj, 'value');

    featureNames = get(handleObj, 'string');

    listData = get(handleObj, 'userdata');

    modalities = listData.modalities;

    files = str2mat(dataInfo(1).feature(getValue).files.name);

    files = deblank(files(1, :));

    if (strcmpi(modalities{getValue}, 'eeg') || strcmpi(modalities{getValue}, 'gene') || strcmpi(modalities{getValue}, 'behavioral'))
        data = ica_fuse_loadData(files);
        if isempty(listData.indices(getValue).ind)
            listData.indices(getValue).ind = ['1:', num2str(size(data, 1))];
        end
    end

    okHandle = findobj(handles, 'tag', 'OK');

    okData = get(okHandle, 'userdata');

    figureData = get(handles, 'userdata');

    if isempty(okData)
        okData = cell(1, length(featureNames));
    end


    if isempty(okData{getValue}) || strcmpi(figureData.changeFiles, 'yes')
        if (strcmpi(modalities{getValue}, 'eeg') || strcmpi(modalities{getValue}, 'gene') || strcmpi(modalities{getValue}, 'behavioral'))
            prompt = {['Enter ', modalities{getValue}, ' indices to include for feature: ', featureNames{getValue}]};
            def = {listData.indices(getValue).ind};
            dlgTitle = [upper(modalities{getValue}), ' Mask']; %'EEG Mask';
            lineNo = 1;
            answer = ica_fuse_inputdlg2(prompt, dlgTitle, lineNo, def);
            if ~isempty(answer)
                val = str2num(answer{1});
                if max(val) > size(data, 1)
                    error(['Indices exceed the maximum limit']);
                elseif min(val) < 1
                    error(['Indices need to be positive integers']);
                end
                okData{getValue} = answer{:};
                listData.indices(getValue).ind = answer{:};
            end
        elseif strcmpi(modalities{getValue}, 'fmri') | strcmpi(modalities{getValue}, 'smri')
            [P] = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', MRI_DATA_FILTER, ...
                'title', ['Select mask for feature ', featureNames{getValue}], 'startPath', startPath);
            okData{getValue} = P;
            if ~isempty(P)
                figureData.startPath = fileparts(deblank(P(1, :)));
            end

        end
    end

    figureData.changeFiles = 'no';

    set(handleObj, 'userdata', listData);

    % Set the user data for OK button
    set(okHandle, 'userdata', okData);

    set(handles, 'userdata', figureData);

catch

    disp(lasterr);
    ica_fuse_errorDialog(lasterr, 'Error Found', 'modal');

end


function okCallback(hObject, event_data, handles)

listH = findobj(handles, 'tag', 'data_listbox');

modalities = get(listH, 'userdata');
featureNames = get(listH, 'string');

okData = get(hObject, 'userdata');

try

    for nn = 1:length(featureNames)

        currentFeatureName = featureNames{nn};

        if isempty(okData{nn})
            error(['Please select the mask for feature ', currentFeatureName]);
        end

    end

    setappdata(0, 'okAppData', okData);
    delete(handles);

catch
    disp(lasterr);
    ica_fuse_errorDialog(lasterr, 'Mask selection', 'modal');
end



function changeCallback(hObject, event_data, handles)
% change the selected files

figH = get(hObject, 'parent'); %get the figure handle
figureData = get(figH, 'userdata'); % get the figure data
figureData.changeFiles = 'yes'; % pass the flag to change the files
set(figH, 'userdata', figureData); % set the figure data

listCallback(handles, [], figH); % use listbox callback to select the files


% 3. View Files callback
function viewFilesCallback(hObject, event_data, handles)

% View files
figH = get(hObject, 'parent'); %get the figure handle
okHandle = findobj(figH, 'tag', 'OK');

okData = get(okHandle, 'userdata');

listValue = get(handles, 'value');
listString = get(handles, 'string');
selectedString = listString{listValue};
textBody = okData{listValue};

if isempty(textBody)
    textBody = '';
end

titleFig = ['Files for ', selectedString];
% open dialog box
figHandle = ica_fuse_dialogBox('title', titleFig, 'textBody', textBody, 'textType', 'large');
waitfor(figHandle);

% 4. Help callback
function helpCallback(hObject, event_data, handles)

getUserdata = get(hObject, 'userdata');

% open dialog box
figHandle = ica_fuse_dialogBox('title', getUserdata.title, 'textBody', getUserdata.string, 'textType', 'large');

waitfor(figHandle);
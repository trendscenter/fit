function [controlPara, displayParameters] = ica_fuse_get_display_parameters(displayParameters, controlPara, ...
    handle_visibility)
% Plot display defaults

% Load defaults
ica_fuse_defaults;
global ANATOMICAL_FILE; % Anatomical file location

if ~exist('handle_visibility', 'var')
    handle_visibility = 'on';
end

% Plot answer string
if isempty(controlPara)

    % Set the answer string
    numPara = 1;
    controlPara(numPara).uiType = 'popup';
    controlPara(numPara).promptString = 'Convert To Z-Scores';
    controlPara(numPara).value = 1;
    controlPara(numPara).width = 0.2;
    controlPara(numPara).tag = 'convert_to_z';
    controlPara(numPara).enable = 'on';
    controlPara(numPara).answerType = 'string';
    controlPara(numPara).callback = '';

    numPara = numPara + 1;
    controlPara(numPara).uiType = 'edit';
    controlPara(numPara).promptString = 'Threshold';
    controlPara(numPara).value = 0;
    controlPara(numPara).width = 0.2;
    controlPara(numPara).tag = 'z_threshold';
    controlPara(numPara).enable = 'on';
    controlPara(numPara).answerType = 'numeric';
    controlPara(numPara).callback = '';

    numPara = numPara + 1;
    controlPara(numPara).uiType = 'popup';
    controlPara(numPara).promptString = 'Image Values';
    controlPara(numPara).value = 1;
    controlPara(numPara).width = 0.4;
    controlPara(numPara).tag = 'image_values';
    controlPara(numPara).enable = 'on';
    controlPara(numPara).answerType = 'string';
    controlPara(numPara).callback = '';

    numPara = numPara + 1;
    controlPara(numPara).uiType = 'popup';
    controlPara(numPara).promptString = 'Components per figure';
    controlPara(numPara).value = 1;
    controlPara(numPara).width = 0.2;
    controlPara(numPara).tag = 'images_per_figure';
    controlPara(numPara).enable = 'on';
    controlPara(numPara).answerType = 'numeric';
    controlPara(numPara).callback = '';

    numPara = numPara + 1;
    controlPara(numPara).uiType = 'popup';
    controlPara(numPara).promptString = 'Anatomical Plane';
    controlPara(numPara).value = 1;
    controlPara(numPara).width = 0.4;
    controlPara(numPara).tag = 'anatomical_plane';
    controlPara(numPara).enable = 'on';
    controlPara(numPara).answerType = 'string';
    controlPara(numPara).callback = '';
    %controlPara(numPara).callback = {@anatomicalPopupCallback, displayHandle};

    numPara = numPara + 1;
    controlPara(numPara).uiType = 'edit';
    controlPara(numPara).promptString = 'Slices (in mm)';
    controlPara(numPara).value = 0;
    controlPara(numPara).width = 0.4;
    controlPara(numPara).tag = 'slices_in_mm';
    controlPara(numPara).enable = 'on';
    controlPara(numPara).answerType = 'numeric';
    controlPara(numPara).callback = '';


    %%%%%%%%%%%%%%%%%%%% Set answer string %%%%%%%%%%%%%%%%%%
    % Set the answer string and their value
    for nn = 1:numPara
        % Get the answer string and value
        [answerStr, answerVal] = ica_fuse_get_displayDefaults(controlPara(nn).tag);
        controlPara(nn).answerString = answerStr;
        controlPara(nn).value = answerVal;
        if strcmpi(controlPara(nn).tag, 'anatomical_plane')
            anatomicalView = deblank(answerStr(answerVal, :));
        end
    end

    % Get slices in mm
    structVol = ica_fuse_get_vol_nifti(ANATOMICAL_FILE);
    [slices_in_mm, displayParameters] = ica_fuse_get_slices_in_mm(displayParameters, anatomicalView);


    slicesIndex = strmatch('slices_in_mm', str2mat(controlPara.tag), 'exact');
    controlPara(slicesIndex).answerString = slices_in_mm;

    %%%%%%%%%%%%%%%% End for setting answer string %%%%%%%%%%%%

end
% end for plotting answer string



okTag = 'OK';
cancelTag = 'Cancel';

[displayHandle] = ica_fuse_plot_controls_fig(controlPara, 'Display Parameters', handle_visibility, okTag, cancelTag);

anatomicalPlaneH = findobj(displayHandle, 'tag', 'anatomical_plane');
okH = findobj(displayHandle, 'tag', okTag);
cancelH = findobj(displayHandle, 'tag', cancelTag);

handles_data.controlPara = controlPara;
handles_data.displayParameters = displayParameters;

set(displayHandle, 'userdata', handles_data);

set(anatomicalPlaneH , 'callback', {@anatomicalPopupCallback, displayHandle});
set(okH , 'callback', {@okCallback, displayHandle});
set(cancelH , 'callback', {@cancelCallback, displayHandle});

if strcmpi(handle_visibility, 'off')
    okCallback(okH, [], displayHandle);
else

    try
        waitfor(displayHandle);
    catch
        delete(displayHandle);
    end

end

% Get the display parameters
if isappdata(0, 'displayDefaultsData')
    okData = getappdata(0, 'displayDefaultsData');
    controlPara = okData.controlPara;
    displayParameters = okData.displayParameters;
    rmappdata(0, 'displayDefaultsData');
end


function anatomicalPopupCallback(hObject, event_data, handles)

% Set string for slices in mm
ica_fuse_setString_slices(hObject, handles);

function okCallback(hObject, event_data, handles)
% Ok callback

try
    % Get the user data of the figure
    handles_data = get(handles, 'userdata');

    controlPara = handles_data.controlPara;
    displayParameters = handles_data.displayParameters;
    
    [displayParameters, controlPara] = ica_fuse_get_values_display_controls(handles, [], controlPara, displayParameters, 'no');  
    okData.controlPara = controlPara;
    okData.displayParameters = displayParameters;

    setappdata(0, 'displayDefaultsData', okData);
    delete(handles);

catch

    if ishandle(handles)
        delete(handles);
    end

    disp(lasterr);
end


function cancelCallback(hObject, event_data, handles)

delete(handles);
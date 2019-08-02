function [displayParameters, controlPara] = ica_fuse_get_values_display_controls(handles, complistH, controlPara, displayParameters, flagSetData)
% Get values from controls

if ~exist('controlPara', 'var')
    controlPara = get(complistH, 'userdata');
end

if ~exist('displayParameters', 'var')
    displayParameters = get(handles, 'userdata');
end

if ~exist('flagSetData', 'var')
    flagSetData = 'set_data';
end

for nn = 1:length(controlPara)
    
    currentTag = controlPara(nn).tag;
    currentHandle = findobj(handles, 'tag', currentTag);
    
    % Get the current control style
    getStyle = get(currentHandle, 'style');
    
    % Current value and string
    currentVal = get(currentHandle, 'value');
    currentStr = get(currentHandle, 'string');
    
    % Store current values and string in control para
    controlPara(nn).value = currentVal;
    controlPara(nn).answerString = currentStr;
    
    % Get the selected string
    if ica_fuse_findstr(getStyle, 'popup') | strcmpi(getStyle, 'listbox')
        currentStr = deblank(currentStr(currentVal, :));
    end
    
    if strcmpi(controlPara(nn).answerType, 'numeric')
        currentStr = str2num(currentStr);
    end
    
    if strcmpi(currentTag, 'image_values')
        [imageValStr, currentStr] = ica_fuse_imageValues(currentStr);
    end
    
    if strcmpi(currentTag, 'convert_to_z')
        if strcmpi(currentStr, 'yes')
            currentStr = 1;
        else
            currentStr = 0;
        end
    end
    
    % Set the field to the structure display parameters
    displayParameters = setfield(displayParameters, currentTag, currentStr);
    
end

if strcmpi(flagSetData, 'set_data')
    set(handles, 'userdata', displayParameters);
end
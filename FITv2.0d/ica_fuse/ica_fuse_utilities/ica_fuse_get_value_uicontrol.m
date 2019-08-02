function [answerString, value] = ica_fuse_get_value_uicontrol(handles, tag)
% Function to get the string and values from uicontrol

if ischar(tag)
    objectH = findobj(handles, 'tag', tag);
else
    objectH = tag;
end

handle_type = get(objectH, 'Type');

% Initialise output
answerString = '';
value = 0;

% Check if the handle type is uicontrol or not
if strcmpi(handle_type, 'uicontrol')
   
    % Get answer string and value
    answerString = get(objectH, 'string');
    value = get(objectH, 'value');
        
else
    
    error('Handle is not uicontrol');
    
end
function ica_fuse_setString_slices(hObject, handles)
% Set the string for slices in mm

[answerString, answerVal] = ica_fuse_get_value_uicontrol(handles, hObject);

anatomicalView = deblank(answerString(answerVal, :));

displayParameters = get(handles, 'userdata');

[slices_in_mm] = ica_fuse_get_slices_in_mm(displayParameters, anatomicalView);

slices_in_mmH = findobj(handles, 'tag', 'slices_in_mm');
set(slices_in_mmH, 'string', slices_in_mm);
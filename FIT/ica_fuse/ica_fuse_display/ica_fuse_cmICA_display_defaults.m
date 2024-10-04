function  dispParameters = ica_fuse_cmICA_display_defaults
%% Get display defaults
%

getDispDefaults;

appStr = 'disp_para_data';
if (isappdata(0, appStr))
    dispParameters = getappdata(0, appStr);
    rmappdata(0, appStr);
else
    error('Input parameters window was quit');
end


function getDispDefaults
ica_fuse_defaults;
global UI_FONT_SIZE;
global UI_FONT_NAME;
global CONVERT_TO_Z;
global IMAGE_VALUES;
global Z_THRESHOLD;
global ANATOMICAL_PLANE;

%% Draw graphics
figureTag = 'setup_image_viewer_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

% Setup figure for GUI
InputHandle = ica_fuse_getGraphics('Image viewer Params', 'displaygui', figureTag, 'on');
set(InputHandle, 'menubar', 'none');

yPos = 0.92;
promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight;
okWidth = 0.12; okHeight = promptHeight;

%% Select display type
promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;

%% Z-scores
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Convert To Z-scores?', 'tag', 'prompt_z', ...
    'fontsize', UI_FONT_SIZE - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
zOptions = char('Yes', 'No');
chkInd = strmatch(lower(CONVERT_TO_Z), lower(zOptions), 'exact');
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', zOptions, 'tag', 'convert_to_z', 'fontsize', UI_FONT_SIZE - 1, ...
    'value', chkInd);

%% Enter threshold
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter threshold', 'tag', 'prompt_threshold', 'fontsize', ...
    UI_FONT_SIZE - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
threshStr = '1.0';
editH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', threshStr, 'tag', 'threshold', 'fontsize', UI_FONT_SIZE - 1);


%% Image values
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select Image Values', 'tag', 'prompt_display', ...
    'fontsize', UI_FONT_SIZE - 1);
ica_fuse_wrapStaticText(textH);


imOptions = char('Positive and Negative', 'Positive', 'Absolute value', 'Negative');
returnValue = 2;
%returnValue = strmatch(lower(IMAGE_VALUES), lower(imOptions), 'exact');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', imOptions, 'tag', 'image_values', 'fontsize', UI_FONT_SIZE - 1, ...
    'value', returnValue);


%% Slice plane
planeOptions = {'Axial', 'Sagittal', 'Coronal'};
planeVal = 1;
%planeVal = strmatch(lower(ANATOMICAL_PLANE), lower(planeOptions));

promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select anatomical plane', 'tag', 'prompt_slice_plane', ...
    'fontsize', UI_FONT_SIZE - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', planeOptions, 'tag', ...
    'slice_plane', 'fontsize', UI_FONT_SIZE - 1, 'value', planeVal, 'callback', {@updateSlicesInMM, InputHandle});


%% Slices in mm
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter slices in mm', 'tag', 'prompt_slice_range', ...
    'fontsize', UI_FONT_SIZE - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', '', 'tag', ...
    'slices_in_mm', 'fontsize', UI_FONT_SIZE - 1, 'value', 1);


%% Load anatomical

startPath = fileparts(which('fusion.m'));
startPath = fullfile(startPath, 'ica_fuse_templates');
structFile = fullfile(startPath, 'ch2bet_3x3x3.nii');

anatWidth = 0.2;
anatHeight = 0.05;
anatPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, anatWidth, anatHeight];
anatPos(2) = anatPos(2) - 0.5*anatPos(4);
ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', anatPos, 'string', 'Select Anatomical', 'tag', 'anat_button', 'fontsize',...
    UI_FONT_SIZE - 1, 'callback', {@selectAnatomical, InputHandle}, 'userdata', structFile);


%% Add cancel, save and run buttons
okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Ok', 'tag', 'run_button', 'fontsize',...
    UI_FONT_SIZE - 1, 'callback', {@getDefs, InputHandle});

slicePlaneH = findobj(InputHandle, 'tag', 'slice_plane');

updateSlicesInMM(slicePlaneH, [], InputHandle);

waitfor(InputHandle);

function getDefs(hObject, event_data, handles)
%% get defs

% convert to z
convertToZH = findobj(handles, 'tag', 'convert_to_z');
zStr = cellstr(get(convertToZH, 'string'));
zVal = get(convertToZH, 'value');
convert_to_z = lower(zStr{zVal});
dispParameters.convert_to_z = convert_to_z;

% Threshold
dispParameters.threshold = str2num(get(findobj(handles, 'tag', 'threshold'), 'string'));

% Image values
imH = findobj(handles, 'tag', 'image_values');
imStr = cellstr(get(imH, 'string'));
imVal = get(imH, 'value');
image_values = lower(imStr{imVal});
dispParameters.image_values = image_values;

% Slice plane
spH = findobj(handles, 'tag', 'slice_plane');
spStr = cellstr(get(spH, 'string'));
spVal = get(spH, 'value');
slice_plane = lower(spStr{spVal});
dispParameters.slice_plane = slice_plane;

% Slices in mm
dispParameters.slices_in_mm = str2num(get(findobj(handles, 'tag', 'slices_in_mm'), 'string'));

% anatomical file
anatH = findobj(handles, 'tag', 'anat_button');
structFile = get(anatH, 'userdata');
dispParameters.structFile = structFile;

setappdata(0, 'disp_para_data', dispParameters);

delete(handles);

drawnow;


function updateSlicesInMM(slicePlaneH, event_data, figH)
%% Update slices in mm
%

anatH = findobj(figH, 'tag', 'anat_button');
structFile = get(anatH, 'userdata');

% slice plane
%slicePlaneH = findobj(figH, 'tag', 'slice_plane');
sliceOptions = cellstr(get(slicePlaneH, 'string'));
sliceOptionsVal = get(slicePlaneH, 'value');
slicePlane = lower(deblank(sliceOptions{sliceOptionsVal}));

drawnow;

% Compute slices in mm
imagVol = ica_fuse_spm_vol(structFile);
% get the slices in mm for the corresponding plane
[sliceParameters] = ica_fuse_get_slice_def(imagVol, slicePlane);
% get the slices in mm
slices_in_mm = sliceParameters.slices;
clear sliceParameters;
% construct string
slices_in_mm = ica_fuse_constructString(slices_in_mm);

% slices in mm
sliceInMMH = findobj(figH, 'tag', 'slices_in_mm');
set(sliceInMMH, 'string', slices_in_mm);


function selectAnatomical(hObject, event_data, figH)
%% Anatomical callback
%

startPath = fileparts(which('fusion.m'));
startPath = fullfile(startPath, 'ica_fuse_templates');

oldDir = pwd;

if (~exist(startPath, 'dir'))
    startPath = pwd;
end

% get the structural file
structFile = ica_fuse_selectEntry('typeEntity', 'file', 'title', 'Select Structural File', 'filter', ...
    '*.img;*.nii', 'fileType', 'image', 'fileNumbers', 1, 'startpath', startPath);

drawnow;

cd(oldDir);

if (~isempty(structFile))
    set(hObject, 'userdata', structFile);
    sliceInMMH = findobj(figH, 'tag', 'slice_plane');
    updateSlicesInMM(sliceInMMH, [], figH);
end
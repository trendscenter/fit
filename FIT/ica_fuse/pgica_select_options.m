function varargout = pgica_select_options(varargin)
% PGICA_SELECT_OPTIONS MATLAB code for pgica_select_options.fig
%      PGICA_SELECT_OPTIONS, by itself, creates a new PGICA_SELECT_OPTIONS or raises the existing
%      singleton*.
%
%      H = PGICA_SELECT_OPTIONS returns the handle to a new PGICA_SELECT_OPTIONS or the handle to
%      the existing singleton*.
%
%      PGICA_SELECT_OPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PGICA_SELECT_OPTIONS.M with the given input arguments.
%
%      PGICA_SELECT_OPTIONS('Property','Value',...) creates a new PGICA_SELECT_OPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pgica_select_options_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pgica_select_options_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pgica_select_options

% Last Modified by GUIDE v2.5 05-Jun-2022 23:00:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @pgica_select_options_OpeningFcn, ...
    'gui_OutputFcn',  @pgica_select_options_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pgica_select_options is made visible.
function pgica_select_options_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pgica_select_options (see VARARGIN)

% Choose default command line output for pgica_select_options
handles.output = hObject;

if (length(varargin) >= 1)
    pgicaInfo = varargin{1};
    try
        outputDir = pgicaInfo.outputDir;
    catch
        outputDir = pwd;
    end
    handles.outputDir = outputDir;
    outputPrefix = '';
    try
        outputPrefix = pgicaInfo.prefix;
    catch
    end
    try
        num_features = pgicaInfo.num_features;
    catch
        num_features = length(pgicaInfo.modalities);
    end
    pgicaInfo.num_features = num_features;
    if (~isfield(pgicaInfo, 'featuresInfo'))
        featuresInfo = repmat(struct('name', '', 'modality_type', '', 'components', '', 'files', '', 'mask', ''), 1, num_features);
    else
        featuresInfo = pgicaInfo.featuresInfo;
    end
    feature_str = cellstr([repmat('Feature ', num_features, 1), num2str((1:num_features)')]);
    feature_str{2} = [ feature_str{2}, ' (fMRI)'];
    featuresInfo(2).modality_type = 'fmri';
    set(handles.prefix, 'string', outputPrefix);
    set(handles.pgica_features, 'string', feature_str);
    set(handles.pgica_features, 'userdata', featuresInfo);
    set(handles.pgica_features, 'value', 1);
    pgica_features_Callback(handles.pgica_features, [], handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pgica_select_options wait for user response (see UIRESUME)
% uiwait(handles.pgica_feature_options);


% --- Outputs from this function are returned to the command line.
function varargout = pgica_select_options_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

waitfor(hObject);

results = [];
appName = 'pgOptsData';
if (isappdata(0, appName))
    results = getappdata(0, appName);
    rmappdata(0, appName);
end

varargout{1} = results;



% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

output_prefix = deblank(get(handles.prefix, 'string'));
featuresInfo = get(handles.pgica_features, 'userdata');
outputDir = handles.outputDir;

for i = 1:length(featuresInfo)
    if isempty(featuresInfo(i).files)
        error(['Feature ', num2str(i), ' files are not selected']);
    end
    if isempty(featuresInfo(i).components)
        error(['Feature ', num2str(i), ' components are not selected']);
    end
end

components = [featuresInfo.components];

ica_options = ica_fuse_paraICAOptions(min(components), 'on');

drawnow;

if (isempty(ica_options))
    error('Options window was quit');
end

pgicaInfo.outputDir = outputDir;
pgicaInfo.featuresInfo = featuresInfo;
pgicaInfo.ica_options = ica_options;
pgicaInfo.prefix = output_prefix;
pgicaInfo.param_file = [output_prefix, '_pgica_ica_param_info.mat'];
pgicaInfo.modalities = cellstr(char(featuresInfo.modality_type));
pgicaInfo.num_features = length(pgicaInfo.modalities);

fileName = fullfile(pgicaInfo.outputDir, pgicaInfo.param_file);
disp(['Saving PGICA fusion parameters in file ', fileName]);
save(fileName, 'pgicaInfo');
disp('Done');
fprintf('\n');

setappdata(0, 'pgOptsData', pgicaInfo);

delete(gcbf);

drawnow;


function feature_text_Callback(hObject, eventdata, handles)
% hObject    handle to feature_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of feature_text as text
%        str2double(get(hObject,'String')) returns contents of feature_text as a double

feature_name = deblank(get(hObject, 'string'));

listval = get(handles.pgica_features, 'value');
featuresInfo = get(handles.pgica_features, 'userdata');
featuresInfo(listval).name = feature_name;
set(handles.pgica_features, 'userdata', featuresInfo);

% --- Executes during object creation, after setting all properties.
function feature_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to feature_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in modality_type.
function modality_type_Callback(hObject, eventdata, handles)
% hObject    handle to modality_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modality_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modality_type

featuresInfo = get(handles.pgica_features, 'userdata');
listval = get(handles.pgica_features, 'value');
mstr = get(hObject, 'string');
mval = get(hObject, 'value');
featuresInfo(listval).modality_type = mstr{mval};
set(handles.pgica_features, 'userdata', featuresInfo);

% --- Executes during object creation, after setting all properties.
function modality_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modality_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mask_type.
function mask_type_Callback(hObject, eventdata, handles)
% hObject    handle to mask_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mask_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mask_type


featuresInfo = get(handles.pgica_features, 'userdata');
listval = get(handles.pgica_features, 'value');
val = get(hObject, 'value');
str=get(hObject, 'string');

if (~strcmpi(str(val, :), 'default mask'))
    maskFile = ica_fuse_selectEntry('typeEntity', 'file', 'title', 'Select binary mask file ...', 'typeSelection', ...
        'single', 'filter', '*.nii;*.asc;*.txt');
    
    drawnow;
    if (isempty(maskFile))
        set(hObject, 'value', 1);
    end
    featuresInfo(listval).mask = maskFile;
else
    featuresInfo(listval).mask = '';
end

set(handles.pgica_features, 'userdata', featuresInfo);

% --- Executes during object creation, after setting all properties.
function mask_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_files.
function browse_files_Callback(hObject, eventdata, handles)
% hObject    handle to browse_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

featuresInfo = get(handles.pgica_features, 'userdata');
listval = get(handles.pgica_features, 'value');
feature_name = featuresInfo(listval).name;
filterPattern = '*.asc;*.dat;*.txt';

if ~isempty(featuresInfo(listval).files)
    question_val = ica_fuse_questionDialog('title', 'Files already selected!!! Do you want to change files?', 'textbody', ...
        featuresInfo(listval).files, 'okstring', 'Change', 'cancelstring', 'Cancel');
    if (isempty(question_val) || (question_val == 0))
        return;
    end
end

modality_type = featuresInfo(listval).modality_type;
if (isempty(modality_type))
    mval = get(handles.modality_type, 'value');
    mstr = get(handles.modality_type, 'string');
    modality_type = deblank(mstr{mval});
end

featuresInfo(listval).modality_type = modality_type;

if (strcmpi(featuresInfo(listval).modality_type, 'smri') || strcmpi(featuresInfo(listval).modality_type, 'fmri') || ...
        strcmpi(featuresInfo(listval).modality_type, 'dti'))
    filterPattern = '*.nii;*.img';
end
files = ica_fuse_selectEntry('typeEntity', 'file', 'title', ['Select ', feature_name, ' files for all subjects ...'], ...
    'typeSelection', 'multiple', 'filter', filterPattern);

featuresInfo(listval).files = files;
if (~isempty(files))
    set(hObject, 'foregroundcolor', [0, 1, 0]);
else
    set(hObject, 'foregroundcolor', [1, 1, 1]);
end
set(handles.pgica_features, 'userdata', featuresInfo);


function components_Callback(hObject, eventdata, handles)
% hObject    handle to components (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of components as text
%        str2double(get(hObject,'String')) returns contents of components as a double

featuresInfo = get(handles.pgica_features, 'userdata');
listval = get(handles.pgica_features, 'value');
featuresInfo(listval).components = str2num(get(hObject, 'string'));
set(handles.pgica_features, 'userdata', featuresInfo);

% --- Executes during object creation, after setting all properties.
function components_CreateFcn(hObject, eventdata, handles)
% hObject    handle to components (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pgica_features.
function pgica_features_Callback(hObject, eventdata, handles)
% hObject    handle to pgica_features (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pgica_features contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pgica_features

liststr = get(hObject, 'string');
listval = get(hObject, 'value');
current_selection = deblank(liststr{listval});
if (strcmpi(current_selection, 'feature 2 (fmri)'))
    set(handles.modality_type, 'enable', 'off');
else
    set(handles.modality_type, 'enable', 'on');
end

featuresInfo = get(hObject, 'userdata');
try
    feature_name = featuresInfo(listval).name;
catch
    feature_name = '';
end

mask = '';
try
    mask = featuresInfo(listval).mask;
catch
end

files = '';
try
    files = featuresInfo(listval).files;
catch
end

set(handles.browse_files, 'foregroundcolor', [1, 1, 1]);
if (~isempty(files))
    set(handles.browse_files, 'foregroundcolor', [0, 1, 0]);
end

if (~isempty(mask))
    set(handles.mask_type, 'value', 2);
else
    set(handles.mask_type, 'value', 1);
end

set(handles.components, 'string', num2str(featuresInfo(listval).components));

set(handles.feature_text, 'string', feature_name);

modality_type = featuresInfo(listval).modality_type;
if (listval == 2)
    modality_type = 'fmri';
end

chk = strmatch(lower(modality_type), lower(get(handles.modality_type, 'string')), 'exact');
if (isempty(chk))
    chk = 1;
end

set(handles.modality_type, 'value', chk);

%set(handles.pgica_features, 'userdata', featuresInfo);


% --- Executes during object creation, after setting all properties.
function pgica_features_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pgica_features (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prefix_Callback(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefix as text
%        str2double(get(hObject,'String')) returns contents of prefix as a double


% --- Executes during object creation, after setting all properties.
function prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

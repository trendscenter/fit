function varargout = CCA_ICA_fusion(varargin)
% CCA_ICA_FUSION M-file for CCA_ICA_FUSION.fig
%      CCA_ICA_FUSION, by itself, creates a new CCA_ICA_FUSION or raises the existing
%      singleton*.
%
%      H = CCA_ICA_FUSION returns the handle to a new CCA_ICA_FUSION or the handle to
%      the existing singleton*.
%
%      CCA_ICA_FUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CCA_ICA_FUSION.M with the given input arguments.
%
%      CCA_ICA_FUSION('Property','Value',...) creates a new CCA_ICA_FUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CCA_ICA_FUSION_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CCA_ICA_FUSION_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CCA_ICA_FUSION

% Last Modified by GUIDE v2.5 28-Jan-2008 11:02:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CCA_ICA_fusion_OpeningFcn, ...
    'gui_OutputFcn',  @CCA_ICA_fusion_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CCA_ICA_fusion is made visible.
function CCA_ICA_fusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CCA_ICA_fusion (see VARARGIN)

% Choose default command line output for CCA_ICA_fusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CCA_ICA_fusion wait for user response (see UIRESUME)
% uiwait(handles.figure1);

movegui(hObject, 'center');

if isempty(which('ica_fuse_run_analysis.m'))
    addpath(genpath(fileparts(which('fusion.m'))), '-end');
end

% --- Outputs from this function are returned to the command line.
function varargout = CCA_ICA_fusion_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function utilities_CreateFcn(hObject, eventdata, handles)
% hObject    handle to utilities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in pushbutton1.
function setup_analysis_Callback(hObject, eventdata, handles, optional)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_defaults;
global FUSION_INFO_MAT_FILE;

if (~exist('optional', 'var'))
    optional = 'new';
end

if strcmpi(optional, 'new')
    ica_fuse_setup_analysis_cca_ica;
else        
    
    fusionFile = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        ['*', FUSION_INFO_MAT_FILE, '*.mat'], 'title', 'Select joint ICA fusion information file');
    
    load(fusionFile);
    
    if ~exist('fusionInfo', 'var')
        error(['Selected file: ', fusionFile, ' is not a valid joint ICA fusion information file']);
    end
    
    % open fusion file
    fusionFile = ica_fuse_setup_analysis_cca_ica([], fusionFile);
    
end


% --- Executes on button press in pushbutton2.
function run_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_run_analysis;

% --- Executes on button press in analysis_info.
function fusion_info_Callback(hObject, eventdata, handles)
% hObject    handle to analysis_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_fusion_info;


% --- Executes on button press in pushbutton4.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_display;

% % --- Executes on selection change in toolboxes.
function utilities_Callback(hObject, eventdata, handles)
% hObject    handle to toolboxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns toolboxes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from toolboxes


% Utilities callback
selectedStr = get(hObject, 'string');
selectedVal = get(hObject, 'value');

if iscell(selectedStr)
    % get the selected string
    selectedStr = lower(selectedStr{selectedVal});
else
    % get the selected string
    selectedStr = lower(deblank(selectedStr(selectedVal, :)));
end
ica_fuse_utilities(selectedStr);


% --- Executes on button press in about.
function about_fusion_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_titleDialog;


function fusion_help_Callback(hObject, eventdata, handles)
% Open HTML help manual

ica_fuse_openHelp;


% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));

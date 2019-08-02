function varargout = mCCA_fusion(varargin)
% MCCA_FUSION M-file for mCCA_fusion.fig
%      MCCA_FUSION, by itself, creates a new MCCA_FUSION or raises the existing
%      singleton*.
%
%      H = MCCA_FUSION returns the handle to a new MCCA_FUSION or the handle to
%      the existing singleton*.
%
%      MCCA_FUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MCCA_FUSION.M with the given input arguments.
%
%      MCCA_FUSION('Property','Value',...) creates a new MCCA_FUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before jointICA_fusion_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mCCA_fusion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mCCA_fusion

% Last Modified by GUIDE v2.5 15-May-2017 15:32:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mCCA_fusion_OpeningFcn, ...
    'gui_OutputFcn',  @mCCA_fusion_OutputFcn, ...
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


% --- Executes just before mCCA_fusion is made visible.
function mCCA_fusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mCCA_fusion (see VARARGIN)

% Choose default command line output for mCCA_fusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mCCA_fusion wait for user response (see UIRESUME)
% uiwait(handles.figure1);

movegui(hObject, 'center');

if isempty(which('ica_fuse_run_analysis.m'))
    addpath(genpath(fileparts(which('fusion.m'))), '-end');
end

% --- Outputs from this function are returned to the command line.
function varargout = mCCA_fusion_OutputFcn(hObject, eventdata, handles)
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
function setup_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_setup_analysis_mcca;


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

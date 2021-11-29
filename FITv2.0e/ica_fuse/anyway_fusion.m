function varargout = anyway_fusion(varargin)
% ANYWAY_FUSION MATLAB code for anyway_fusion.fig
%      ANYWAY_FUSION, by itself, creates a new ANYWAY_FUSION or raises the existing
%      singleton*.
%
%      H = ANYWAY_FUSION returns the handle to a new ANYWAY_FUSION or the handle to
%      the existing singleton*.
%
%      ANYWAY_FUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANYWAY_FUSION.M with the given input arguments.
%
%      ANYWAY_FUSION('Property','Value',...) creates a new ANYWAY_FUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before anyway_fusion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to anyway_fusion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help anyway_fusion

% Last Modified by GUIDE v2.5 22-Nov-2021 23:50:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @anyway_fusion_OpeningFcn, ...
    'gui_OutputFcn',  @anyway_fusion_OutputFcn, ...
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


% --- Executes just before anyway_fusion is made visible.
function anyway_fusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to anyway_fusion (see VARARGIN)

% Choose default command line output for anyway_fusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes anyway_fusion wait for user response (see UIRESUME)
% uiwait(handles.anyway_fusion);


% --- Outputs from this function are returned to the command line.
function varargout = anyway_fusion_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_anyway_fusion.
function setup_anyway_fusion_Callback(hObject, eventdata, handles)
% hObject    handle to setup_anyway_fusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dynamic_coherence_setup;
ica_fuse_anyway_setup_analysis;


% --- Executes on button press in display_anyway.
function display_anyway_Callback(hObject, eventdata, handles)
% hObject    handle to display_anyway (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_anyway_display;

% --- Executes on button press in exit_dyn_coh.
function exit_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to exit_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));


% --- Executes on button press in run_anyway_fuson.
function run_anyway_fuson_Callback(hObject, eventdata, handles)
% hObject    handle to run_anyway_fuson (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_run_anyway_fusion;

function varargout = joint_cmica_fusion(varargin)
% JOINT_CMICA_FUSION MATLAB code for joint_cmica_fusion.fig
%      JOINT_CMICA_FUSION, by itself, creates a new JOINT_CMICA_FUSION or raises the existing
%      singleton*.
%
%      H = JOINT_CMICA_FUSION returns the handle to a new JOINT_CMICA_FUSION or the handle to
%      the existing singleton*.
%
%      JOINT_CMICA_FUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JOINT_CMICA_FUSION.M with the given input arguments.
%
%      JOINT_CMICA_FUSION('Property','Value',...) creates a new JOINT_CMICA_FUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before joint_cmica_fusion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to joint_cmica_fusion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help joint_cmica_fusion

% Last Modified by GUIDE v2.5 01-Jul-2024 21:11:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @joint_cmica_fusion_OpeningFcn, ...
    'gui_OutputFcn',  @joint_cmica_fusion_OutputFcn, ...
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


% --- Executes just before joint_cmica_fusion is made visible.
function joint_cmica_fusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to joint_cmica_fusion (see VARARGIN)

% Choose default command line output for joint_cmica_fusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes joint_cmica_fusion wait for user response (see UIRESUME)
% uiwait(handles.joint_cmica_fusion);


% --- Outputs from this function are returned to the command line.
function varargout = joint_cmica_fusion_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_joint_cmica.
function setup_joint_cmica_Callback(hObject, eventdata, handles)
% hObject    handle to setup_joint_cmica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_setup_joint_cmica;


% --- Executes on button press in display_joint_cmica.
function display_joint_cmica_Callback(hObject, eventdata, handles)
% hObject    handle to display_joint_cmica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_cmICA_fusion_summary;

% --- Executes on button press in exit_dyn_coh.
function exit_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to exit_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));


% --- Executes on button press in run_joint_cmica.
function run_joint_cmica_Callback(hObject, eventdata, handles)
% hObject    handle to run_joint_cmica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_run_joint_cmICA;

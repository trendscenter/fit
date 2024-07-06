function varargout = neural_net_fusion(varargin)
% NEURAL_NET_FUSION MATLAB code for neural_net_fusion.fig
%      NEURAL_NET_FUSION, by itself, creates a new NEURAL_NET_FUSION or raises the existing
%      singleton*.
%
%      H = NEURAL_NET_FUSION returns the handle to a new NEURAL_NET_FUSION or the handle to
%      the existing singleton*.
%
%      NEURAL_NET_FUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURAL_NET_FUSION.M with the given input arguments.
%
%      NEURAL_NET_FUSION('Property','Value',...) creates a new NEURAL_NET_FUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before neural_net_fusion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to neural_net_fusion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help neural_net_fusion

% Last Modified by GUIDE v2.5 11-Nov-2019 17:09:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @neural_net_fusion_OpeningFcn, ...
    'gui_OutputFcn',  @neural_net_fusion_OutputFcn, ...
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


% --- Executes just before neural_net_fusion is made visible.
function neural_net_fusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to neural_net_fusion (see VARARGIN)

% Choose default command line output for neural_net_fusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes neural_net_fusion wait for user response (see UIRESUME)
% uiwait(handles.neural_net_fusion);


% --- Outputs from this function are returned to the command line.
function varargout = neural_net_fusion_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_neural_net.
function setup_neural_net_Callback(hObject, eventdata, handles)
% hObject    handle to setup_neural_net (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dynamic_coherence_setup;
ica_fuse_setup_neural_net_fusion;


% --- Executes on button press in Stats_neural_net.
function Stats_neural_net_Callback(hObject, eventdata, handles)
% hObject    handle to Stats_neural_net (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_stats_neural_net_fusion;

% --- Executes on button press in exit_dyn_coh.
function exit_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to exit_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));

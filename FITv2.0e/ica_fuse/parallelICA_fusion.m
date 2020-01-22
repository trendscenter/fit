function varargout = parallelICA_fusion(varargin)
% PARALLELICA_FUSION M-file for parallelICA_fusion.fig
%      PARALLELICA_FUSION, by itself, creates a new PARALLELICA_FUSION or raises the existing
%      singleton*.
%
%      H = PARALLELICA_FUSION returns the handle to a new PARALLELICA_FUSION or the handle to
%      the existing singleton*.
%
%      PARALLELICA_FUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARALLELICA_FUSION.M with the given input arguments.
%
%      PARALLELICA_FUSION('Property','Value',...) creates a new PARALLELICA_FUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before parallelICA_fusion_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to parallelICA_fusion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help parallelICA_fusion

% Last Modified by GUIDE v2.5 07-Mar-2017 23:53:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @parallelICA_fusion_OpeningFcn, ...
    'gui_OutputFcn',  @parallelICA_fusion_OutputFcn, ...
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


% --- Executes just before parallelICA_fusion is made visible.
function parallelICA_fusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to parallelICA_fusion (see VARARGIN)

% Choose default command line output for parallelICA_fusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

if isempty(which('ica_fuse_run_parallelICA.m'))
    addpath(genpath(fileparts(which('fusion.m'))), '-end');
end

% UIWAIT makes parallelICA_fusion wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = parallelICA_fusion_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%% Define Function Callbacks %%%%%%%%
function preprocess_snps_Callback(handleObj, event_data, handles)
% Preprocess SNPS

ica_fuse_preprocSNP;

function setup_analysis_Callback(handleObj, event_data, handles)
% Setup Analysis

ica_fuse_setup_parallelICA;

function run_analysis_Callback(handleObj, event_data, handles)
% Run Analysis

ica_fuse_run_parallelICA;

function openCallback(handleObj, event_data, handles)
% Open fMRI gene fusion file

ica_fuse_defaults;
global PARALLEL_ICA_INFO_MAT_FILE ;


fusionFile = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
    ['*', PARALLEL_ICA_INFO_MAT_FILE , '*.mat'], 'title', 'Select parallel ICA information file');

load(fusionFile);

if ~exist('paraICAInfo', 'var')
    error(['Selected file: ', fusionFile, ' is not a valid fmri-gene fusion info file']);
end

% open fusion file
ica_fuse_setup_parallelICA([], fusionFile);

function display_Callback(handleObj, event_data, handles)
% Display Callback

ica_fuse_parallelICA_displayGUI;

% --- Executes on button press in about.
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ica_fuse_titleDialog;


function help_Callback(hObject, eventdata, handles)
% Open HTML help manual

ica_fuse_openHelp;

% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));


% --------------------------------------------------------------------
function snp_estimation_Callback(hObject, eventdata, handles)
% hObject    handle to snp_estimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ica_fuse_snp_est_main;


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


% --- Executes on selection change in utilities.
function utilities_Callback(hObject, eventdata, handles)
% hObject    handle to utilities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns utilities contents as cell array
%        contents{get(hObject,'Value')} returns selected item from utilities


val = get(hObject, 'value');
val = val - 1;

if (val == 1)
    %% SNP dimensionality estimation
    ica_fuse_snp_est_main;
elseif (val == 2)
    %% Permutation test
    ica_fuse_parallelICA_permutationtest;
elseif (val == 3)
    %% Preprocess SNps
    ica_fuse_preprocSNP;
elseif (val == 4)
    %% Leave one out evaluation
    ica_fuse_leave_one_out;
end

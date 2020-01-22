function varargout = fusion(varargin)
% FUSION M-file for FUSION.fig
%      FUSION, by itself, creates a new FUSION or raises the existing
%      singleton*.
%
%      H = FUSION returns the handle to a new FUSION or the handle to
%      the existing singleton*.
%
%      FUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FUSION.M with the given input arguments.
%
%      FUSION('Property','Value',...) creates a new FUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CCA_ICA_FUSION_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fusion_openingfcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FUSION

% Last Modified by GUIDE v2.5 04-Nov-2019 13:12:36

ica_fuse_delete_gui({'fusion', 'jointICA_fusion', 'parallel_ica_toolbox', 'CCA_ICA_fusion', 'tIVA_fusion', 'mCCA_fusion', 'pgica_ica', 'neural_net_fusion'});


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @fusion_OpeningFcn, ...
    'gui_OutputFcn',  @fusion_OutputFcn, ...
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


% --- Executes just before fusion is made visible.
function fusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fusion (see VARARGIN)


%ica_fuse_delete_gui({'jointICA_fusion', 'parallel_ica_toolbox', 'CCA_ICA_fusion'});

% Choose default command line output for fusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fusion wait for user response (see UIRESUME)
% uiwait(handles.figure1);

movegui(hObject, 'center');

try
    % add this statement to fix the button color
    feature('JavaFigures', 0);
catch
end

if isempty(which('ica_fuse_run_analysis.m'))
    addpath(genpath(fileparts(which('fusion.m'))), '-end');
end


if length(varargin) == 1
    fusionType = varargin{1};
    %findobj(
    if strcmpi(fusionType, 'jointica')
        % Open joint ica
        quit_Callback(handles.quit, [], hObject);
        jointICA_fusion;
    elseif strcmpi(fusionType, 'paraica')
        % Open parallel ica
        quit_Callback(handles.quit, [], hObject);
        parallelICA_fusion;
    elseif strcmpi(fusionType, 'ccaica')
        % Open CCA + Joint ICA
        quit_Callback(handles.quit, [], hObject);
        CCA_ICA_fusion;
    elseif strcmpi(fusionType, 'mcca')
        % Open mcca
        quit_Callback(handles.quit, [], hObject);
        mCCA_fusion;
    elseif strcmpi(fusionType, 'tiva')
        % Open tiva
        quit_Callback(handles.quit, [], hObject);
        tIVA_fusion;
    elseif strcmpi(fusionType, 'pgica')
        % Open pgICA
        quit_Callback(handles.quit, [], hObject);
        pgica_ica;
    elseif strcmpi(fusionType, 'neuralnet')
        % Open neural net fusion
        quit_Callback(handles.quit, [], hObject);
        neural_net_fusion;
    elseif strcmpi(fusionType, 'exit') || strcmpi(fusionType, 'quit')
        % quit
        quit_Callback(handles.quit, [], hObject);
    elseif strcmpi(fusionType, 'help')
        ica_fuse_openHelp;
    elseif strcmpi(fusionType, 'ver')
        disp('.........................................................................');
        disp('............... Fusion ICA Toolbox version v2.0d ............................');
        disp('.........................................................................');
        disp('Fusion ICA Toolbox: v2.0d ');
        disp('Joint ICA Toolbox: v2.0a');
        disp('Parallel ICA Toolbox: v1.0b');
        disp('CCA + Joint ICA Toolbox: v1.0b');
        disp('MCCA: v1.0a');
        disp('tIVA: v1.0a');
        fprintf('\n\n');
    else
        delete(figHandle);
        error('Unknown modality specified');
    end
end




% --- Outputs from this function are returned to the command line.
function varargout = fusion_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.output;
else
    varargout{1} = [];
end


% --- Executes during object creation, after setting all properties.
function pgica_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pgica (see GCBO)
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
function parallel_ica_Callback(hObject, eventdata, handles, optional)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parallelICA_fusion;

% --- Executes on button press in pushbutton2.
function joint_ica_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

jointICA_fusion;

% --- Executes on button press in analysis_info.
function tiva_Callback(hObject, eventdata, handles)
% hObject    handle to analysis_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tIVA_fusion;


% --- Executes on button press in pushbutton4.
function mcca_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mCCA_fusion;


% % --- Executes on selection change in toolboxes.
function pgica_Callback(hObject, eventdata, handles)
% hObject    handle to toolboxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns toolboxes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from toolboxes

pgica_ica;


function fusion_help_Callback(hObject, eventdata, handles)
% Open HTML help manual

ica_fuse_openHelp;


% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));

return;


% --- Executes on button press in cca_jointica.
function cca_jointica_Callback(hObject, eventdata, handles)
% hObject    handle to cca_jointica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CCA_ICA_fusion;

% --- Executes on button press in about.
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ica_fuse_titleDialog;

% --- Executes on button press in neural_net_fusion.
function neural_net_fusion_Callback(hObject, eventdata, handles)
% hObject    handle to neural_net_fusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

neural_net_fusion;

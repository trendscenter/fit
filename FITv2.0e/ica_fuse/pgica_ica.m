function varargout = pgica_ica(varargin)
% PGICA_ICA MATLAB code for pgica_ica.fig
%      PGICA_ICA, by itself, creates a new PGICA_ICA or raises the existing
%      singleton*.
%
%      H = PGICA_ICA returns the handle to a new PGICA_ICA or the handle to
%      the existing singleton*.
%
%      PGICA_ICA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PGICA_ICA.M with the given input arguments.
%
%      PGICA_ICA('Property','Value',...) creates a new PGICA_ICA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pgica_ica_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pgica_ica_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pgica_ica

% Last Modified by GUIDE v2.5 03-Oct-2019 15:39:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @pgica_ica_OpeningFcn, ...
    'gui_OutputFcn',  @pgica_ica_OutputFcn, ...
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


% --- Executes just before pgica_ica is made visible.
function pgica_ica_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pgica_ica (see VARARGIN)

% Choose default command line output for pgica_ica
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes pgica_ica wait for user response (see UIRESUME)
% uiwait(handles.pgica_ica);


% --- Outputs from this function are returned to the command line.
function varargout = pgica_ica_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_pgica.
function setup_pgica_Callback(hObject, eventdata, handles)
% hObject    handle to setup_pgica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dynamic_coherence_setup;
ica_fuse_setup_pgica_ica;


% --- Executes on button press in display_pgica.
function display_pgica_Callback(hObject, eventdata, handles)
% hObject    handle to display_pgica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Select results file
param_file = ica_fuse_selectEntry('typeEntity', 'file', 'title', 'Select parallel group ICA + ICA fusion param file', 'filter', '*pgica*ica*mat');

drawnow;

if (isempty(param_file))
    error('PGICA-ICA parameter file is not selected');
end

load(param_file);

if (~exist('pgicaInfo', 'var'))
    error('Selected file is not a valid PGICA-ICA parameter file');
end

formatName = questdlg('Select results format', 'Results format', 'HTML', 'PDF', 'HTML');

if (~isempty(formatName))
    
    subjectLabels = cellstr([repmat('Subject ', pgicaInfo.num_subjects, 1), num2str((1:pgicaInfo.num_subjects)')]);
    
    chkGroups = questdlg('Do You Want to Specify Groups Information?', 'Groups Info', 'Yes', 'No', 'Yes');
    
    groupsInfo = [];
    if (strcmpi(chkGroups, 'yes'))
        [groupName1, groupVal1] = ica_fuse_select_groups_gui(subjectLabels, 'Group 1 Info','groups_gui');
        if (isempty(groupVal1))
            error('Figure window was quit');
        end
        groupsInfo(1).name = groupName1;
        groupsInfo(1).value = groupVal1;
        [groupName2, groupVal2] = ica_fuse_select_groups_gui(subjectLabels, 'Group 2 Info','groups_gui');
        if (isempty(groupVal2))
            error('Figure window was quit');
        end
        groupsInfo(2).name = groupName2;
        groupsInfo(2).value = groupVal2;
    end
    
    pgicaOpts.stats.groupsInfo = groupsInfo;
    
    outDir = fullfile(fileparts(param_file), [pgicaInfo.prefix, '_pgica_results']);
    if (exist(outDir, 'dir') ~= 7)
        mkdir(outDir);
    end
    opts.outputDir = outDir;
    opts.showCode = false;
    opts.useNewFigure = false;
    opts.format = lower(formatName);
    opts.createThumbnail = true;
    if (strcmpi(opts.format, 'pdf'))
        opts.useNewFigure = false;
    end
    assignin('base', 'param_file', param_file);
    assignin('base', 'pgicaOpts', pgicaOpts);
    opts.codeToEvaluate = 'ica_fuse_pgica_summary(param_file, pgicaOpts);';
    %publish('icatb_gica_html_report', 'outputDir', outDir, 'showCode', false, 'useNewFigure', false);
    disp('Generating reults summary. Please wait ....');
    drawnow;
    publish('ica_fuse_pgica_summary', opts);
    
    close all;
    
    if (strcmpi(opts.format, 'html'))
        ica_fuse_openHTMLHelpFile(fullfile(outDir, 'ica_fuse_pgica_summary.html'));
    else
        open(fullfile(outDir, 'ica_fuse_pgica_summary.pdf'));
    end
    
    disp('Done');
    
end


% --- Executes on button press in exit_dyn_coh.
function exit_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to exit_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));

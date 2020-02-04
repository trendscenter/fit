function ica_fuse_setup_pgica_ica(param_file)
%% Setup parallel groupica and ica fusion
%

chkFile = ica_fuse_questionDialog('title', 'Select Param file', 'textbody', 'Do you have an existing PGICA-ICA parameter file?');
drawnow;
if (chkFile)
    param_file =  ica_fuse_selectEntry('typeEntity', 'file', 'title', ...
        'Select PGICA-ICA fusion param file ...', 'typeSelection', 'single', 'filter', '*pgica*ica*param*');
    drawnow;
    if (isempty(param_file))
        error('PGICA-ICA param file is not selected');
    end
end

if (~exist('param_file', 'var'))
    outputDir = ica_fuse_selectEntry('typeEntity', 'directory', 'title', ...
        'Select output directory for parallel group ICA + ICA fusion');
    drawnow;
    if (isempty(outputDir))
        error('Output analysis directory is not selected');
    end
    
    
    dlg_title = 'Enter modality name (first modality)';
    
    numParameters = 1;
    
    inputText(numParameters).promptString = 'Enter first modality name';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'sMRI', 'Gene', 'EEG', 'Behavioral'};
    inputText(numParameters).answerType = 'string';
    inputText(numParameters).tag = 'modality1_name';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    answer = ica_fuse_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');
    
    if (isempty(answer))
        error('Modality name is not selected');
    end
    
    modality1_name = lower(answer{1});
    
    
    filterPattern = '*.asc;*.dat;*.txt';
    if (strcmpi(modality1_name, 'smri'))
        filterPattern = '*.nii;*.img';
    end
    
    prefix = '';
    modality1_files =  ica_fuse_selectEntry('typeEntity', 'file', 'title', ...
        'Select First Modality files for all subjects ...', 'typeSelection', 'multiple', 'filter', filterPattern);
    drawnow;
    if (isempty(modality1_files))
        error('First modality files are not selected');
    end
    pgicaInfo.files{1} = modality1_files;
    
    modality2_files =  ica_fuse_selectEntry('typeEntity', 'file', 'title', ...
        'Select fMRI files (second modality) for all subjects ...', 'typeSelection', 'multiple', 'filter', '*.nii');
    drawnow;
    if (isempty(modality2_files))
        error('Second modality files are not selected');
    end
    pgicaInfo.files{2} = modality2_files;
    
    pgicaInfo.modalities = {modality1_name, 'fmri'};
    
else
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    load(param_file);
    
    if (~exist('pgicaInfo', 'var'))
        error('Selected file is not a valid PGICA parameter file');
    end
    
    prefix = pgicaInfo.prefix;
    
end

cd(outputDir);

numComp1 = '';
if (isfield(pgicaInfo, 'numComp1'))
    numComp1 = num2str(pgicaInfo.numComp1);
end

numComp2 = '';
if (isfield(pgicaInfo, 'numComp2'))
    numComp2 = num2str(pgicaInfo.numComp2);
end


pgicaInfo.outputDir = outputDir;


numParameters = 1;

% Prefix
inputText(numParameters).promptString =  'Enter analysis files output prefix';
inputText(numParameters).answerString =  prefix;
inputText(numParameters).tag = 'prefix';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'string';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;


numParameters = numParameters + 1;

mval = 1;
if (isfield(pgicaInfo, 'maskFile_modality1') && ~isempty(pgicaInfo.maskFile_modality1))
    mval = 2;
end

% Mask Selection
inputText(numParameters).promptString =  'What Mask Do You Want To Use For First Modality?';
inputText(numParameters).answerString =  char('Default Mask', 'Select Mask');
inputText(numParameters).tag = 'maskFile_modality1';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = mval;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;

mval = 1;
if (isfield(pgicaInfo, 'maskFile_modality2') && ~isempty(pgicaInfo.maskFile_modality2))
    mval = 2;
end

numParameters = numParameters + 1;

% Mask Selection
inputText(numParameters).promptString =  'What Mask Do You Want To Use For Second Modality?';
inputText(numParameters).answerString =  char('Default Mask', 'Select Mask');
inputText(numParameters).tag = 'maskFile_modality2';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = mval;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 0;


numParameters = numParameters + 1;

% Number of Independent Components
inputText(numParameters).promptString =  'Number of PC for modality 1';
inputText(numParameters).answerString =  numComp1;
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'numComp1';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;


numParameters = numParameters + 1;

% Number of Independent Components
inputText(numParameters).promptString =  'Number of PC for modality 2. Enter first level PCA followed by second level like [10, 8]';
inputText(numParameters).answerString = numComp2;
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'numComp2';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).read_from_file = 0;
inputText(numParameters).executeCallback = 1;


InputHandle = ica_fuse_plot_controls_fig(inputText, 'PGICA_params', 'on', 'Apply', 'Cancel');

set(InputHandle, 'userdata', pgicaInfo);

set(findobj(InputHandle, 'tag', 'maskFile_modality1'), 'callback', {@maskSelection, InputHandle});

set(findobj(InputHandle, 'tag', 'maskFile_modality2'), 'callback', {@maskSelection, InputHandle});

set(findobj(InputHandle, 'tag', 'Apply'), 'callback',  {@applyCallback, InputHandle});

try
    waitfor(InputHandle);
catch
end


function maskSelection(hObject, event_data, handles)
%% Mask selection

pgicaInfo = get(handles, 'userdata');
tagVal = get(hObject, 'tag');
val = get(hObject, 'value');
str=get(hObject, 'string');

if (~strcmpi(str(val, :), 'default mask'))
    maskFile = ica_fuse_selectEntry('typeEntity', 'file', 'title', 'Select binary mask file ...', 'typeSelection', 'single', 'filter', '*.nii;*.asc;*.txt');
    
    drawnow;
    if (isempty(maskFile))
        set(hObject, 'value', 1);
    end
    pgicaInfo.(tagVal) = maskFile;
else
    pgicaInfo.(tagVal)='';
end

set(handles, 'userdata', pgicaInfo);



function applyCallback(hObject, event_data, handles)
%% Mask selection

pgicaInfo = get(handles, 'userdata');

%% Get params
prefix = get(findobj(handles, 'tag', 'prefix'), 'string');
pgicaInfo.param_file = [prefix, '_pgica_ica_param_info.mat'];
pgicaInfo.prefix = prefix;

if (~isfield(pgicaInfo, 'maskFile_modality1'))
    pgicaInfo.maskFile_modality1 = '';
end

if (~isfield(pgicaInfo, 'maskFile_modality2'))
    pgicaInfo.maskFile_modality2 = '';
end

numComp1 = str2num(get(findobj(handles, 'tag', 'numComp1'), 'string'));
numComp2 = str2num(get(findobj(handles, 'tag', 'numComp2'), 'string'));

if (isempty(numComp1))
    error('Number of independent components is not selected for modality 1');
end


if (isempty(numComp2))
    error('Number of independent components is not selected for modality 2');
end

numComp1 = numComp1(1);
if (length(numComp2) == 1)
    numComp2 = [numComp2, numComp2];
end

if (length(numComp2) > 2)
    numComp2 = numComp2(1:2);
end


pgicaInfo.numComp1 = numComp1;
pgicaInfo.numComp2 = numComp2;

ica_options = ica_fuse_paraICAOptions(min([numComp1, numComp2]), 'on');

drawnow;

if (isempty(ica_options))
    error('Options window was quit');
end

pgicaInfo.ica_options = ica_options;

fileName = fullfile(pgicaInfo.outputDir, pgicaInfo.param_file);
disp(['Saving PGICA fusion parameters in file ', fileName]);
save(fileName, 'pgicaInfo');
disp('Done');
fprintf('\n');


delete (handles);

drawnow;

ica_fuse_run_pgica_ica(pgicaInfo);

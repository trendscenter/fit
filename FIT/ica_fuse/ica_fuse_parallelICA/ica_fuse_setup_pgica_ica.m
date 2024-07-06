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
    answer = ica_fuse_inputdlg2({'Enter number of modalities'}, 'No. of Modalities', 1, {'2'});
    if (isempty(answer) || isempty(answer{1}))
        error ('No of modalities not entered');
    end
    num_features = str2num(answer{1});
    if ((num_features ~= 2) && (num_features ~= 3))
        error('Number of modalities should be either 2 or 3');
    end
else
    
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    load(param_file);
    
    if (~exist('pgicaInfo', 'var'))
        error('Selected file is not a valid PGICA parameter file');
    end
    
    num_features = length(pgicaInfo.modalities);
    
    if (~isfield(pgicaInfo, 'featuresInfo'))
        featuresInfo = repmat(struct('name', '', 'modality_type', '', 'components', '', 'files', '', 'mask', ''), 1, num_features);
        
        for nF = 1:num_features
            featuresInfo(nF).name = ['Feature ', num2str(nF)];
            featuresInfo(nF).modality_type = pgicaInfo.modalities{nF};
            featuresInfo(nF).components = pgicaInfo.(['numComp', num2str(nF)]);
            featuresInfo(nF).mask = pgicaInfo.(['maskFile_modality', num2str(nF)]);
            featuresInfo(nF).files = pgicaInfo.files{nF};
        end
        pgicaInfo.featuresInfo = featuresInfo;
    end
end

pgicaInfo.num_features = num_features;
pgicaInfo.outputDir = outputDir;

%% Select options
pgicaInfo = pgica_select_options(pgicaInfo);

%% Run Analysis
ica_fuse_run_pgica_ica(pgicaInfo);


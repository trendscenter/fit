function ica_fuse_run_analysis(fusionFile, resume_analysis)
%% Run fusion analysis
%
% Inputs:
% 1. fusionFile - full file path of the fusion file.
% 2. resume_analysis - Resume analysis. Options are 0 and 1. A value of 0
% will restart the analysis.
%
% Analysis steps:
% 1. All steps
% 2. PCA
% 3. ICA
% 4. Back reconstruct
% 5. Scale components (Only for the first combination)
%

try
    
    %% Load defaults
    ica_fuse_defaults;
    
    global FUSION_INFO_MAT_FILE; % Fusion information file
    global JOINT_COMPONENT_NAMING; % JOINT COMPONENT IMAGE FILE NAMING
    global FLIP_ANALYZE_IM; % Flip parameter for analyze images
    global OPTIMIZE_FEATURES;
    global STACK_ALL_FEATURES;
    global OPEN_DISPLAY_WINDOW;
    
    %% Open Figure
    isGUI = 0;
    if (~exist('fusionFile', 'var'))
        isGUI = 1;
        fusionFile = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
            ['*', FUSION_INFO_MAT_FILE, '*.mat'], 'title', 'Select fusion information file for running fusion analysis');
        drawnow;
    end
    
    load(fusionFile);
    
    if (~exist('fusionInfo', 'var'))
        error('Error:Fusion', 'Selected file %s is not a valid fusion file\n', fusionFile);
    end
    
    %% Check setup analysis
    if (~isfield(fusionInfo, 'setup_analysis'))
        error('Error:SetupAnalysis', 'Please setup analysis using Setup ICA\n');
    end
    
    %% Get all possible combinations
    if isfield(fusionInfo.setup_analysis, 'type_pca') && (strcmpi(fusionInfo.setup_analysis.type_pca, 'cca') ||  strcmpi(fusionInfo.setup_analysis.type_pca, 'mcca') || strcmpi(fusionInfo.setup_analysis.type_pca, 'mccar'))
        disp('Optimal features utility is turned off when using CCA or MCCA data reduction strategy');
        disp('');
        OPTIMIZE_FEATURES = 'no';
    end
    
    ica_algo = ica_fuse_icaAlgorithm;
    ica_algo = cellstr(ica_algo);
    algorithmName = ica_algo{fusionInfo.setup_analysis.algorithm};
    
    if (strcmpi(algorithmName, 'iva-g') || strcmpi(algorithmName, 'iva-ggd'))
        OPTIMIZE_FEATURES = 'no';
    end
    
    optimalFeatures = strcmpi(OPTIMIZE_FEATURES, 'yes');
    numFeatures = fusionInfo.setup_analysis.numFeatures;
    all_comb =  ica_fuse_get_combinations(numFeatures, STACK_ALL_FEATURES, optimalFeatures);
    
    if (~isfield(fusionInfo.run_analysis, 'all_comb'))
        fusionInfo.run_analysis.all_comb = all_comb;
    end
    
    %% Check if there are any combination/combinations left to be run
    [combnsToRun, stepsToRun] = ica_fuse_check_analysis(fusionInfo);
    
    if (~isGUI)
        
        %% Start fresh analysis if not specified.
        if (~exist('resume_analysis', 'var'))
            resume_analysis = 0;
        end
        
        if (resume_analysis == 1) && (isempty(combnsToRun))
            disp(sprintf(['No combinations and analysis steps left to be run.\nPlease use GUI or command ica_fuse_run_analysis(fusionFile) ', ...
                'to run a new analysis.\n']));
            return;
        end
        
        if (resume_analysis == 2)
            return;
        end
        
    else
        
        resume_analysis = 0;
        
        %% Open question dialog if needed
        if (length( fusionInfo.run_analysis.all_comb) == length(combnsToRun)) && (length(stepsToRun) == 4)
            resume_analysis = 0;
        else
            if (~isempty(combnsToRun))
                % Get the analysis type
                resume_analysis = getAnalysisType;
                if (resume_analysis == 2)
                    return;
                end
            end
        end
        
    end
    
    drawnow;
    
    %% Output directory
    outputDir = fileparts(fusionFile);
    if isempty(outputDir)
        outputDir = pwd;
    end
    
    cd(outputDir);
    
    %% Remove the field run_analysis from fusionInfo variable for fresh
    % analysis
    if (~resume_analysis)
        if (isGUI)
            stepsToRun = selectAnalysisSteps;
        else
            stepsToRun = 1;
        end
        stepsToRun = sort(stepsToRun);
        if (any(stepsToRun == 1) || (length(intersect(stepsToRun, (2:5))) == 4))
            stepsToRun = 1;
            fusionInfo = rmfield(fusionInfo, 'run_analysis');
        end
        
        fusionInfo.run_analysis.all_comb = all_comb;
        %% Run all combination/combinations and all steps for fresh
        % analysis
        combnsToRun = find(ica_fuse_good_cells(fusionInfo.run_analysis.all_comb) ~= 0);
        %% Copy required fields from fusionInfo.setup_analysis to
        % fusionInfo.run_analysis
        fusionInfo = ica_fuse_copyfields_to_run(fusionInfo);
    end
    
    if isfield(fusionInfo.setup_analysis, 'type_pca') && (strcmpi(fusionInfo.setup_analysis.type_pca, 'cca') || strcmpi(fusionInfo.setup_analysis.type_pca, 'mcca'))
        fusionInfo.run_analysis.cca_opts = fusionInfo.setup_analysis.cca_opts;
    end
    
    try
        fusionInfo.run_analysis.numICARuns = fusionInfo.setup_analysis.numICARuns;
    catch
    end
    
    try
        fusionInfo.run_analysis.type_ica = fusionInfo.setup_analysis.type_ica;
    catch
    end
    
    %% Output directory
    fusionInfo.run_analysis.outputDir = outputDir;
    
    %% Set initialize field
    fusionInfo.run_analysis.isInitialized = 1;
    
    %% Store fusion file
    fusionInfo.run_analysis.fusionFile = fusionFile;
    
    %% Flip parameter for analyze images
    fusionInfo.run_analysis.flip_analyze_images = FLIP_ANALYZE_IM;
    
    %% Handle interpolated EEG data
    if (~isfield(fusionInfo.setup_analysis, 'newDims'))
        fusionInfo.setup_analysis.newDims = ica_fuse_compute_new_dims(fusionInfo.run_analysis.dims, ...
            cellstr(str2mat(fusionInfo.setup_analysis.dataInfo(1).feature.modality)), fusionInfo.run_analysis.mask_ind);
    end
    fusionInfo.run_analysis.newDims = fusionInfo.setup_analysis.newDims;
    
    
    all_comb = fusionInfo.run_analysis.all_comb;
    fusionPrefix = fusionInfo.run_analysis.prefix;
    
    %% Record analysis information in a log file
    diaryFile = fullfile(outputDir, [fusionPrefix, '_results.log']);
    diary(diaryFile);
    
    tic;
    
    fprintf('\n');
    if (~resume_analysis)
        disp('Starting Analysis ...');
    else
        disp('Resuming Analysis ...');
    end
    
    fprintf('\n');
    
    %% Initialise normalization and standard deviation parameters
    featureDataLength = repmat(struct('Length', []), 1, length(all_comb));
    featureNormPara = zeros(1, numFeatures);
    stdParameters = repmat(struct('stdValues', []), 1, numFeatures);
    
    if (resume_analysis == 1)
        try
            featureNormPara = fusionInfo.run_analysis.featureNormPara;
            stdParameters = fusionInfo.run_analysis.stdParameters;
            featureDataLength = fusionInfo.run_analysis.featureDataLength;
        catch
        end
    end
    
    
    %% LOOP OVER NUMBER OF COMBINATIONS
    for nComb = combnsToRun
        
        % Current combination number
        combNumber = all_comb{nComb};
        
        if isempty(combNumber)
            continue;
        end
        
        %% Prepare data for joint ICA fusion analysis
        if (nComb == 1)
            % Return mean data for each feature and mean data for each
            % feature and group
            stackInfo = ica_fuse_prepare_data('dataInfo', fusionInfo.run_analysis.dataInfo, 'mask_ind', fusionInfo.run_analysis.mask_ind, ...
                'normalize_scheme', fusionInfo.run_analysis.normalize, 'voxels', fusionInfo.run_analysis.voxels, 'sel_groups', ...
                (1:fusionInfo.run_analysis.numGroups), 'sel_features', combNumber, 'compute_mean', 1, 'num_subjects', ...
                fusionInfo.run_analysis.numSubjects);
        else
            if (stepsToRun(1) == 5)
                break;
            end
            stackInfo = ica_fuse_prepare_data('dataInfo', fusionInfo.run_analysis.dataInfo, 'mask_ind', fusionInfo.run_analysis.mask_ind, ...
                'normalize_scheme', fusionInfo.run_analysis.normalize, 'voxels', fusionInfo.run_analysis.voxels, 'sel_groups', ...
                (1:fusionInfo.run_analysis.numGroups), 'sel_features', combNumber);
        end
        
        % Store the normalization and standard deviation parameters
        featureNormPara(combNumber) = stackInfo.featureNormPara;
        stdParameters(combNumber) = stackInfo.stdParameters;
        
        % Get information from stackInfo
        combinationName = stackInfo.combName;
        meanData = stackInfo.meanData; stackInfo.meanData = [];
        meanDataGroups = stackInfo.meanDataGroups; stackInfo.meanDataGroups = [];
        timeAxis = stackInfo.timeAxis; stackInfo.timeAxis = [];
        fileInfo = stackInfo.fileInfo; stackInfo.fileInfo = [];
        
        featureDataLength(nComb).Length = stackInfo.dataLength;
        
        %% Store the information in increments
        fusionInfo.run_analysis.featureNormPara = featureNormPara;
        fusionInfo.run_analysis.stdParameters = stdParameters;
        fusionInfo.run_analysis.featureDataLength = featureDataLength;
        
        %% Add current combination number and name
        fusionInfo.run_analysis.currentComb = nComb;
        fusionInfo.run_analysis.currentCombName = combinationName;
        
        %% Loop over steps to run
        for sT = stepsToRun
            
            if (sT == 1 || sT == 2)
                %% Data reduction step
                fusionInfo = ica_fuse_dataReduction(fusionInfo, stackInfo.data);
            end
            
            if (sT == 1 || sT == 3)
                %% Calculate ICA
                fusionInfo = ica_fuse_calculateICA(fusionInfo);
            end
            
            if (sT == 1 || sT == 4)
                %% Back reconstruction
                fusionInfo = ica_fuse_backReconstruct(fusionInfo, stackInfo.data);
            end
            
            %% Write images for combination one only
            if (nComb == 1)
                if (sT == 1 || sT == 5)
                    
                    %% Component calibration
                    [icasig, A, groups_icasig, fusionInfo] = ica_fuse_scaleComponents(fusionInfo, meanData, meanDataGroups, ...
                        featureDataLength(nComb).Length);
                    
                    outFile = [fusionInfo.run_analysis.prefix, JOINT_COMPONENT_NAMING];
                    
                    disp(['Writing data for component set ', outFile, ' ...']);
                    
                    fprintf('\n');
                    
                    % Write data as images or ascii files
                    outputFiles  =  ica_fuse_write_feature_data(icasig, A, groups_icasig, fileInfo, fusionInfo.run_analysis.mask_ind, ...
                        outFile, outputDir, timeAxis, fusionInfo.run_analysis.numSubjects);
                    fusionInfo.run_analysis.outputFiles = outputFiles;
                    ica_fuse_save(fusionFile, 'fusionInfo');
                    
                    fprintf('\n');
                    
                    disp(['Done writing data for component set ', outFile]);
                    
                    fprintf('\n');
                    clear groups_icasig outputFiles;
                    
                    clear A icasig;
                end
                %% End for writing images for combination one only
            end
            
        end
        %% End loop over steps to run
        
        if (resume_analysis)
            %% Reset steps to run for the next combination for resume step
            stepsToRun = 1;
        end
        
        clear stackInfo;
        
    end
    %% END LOOP OVER NUMBER OF COMBINATIONS
    
    
    %% All combination names
    featureNames = cellstr(str2mat(fusionInfo.run_analysis.dataInfo(1).feature.name));
    allCombNames = repmat({''}, length(all_comb), 1);
    good_inds = find(ica_fuse_good_cells(fusionInfo.run_analysis.all_comb) ~= 0);
    
    for nComb = good_inds
        currentCombNo = all_comb{nComb};
        allCombNames{nComb} = ica_fuse_formatStr(featureNames(currentCombNo), ' & ');
        %% Handle missing information
        if (isempty(featureDataLength(nComb).Length))
            featureDataLength(nComb).Length = fusionInfo.run_analysis.newDims(currentCombNo);
        end
    end
    
    fusionInfo.run_analysis.allCombNames = allCombNames;
    fusionInfo.run_analysis.featureDataLength = featureDataLength;
    
    clear featureDataLength;
    
    fprintf('\n');
    %% Save the information for running the analysis
    ica_fuse_save(fusionFile, 'fusionInfo');
    disp(['Fusion information is saved in file: ', fusionFile]);
    
    t_end = toc;
    
    disp(['Time taken to run the analysis is ', num2str(t_end), ' seconds']);
    fprintf('\n');
    disp(['All the analysis information is stored in log file: ', diaryFile]);
    fprintf('\n');
    
    diary('off');
    
    if (OPEN_DISPLAY_WINDOW && ~isempty(all_comb{1}))
        ica_fuse_display(fusionFile);
    end
    
catch
    
    diary('off');
    %rethrow(lasterror);
    ica_fuse_displayErrorMsg;
end

function analysisType = getAnalysisType

%% Load defaults
ica_fuse_defaults;

global FIG_BG_COLOR; % Figure background color
global FIG_FG_COLOR; % Figure foreground color
global BUTTON_BG_COLOR; % background color for pushbutton
global BUTTON_FG_COLOR; % Foreground color for pushbutton

DefaultFigureColor = get(0, 'DefaultFigureColor');
DefaulTextColor = get(0, 'DefaultTextColor');
DefaultUIBackgroundcolor = get(0, 'DefaultUICONTROLBackgroundcolor');
DefaultUIForegroundcolor = get(0, 'DefaultUICONTROLForegroundcolor');

answerQuestion = '';

try
    
    %% Set colors
    set(0, 'DefaultFigureColor', FIG_BG_COLOR);
    set(0, 'DefaultTextColor', FIG_FG_COLOR);
    set(0, 'DefaultUICONTROLBackgroundcolor', BUTTON_BG_COLOR);
    set(0, 'DefaultUICONTROLFOREGROUNDCOLOR', BUTTON_FG_COLOR);
    
    answerQuestion = questdlg('How do you want to run the analysis?', 'Analysis Type', 'Restart', 'Resume', 'Cancel', ...
        'Resume');
    
catch
    
end

drawnow;

if strcmpi(answerQuestion, 'resume')
    analysisType = 1;
elseif strcmpi(answerQuestion, 'restart')
    analysisType = 0;
else
    analysisType = 2;
end

%% Revert to the previous colors
set(0, 'DefaultFigureColor', DefaultFigureColor);
set(0, 'DefaultTextColor', DefaulTextColor);
set(0, 'DefaultUICONTROLBackgroundcolor', DefaultUIBackgroundcolor);
set(0, 'DefaultUICONTROLFOREGROUNDCOLOR', DefaultUIForegroundcolor);


function stepsToRun = selectAnalysisSteps
%% Select analysis step/steps to run
%

stepsToRun = ica_fuse_listdlg('promptstring', 'Which Step/Steps Of The Analysis You Want To Run', 'liststring', ...
    {'All', 'Data Reduction', 'Calculate ICA', 'Back Reconstruction', 'Scaling Components'}, 'title_fig', 'Run Analysis', 'help', struct('title', 'Run Analysis', 'str', char('You could select step/steps of the analysis to run. These need to be in serial order. For example if you run ICA, subsequent steps need to be run inorder to see the changes. Steps are as follows:', ...
    '', '1. All - All the analysis steps are run. This is same as running from data reduction to scaling components.', ...
    '2. Data Reduction - Data is organized as subjects by samples and reduced in the subject dimension to few components.', ...
    '3. Calculate ICA - ICA is run on the reduced data and maximally independent components are found.', ...
    '4. Back reconstruction - Mixing coefficients of the aggregate components are computed. Also computed are the backreconstructed components of each group.', ...
    '5. Scaling components - Components are scaled to data units or z-scores.')));

if (isempty(stepsToRun))
    error('Analysis step/steps are not selected');
end

if (any(stepsToRun == 1))
    stepsToRun = 1;
end

drawnow;

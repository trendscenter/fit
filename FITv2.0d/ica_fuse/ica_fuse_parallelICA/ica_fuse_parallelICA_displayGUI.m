function ica_fuse_parallelICA_displayGUI(parallel_ica_fusion_file)
% Parallel ICA display GUI

ica_fuse_defaults;
global PARALLEL_ICA_INFO_MAT_FILE;
global UI_FONT_NAME;
global UI_FONT_SIZE;
global ANATOMICAL_FILE;
global HELP_FG_COLOR;

% Display Defaults
global IMAGE_VALUES;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGES_PER_FIGURE;
global ANATOMICAL_PLANE;

if ~exist('parallel_ica_fusion_file', 'var')
    parallel_ica_fusion_file = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        ['*', PARALLEL_ICA_INFO_MAT_FILE, '*.mat'], 'title', 'Select parallel ICA fusion information file');
end

if isempty(parallel_ica_fusion_file)
    error(['Parallel ICA fusion information file is not selected']);
end

outputDir = fileparts(parallel_ica_fusion_file);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

drawnow;

load(parallel_ica_fusion_file);

if ~exist('paraICAInfo', 'var')
    error(['Selected file: ', parallel_ica_fusion_file, ' is not a valid parallel ICA fusion information file']);
end

if isfield(paraICAInfo.run_analysis, 'outputFiles')
    % Output files
    outputFiles = paraICAInfo.run_analysis.outputFiles;
else
    error('outputFiles field doesn''t exist. Please run the analysis');
end

% Output prefix
output_prefix = paraICAInfo.run_analysis.prefix;
% dataInfo
dataInfo = paraICAInfo.run_analysis.dataInfo;
% GroupNames
groupNames = str2mat(paraICAInfo.run_analysis.dataInfo.name);
% Feature names
featureNames = str2mat(paraICAInfo.run_analysis.dataInfo(1).feature.name);
% Number of groups
numGroups = paraICAInfo.run_analysis.numGroups;
% Number of features
numFeatures = paraICAInfo.run_analysis.numFeatures;
% voxel indices
mask_ind = paraICAInfo.run_analysis.mask_ind;
% Number of subjects
numSubjects = paraICAInfo.run_analysis.numSubjects;

% Number of components
numComp = paraICAInfo.run_analysis.numComp;


% Average correlation
avg_correlation = paraICAInfo.run_analysis.average_correlation;

corrIndices = paraICAInfo.run_analysis.corrIndices;

% Output files
outputFiles = paraICAInfo.run_analysis.outputFiles;

% Flip parameter for analyze images
flip_analyze_images = [];
if isfield(paraICAInfo.run_analysis, 'flip_analyze_images')
    flip_analyze_images = paraICAInfo.run_analysis.flip_analyze_images;
end


% Construct display parameters structure
displayParameters = struct('outputDir', outputDir, 'output_prefix', output_prefix, 'dataInfo', dataInfo, 'groupNames', groupNames, ...
    'featureNames', featureNames, 'numGroups', numGroups, 'numFeatures', numFeatures, 'mask_ind', mask_ind, ...
    'numSubjects', numSubjects, 'numComp', numComp, 'avg_correlation', avg_correlation, 'corrIndices', corrIndices, ...
    'outputFiles', outputFiles, 'parallel_ica_fusion_file', parallel_ica_fusion_file, 'selectedFeatureVal', 1, 'flip_analyze_images', ...
    flip_analyze_images);

% Draw a figure with two listboxes, a check box, Convert to z scores,
% Threshold, Slices, Anatomical Plane
displayTag = 'DisplayGUI_para_ica_fusion';

%%%%%% Delete any previous figures of display GUI %%%%%
displayH = findobj('tag', displayTag);

for ii = 1:length(displayH)
    delete(displayH(ii));
end
%%%%%% end for deleting previous figures of display GUI %%%%

% Display GUI figure
[displayHandle] = ica_fuse_getGraphics('Display GUI For Parallel ICA', 'displaygui', displayTag, 'off');
set(displayHandle, 'menubar', 'none');

% Plot Title for the figure
% parameters to plot in a menu

% offsets
xOffset = 0.04; yOffset = 0.03;


%%%%%%%%%%%%% Draw Title here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title color
titleColor = [0 0.9 0.9];
% fonts
titleFont = 13;
titleAxesH = axes('Parent', displayHandle, 'position', [0 0 1 1], 'visible', 'off');
axis(titleAxesH, 'off');
xPos = 0.5; yPos = 0.97;
text(xPos, yPos, 'Display GUI For Parallel ICA', 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', ...
    'fontsize', titleFont, 'HorizontalAlignment', 'center', 'FontName', UI_FONT_NAME);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot display button
buttonWidth = 0.2; buttonHeight = 0.05;
displayButtonPos = [0.75 yOffset buttonWidth buttonHeight];


displayButtonH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', ...
    'pushbutton', 'position', displayButtonPos, 'string', 'Display', ...
    'tag', 'display_button', 'callback', {@displayCallback, displayHandle});

%%%%%%%%%%%%%% Component Listbox %%%%%%%%%%%%%%%%%%%%

compListWidth = 0.20; listboxWidth = 0.4; yPos = 0.92;

% Plot textbox
textBoxWidth = compListWidth; textboxHeight = 0.05;
textboxPos = [xOffset, yPos - yOffset, textBoxWidth, textboxHeight];

% Plot component text
compTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Comp No.');


featureListbox = repmat(struct('string', []), 1, length(numComp));
for nF = 1:length(numComp)
    % Component string structure
    compStr = repmat(struct('name', ''), 1, numComp(nF));
    
    % form component string
    for nComp = 1:numComp(nF)
        compStr(nComp).name = num2str(nComp);
    end
    featureListbox(nF).string = str2mat(compStr.name);
    clear compStr;
end

displayParameters.featureListbox = featureListbox;
clear featureListbox;

listboxPos = textboxPos;
listboxPos(4) = 0.22;
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);
% Plot component listbox
compListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', displayParameters.featureListbox(1).string, 'value', [1:size(displayParameters.featureListbox(1).string, 1)], ...
    'min', 0, 'max', 2,  'tag', 'selComp');
%%%%%%%%%%%%% End for plotting component listbox and text %%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% Features Text and Listbox %%%%%%%%%%%%%%%%%%%%%%
textboxPos(1) = 1 - xOffset - listboxWidth;
textboxPos(3) = listboxWidth;

% Plot Feature text
featureTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', 'Feature');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
% Plot feature listbox
featureListH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', displayParameters.featureNames, 'value', 1, 'min', 0, 'max', 1, 'tag', 'selFeature', ...
    'callback', {@featureListCallback, displayHandle});

% Textbox width
clear textboxPos;
%%%%%%%%%%%%%%%%%%%%%%% End for plotting features listbox and text %%%%%%%

% Plot drop down box for sorting joint components
textPos = [xOffset, listboxPos(2) - yOffset, 0.3, 0.05];
textToPlot = {'How do you want to sort components?'};

sortTextH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textPos, 'String', textToPlot);
[textToPlot, newPos] = textwrap(sortTextH, textToPlot);
textPos(4) = newPos(4);
textPos(2) = listboxPos(2) - yOffset - textPos(4);

set(sortTextH, 'string', textToPlot);
set(sortTextH, 'position', textPos);

popupPos = textPos;
popupWidth = 0.45;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = popupWidth;
popupPos(4) = 0.05;
popupPos(2) = textPos(2) + 0.5*textPos(4) - 0.5*popupPos(4);
sortPopupH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupPos, 'String', char('None', 'Correlation between modalities', 'Two sample t-test on mixing coefficients'), 'tag', 'sort_components', 'callback', ...
    {@sortPopupCallback, displayHandle});

helpButtonPos = popupPos;
helpButtonPos(1) = helpButtonPos(1) + helpButtonPos(3) + xOffset;
helpButtonPos(3) = 0.06;

helpButtonH = ica_fuse_uicontrol('parent', displayHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', helpButtonPos, 'String', '?', 'tag', 'help_sort_components', 'callback', ...
    {@helpCallback, displayHandle}, 'foregroundcolor', HELP_FG_COLOR);

% plot display defaults
display_defaultsMenu = uimenu('parent', displayHandle, 'label', 'Display Defaults', 'callback', ...
    {@display_defaults_callback, displayHandle});

[controlPara, displayParameters] = ica_fuse_get_display_parameters(displayParameters, [], 'off');

set(compListH, 'userdata', controlPara);

set(displayHandle, 'userdata', displayParameters);

% Plot options menu
optionsMenuH = uimenu('parent', displayHandle, 'label', 'Options');
resultsSummaryH = uimenu('parent', optionsMenuH, 'label', 'Results Summary', 'callback', ...
    {@resultsSummaryCallback, displayHandle});
selectLocusFileH = uimenu('parent', optionsMenuH, 'label', 'Select Locus File', 'callback', ...
    {@selectLocusFile, displayHandle, 1}, 'tag', 'locusFile');

% help on display GUI
fitHelpTitle = uimenu('parent', displayHandle, 'label', 'FIT-Help');
fitHelpMenu = uimenu(fitHelpTitle, 'label', 'Display GUI For Parallel ICA', 'callback', ...
    'ica_fuse_openHTMLHelpFile(''fit_paraICA_display.htm'');');

anatomicalPlaneH = findobj(displayHandle, 'tag', 'anatomical_plane');
set(anatomicalPlaneH , 'callback', {@anatomicalPopupCallback, displayHandle});

% Make the graphics visible after plotting all controls
set(displayHandle, 'visible', 'on');

sortPopupCallback(findobj(displayHandle, 'tag', 'sort_components'), [], displayHandle);



function displayCallback(hObject, event_data, handles)
% hObject - Display push button
% handles - Display GUI figure
% Purpose: Display the components by features

set(handles, 'pointer', 'watch');
ica_fuse_defaults;
global ANATOMICAL_FILE;
global SNP_Z_THRESHOLD;

try
    
    % Component listbox
    compListH = findobj(handles, 'tag', 'selComp');
    compString = get(compListH, 'string');
    
    selectLocusFileH = findobj(handles, 'tag', 'locusFile');
    
    displayParameters = get(handles, 'userdata');
    
    % Parallel ICA fusion file
    parallel_ica_fusion_file = displayParameters.parallel_ica_fusion_file;
    
    % Number of groups and features
    numGroups = displayParameters.numGroups;
    numFeatures = displayParameters.numFeatures;
    % Number of components
    numComp = displayParameters.numComp;
    % Number of subjects
    numSubjects = displayParameters.numSubjects;
    
    flip_analyze_images = displayParameters.flip_analyze_images;
    
    % Output directory
    outputDir = displayParameters.outputDir;
    
    cd(outputDir);
    
    % Output prefix
    output_prefix = displayParameters.output_prefix;
    
    % data info
    dataInfo = displayParameters.dataInfo;
    % Output files
    outputFiles = displayParameters.outputFiles;
    
    % Voxels
    mask_ind = displayParameters.mask_ind;
    
    % corr indices
    corrIndices = displayParameters.corrIndices;
    
    % Modalities
    modalities = cellstr(str2mat(dataInfo(1).feature.modality));
    
    % Feature names
    featureNames = cellstr(str2mat(dataInfo(1).feature.name));
    
    groupNames = char(dataInfo.name);
    
    [dims, voxels] = ica_fuse_getFeatureDIM(mask_ind);
    
    load(displayParameters.parallel_ica_fusion_file);
    
    locusFileName = '';
    if isfield(paraICAInfo.run_analysis, 'locusFileName')
        locusFileName = paraICAInfo.run_analysis.locusFileName;
    end
    
    type_parallel_ica = paraICAInfo.run_analysis.type_parallel_ica;
    
    % Get the display fields
    % Selected components
    [answerString, selComp] = ica_fuse_get_value_uicontrol(handles, 'selComp');
    displayParameters.selectedComp = selComp;
    
    % Selected Features
    [selectedFeature, selectedFeatureVal] = ica_fuse_get_value_uicontrol(handles, 'selFeature');
    selectedFeature = deblank(selectedFeature(selectedFeatureVal, :));
    displayParameters.selectedFeature = selectedFeature;
    displayParameters.selectedFeatureVal = selectedFeatureVal;
    
    % Retrieve display parameters
    convertToZ = displayParameters.convert_to_z;
    imageValues = displayParameters.image_values;
    images_per_figure = displayParameters.images_per_figure;
    threshold = abs(displayParameters.z_threshold);
    anatomical_plane = lower(displayParameters.anatomical_plane);
    slices_in_mm = displayParameters.slices_in_mm;
    
    if (min(threshold) == max(threshold))
        threshold = threshold(1);
    end
    
    
    % Sort components handle
    sortCompH = findobj(handles, 'tag', 'sort_components');
    sortCompVal = get(sortCompH, 'value');
    sortCompStr = get(sortCompH, 'string');
    
    %sort_components = 0;
    sortComponentsStr = deblank(sortCompStr(sortCompVal, :));
    %     if ~strcmpi(sortComponentsStr, 'none')
    %         sort_components = 1;
    %     end
    
    selectedComp(1).comp = selComp;
    
    if (strcmpi(sortComponentsStr, 'correlation between modalities'))
        selectedFeatureVal = (1:length(modalities));
    end
    
    % Output files of selected features
    outputFiles = outputFiles(selectedFeatureVal);
    
    set(handles, 'userdata', displayParameters);
    
    if (strcmpi(sortComponentsStr, 'correlation between modalities'))
        
        % Load parallel ICA fusion file
        load(parallel_ica_fusion_file);
        
        if (length(featureNames) == 2)
            
            % Get the correlation values between modalities
            sorted_values(1).val = displayParameters.avg_correlation;
            
            disp(['Sorting components based on correlation between ', featureNames{1}, ' and ', featureNames{2}, ' components ...']);
            [avg_corr, selComp] = sort(abs(sorted_values(1).val));
            selComp = selComp(end:-1:1);
            sorted_values(1).val = sorted_values(1).val(selComp);
            
            
            selectedComp(1).comp = selComp;
            selectedComp(2).comp = corrIndices(selectedComp(1).comp);
            
            
            numPara = 1;
            varStruct(numPara).tag = [featureNames{1}, ' Component'];
            varStruct(numPara).value = num2str(selectedComp(1).comp(:));
            
            numPara = numPara + 1;
            varStruct(numPara).tag = [featureNames{2}, ' Component'];
            varStruct(numPara).value = num2str(selectedComp(2).comp(:));
            
            numPara = numPara + 1;
            varStruct(numPara).tag = 'Correlation';
            varStruct(numPara).value = num2str(sorted_values(1).val(:), '%0.4f');
            
            corr_results_file  = fullfile(outputDir, [output_prefix, '_para_ica_correlations.txt']);
            ica_fuse_printToFile(corr_results_file, varStruct, ['Correlations between features ', featureNames{1}, ' and ', featureNames{2}], 'column_wise', 'append');
            
            clear varStruct;
            
        else
            corr_results_file  = fullfile(outputDir, [output_prefix, '_para_ica_correlations.txt']);
            fid = fopen(corr_results_file, 'w+');
            if (fid == -1)
                disp(['File ', corr_results_file, ' cannot be opened for writing results']);
            end
            fprintf(fid, [featureNames{1}, ' vs ', featureNames{2}, ' correlation: % 0.4f\n'], displayParameters.avg_correlation(1));
            fprintf(fid, [featureNames{2}, ' vs ', featureNames{3}, ' correlation: % 0.4f\n'], displayParameters.avg_correlation(2));
            fprintf(fid, [featureNames{1}, ' vs ', featureNames{3}, ' correlation: % 0.4f\n'], displayParameters.avg_correlation(3));
            fprintf(fid, [featureNames{1}, ' component no: %d\n'], displayParameters.corrIndices(1));
            fprintf(fid, [featureNames{2}, ' component no: %d\n'], displayParameters.corrIndices(2));
            fprintf(fid, [featureNames{3}, ' component no: %d\n'], displayParameters.corrIndices(3));
            fclose(fid);
            % a = 1;
            
            selectedComp(1).comp = displayParameters.corrIndices(1);
            selectedComp(2).comp = displayParameters.corrIndices(2);
            selectedComp(3).comp = displayParameters.corrIndices(3);
            
            images_per_figure = 1;
            
        end
        
        disp(['Correlation results are saved in file ', corr_results_file]);
        fprintf('\n\n');
        
        geneModality = strmatch('gene', modalities, 'exact');
        
        if (~isempty(geneModality) && ~isempty(locusFileName) && (strcmpi(type_parallel_ica, 'aa') || strcmpi(type_parallel_ica, 'aa-ref')))
            
            origLocusNames = textread(locusFileName, '%s', 'delimiter', '\n');
            
            for nGenes = 1:length(geneModality)
                
                locusNames = origLocusNames;
                
                % Mask indices
                geneMaskInd = displayParameters.mask_ind(geneModality(nGenes)).ind;
                
                % Number of SNPs
                numSNPS = length(geneMaskInd);
                
                % Check locus names
                if numSNPS ~= length(locusNames)
                    if (length(locusNames) < numSNPS)
                        fprintf('\nPlease check the text file %s \nas the number of locus names (%d) is less than the number of SNPs (%d)\n\n', ...
                            locusFileName, length(locusNames), numSNPS);
                        continue;
                    elseif (length(locusNames) > numSNPS)
                        % Update locus names
                        locusNames = locusNames(geneMaskInd);
                    end
                end
                
                
                %%%%%% PRINT DOMINANT SNPS %%%%%%%%%%%%%%%%%%%%%
                
                % Component files
                compFiles = str2mat(outputFiles(geneModality(nGenes)).name);
                compFiles = ica_fuse_fullFile('directory', outputDir, 'files', compFiles);
                
                
                selectedComponents = selectedComp(geneModality(nGenes)).comp;
                
                [selectedComponents] = uniqueVals(selectedComponents);
                
                % Selected Gene components
                compData = ica_fuse_loadData(compFiles, selectedComponents);
                
                compData = squeeze(compData(:, 2, :));
                
                %locusNames = displayParameters.locusNames;
                snps_file_name = fullfile(outputDir, [output_prefix, '_feature_', num2str(geneModality(nGenes)), ...
                    '_dominant_snps.txt']);
                
                % Open file in write mode
                fid = fopen(snps_file_name, 'w+');
                titleString = ['Dominant SNPs (Z-Threshold = ', num2str(SNP_Z_THRESHOLD), ') are as follows:'];
                fprintf(fid, '%s\n', titleString);
                fclose(fid);
                
                % Loop over selected components
                for nComp = 1:length(selectedComponents)
                    snp_comp = ica_fuse_returnFileIndex(selectedComponents(nComp));
                    titleToPrint = [featureNames{geneModality(nGenes)}, ' Comp ', snp_comp, ':'];
                    % Convert SNPs to Z-scores
                    snp_z_values = detrend(compData(:, nComp), 0);
                    snp_z_values = snp_z_values ./ std(snp_z_values);
                    % End for converting SNPs to Z-scores
                    snp_indices = find(abs(snp_z_values) >= SNP_Z_THRESHOLD);
                    snp_z_values = snp_z_values(snp_indices);
                    [dd, sorted_ind] = sort(abs(snp_z_values));
                    sorted_ind = sorted_ind(end:-1:1);
                    % Dominant SNP z-values and indices
                    snp_z_values = snp_z_values(sorted_ind);
                    snp_z_values = snp_z_values(:);
                    snp_indices = snp_indices(sorted_ind);
                    temp_locus_names = str2mat(locusNames(snp_indices));
                    
                    numPara = 1;
                    varStruct(numPara).tag = 'SNP';
                    varStruct(numPara).value = temp_locus_names;
                    
                    numPara = numPara + 1;
                    varStruct(numPara).tag = 'Z Score';
                    varStruct(numPara).value = snp_z_values;
                    
                    % Print to file
                    ica_fuse_printToFile(snps_file_name, varStruct, titleToPrint, 'column_wise', 'append');
                    
                    clear varStruct;
                end
                % End loop over selected components
                
                disp(['Dominant SNPS information for feature ', featureNames{geneModality(nGenes)}, ' is stored in file: ', snps_file_name]);
                
                fprintf('\n');
                
                %%%%%%% END FOR PRINTING DOMINANT SNPS %%%%%%%%%%%
                
            end
            
        end
        
        
    elseif (strcmpi(sortComponentsStr, 'two sample t-test on mixing coefficients'))
        
        if (numGroups == 1)
            error('No. of groups must be atleast two in order to two sample t-test on mixing coefficients');
        end
        
        selectedGroups = displayParameters.selectedGroups;
        
        disp(['Computing two sample t-test on mixing coefficients of feature ', selectedFeature, ' between ', deblank(groupNames(selectedGroups(1), :)), ...
            ' and ', deblank(groupNames(selectedGroups(2), :)), ' ...']);
        
        tmp_loadings = ica_fuse_loadData(fullfile(outputDir, outputFiles(1).loadingCoeffFile));
        
        g1_ind = ica_fuse_get_groupInd(selectedGroups(1), numSubjects);
        g2_ind = ica_fuse_get_groupInd(selectedGroups(2), numSubjects);
        
        p_values = zeros(1, size(tmp_loadings, 2));
        t_values = p_values;
        for nComp = 1:length(p_values)
            [p_values(nComp), t_values(nComp)] = ica_fuse_ttest2(tmp_loadings(g1_ind, nComp), tmp_loadings(g2_ind, nComp));
        end
        
        [p_values, sorted_inds] = sort(p_values);
        
        selectedComp(1).comp = sorted_inds;
        t_values = t_values(sorted_inds);
        
        numPara = 1;
        varStruct(numPara).tag = 'Component';
        varStruct(numPara).value = num2str(sorted_inds(:));
        
        numPara = numPara + 1;
        varStruct(numPara).tag = 'p-value';
        varStruct(numPara).value = num2str(p_values(:), '%0.4f');
        
        numPara = numPara + 1;
        varStruct(numPara).tag = 'T-value';
        varStruct(numPara).value = num2str(t_values(:), '%0.4f');
        
        ttest_results_file  = fullfile(outputDir, [output_prefix, '_para_ica_ttest2_results.txt']);
        
        % Print to file
        ica_fuse_printToFile(ttest_results_file, varStruct, ['Two sample t-test results of feature ', selectedFeature, ' between ', ...
            deblank(groupNames(selectedGroups(1), :)), ' and ', deblank(groupNames(selectedGroups(2), :))], 'column_wise', 'append');
        
        clear varStruct;
        
        disp(['Two sample t-test results are saved in file ', ttest_results_file]);
        fprintf('\n\n');
        
        
    end
    
    
    % Initialise variables
    cm = [];
    minInterval = [];
    maxInterval = [];
    colorbarLim = [];
    
    text_left_right = [];
    % Loop over features
    for nFeature = 1:length(selectedFeatureVal)
        
        checkModality = lower(deblank(modalities(selectedFeatureVal(nFeature), :)));
        
        % Create featureInfo structure
        featureInfo.feature_name = featureNames{selectedFeatureVal(nFeature)};
        featureInfo.groupNames = groupNames;
        featureInfo.numSubjects = numSubjects;
        
        % Component files
        compFiles = str2mat(outputFiles(nFeature).name);
        compFiles = ica_fuse_fullFile('directory', outputDir, 'files', compFiles);
        
        % Get the input files
        for nInputFiles = 1:length(dataInfo)
            
            if nInputFiles == 1
                feature_input_files = str2mat(dataInfo(nInputFiles).feature(selectedFeatureVal(nFeature)).files.name);
            else
                feature_input_files = str2mat(feature_input_files, ...
                    str2mat(dataInfo(nInputFiles).feature(selectedFeatureVal(nFeature)).files.name));
            end
            
        end
        % end for getting feature files
        
        if strcmpi(checkModality, 'fmri') || strcmpi(checkModality, 'smri')
            
            [compData, anatData, meanData, groupCompData, meanDataLegend, groupCompLegend, text_left_right, ...
                plotType, HInfo] = ica_fuse_loadCompData('component_files', compFiles, 'slices_in_mm', slices_in_mm, ...
                'anatomical_file', ANATOMICAL_FILE, 'anatomical_view', anatomical_plane, 'component_numbers', ...
                selectedComp(nFeature).comp, 'interp_message', ['Interpolating components of feature ', featureInfo.feature_name], ...
                'interp_title', ['Interp comp of ', featureInfo.feature_name], 'input_files', feature_input_files, 'voxels', voxels, ...
                'feature_info', featureInfo, 'mask_ind', mask_ind(selectedFeatureVal(nFeature)).ind, 'outputDir', ...
                outputDir, 'flip_analyze_images', flip_analyze_images);
            
            if ~isfield(displayParameters, 'text_left_right')
                displayParameters.text_left_right = text_left_right;
            end
            
            clear feature_input_files;
            clear maxICAIm; clear minICAIm;
            maxICAIm = zeros(1, length(selectedComp(nFeature).comp));
            minICAIm = zeros(1, length(selectedComp(nFeature).comp));
            minInterval = 0;
            maxInterval = 0;
            
            structDIM = HInfo.DIM;
            % Apply display parameters
            compData = ica_fuse_applyDispParameters(plotType, compData, convertToZ, imageValues, threshold);
            
            % create composite image
            [compData, minICAIm, maxICAIm, minInterval, maxInterval] = ica_fuse_create_composite(anatData, compData, ...
                imageValues, anatomical_plane, structDIM, threshold);
            colorbarLim = [minInterval, maxInterval];
            % Get the associated colormap
            cm = ica_fuse_getColormap(imageValues, 1);
            clear anatData;
            
        elseif (strcmpi(checkModality, 'gene') || strcmpi(checkModality, 'behavioral'))
            
            % Selected Gene components
            compData = ica_fuse_loadData(compFiles, selectedComp(nFeature).comp);
            
            compData = squeeze(compData(:, 2, :));
            
        elseif strcmpi(checkModality, 'eeg')
            
            % Selected Gene components
            compData = ica_fuse_loadData(compFiles, selectedComp(nFeature).comp);
            
        end
        
        % Loading coefficents
        loadingCoeffFile = fullfile(outputDir, outputFiles(nFeature).loadingCoeffFile);
        loadingCoeff = ica_fuse_loadData(loadingCoeffFile);
        
        loadingCoeff = detrend(loadingCoeff, 0); % Remove mean
        % Normalize loading coefficients
        for nL = 1:size(loadingCoeff, 2)
            current_data = loadingCoeff(:, nL);
            loadingCoeff(:, nL) = current_data./norm(current_data, 2);
        end
        % End for normalizing loading coefficients
        
        % Current components
        currentComps = selectedComp(nFeature).comp;
        loadingCoeff = loadingCoeff(:, currentComps);
        
        % Loop over selected components
        for nComp = 1:length(currentComps)
            compIndex = ica_fuse_returnFileIndex(currentComps(nComp));
            if strcmpi(checkModality, 'fmri') || strcmpi(checkModality, 'smri')
                plotData(nComp).feature(nFeature).data = squeeze(compData(nComp, :, :));
                plotData(nComp).feature(nFeature).loadingCoeff = loadingCoeff(:, nComp);
                plotData(nComp).feature(nFeature).feature_name = featureNames{selectedFeatureVal(nFeature)};
                plotData(nComp).feature(nFeature).colorbarLim = colorbarLim;
                plotData(nComp).feature(nFeature).minICAIm = minICAIm(nComp);
                plotData(nComp).feature(nFeature).maxICAIm = maxICAIm(nComp);
                plotData(nComp).feature(nFeature).text_left_right = displayParameters.text_left_right;
            elseif strcmpi(checkModality, 'eeg')
                plotData(nComp).feature(nFeature).data = compData(:, :, nComp);
                plotData(nComp).feature(nFeature).loadingCoeff = loadingCoeff(:, nComp);
                plotData(nComp).feature(nFeature).feature_name = featureNames{selectedFeatureVal(nFeature)};
                plotData(nComp).feature(nFeature).text_left_right = [];
                
                if (nComp == 1)
                    meanDataEEG = ica_fuse_loadData(feature_input_files);
                    meanDataEEG = mean(ica_fuse_remove_mean(squeeze(meanDataEEG(:, 2, :))), 2);
                    meanDataEEG = meanDataEEG(mask_ind(selectedFeatureVal(nFeature)).ind);
                    meanDataEEG = detrend(meanDataEEG(:), 0);
                end
                
                plotData(nComp).feature(nFeature).meanData = meanDataEEG;
                
            else
                plotData(nComp).feature(nFeature).data = compData(:, nComp);
                plotData(nComp).feature(nFeature).loadingCoeff = loadingCoeff(:, nComp);
                plotData(nComp).feature(nFeature).feature_name = featureNames{selectedFeatureVal(nFeature)};
                plotData(nComp).feature(nFeature).text_left_right = [];
            end
            
            plotData(nComp).feature(nFeature).modality = checkModality;
            plotData(nComp).feature(nFeature).titleStr = [featureNames{selectedFeatureVal(nFeature)}, ' ', compIndex];
            
            if nFeature == 1
                if (strcmpi(sortComponentsStr, 'correlation between modalities'))
                    if (length(featureNames) == 2)
                        axesTitle = ['Mixing Coefficients (Corr = ', num2str(sorted_values(1).val(nComp), '%0.4f'), ', T-val = ', num2str(r_to_t(sorted_values(1).val(nComp), sum(numSubjects)), '%0.4f'),  ')'];
                    else
                        axesTitle = 'Mixing coefficients';
                    end
                elseif (strcmpi(sortComponentsStr, 'two sample t-test on mixing coefficients'))
                    axesTitle = ['Mixing Coefficients (p-val=', num2str(p_values(nComp), '%0.4f'), ', T-val=', num2str(t_values(nComp), '%0.4f'), ')'];
                else
                    axesTitle = 'Mixing Coefficients';
                end
                
                plotData(nComp).axesTitle = axesTitle;
                plotData(nComp).groupNames = str2mat(dataInfo.name);
                plotData(nComp).numSubjects = numSubjects;
            end
        end
        % End loop over selected components
        
    end
    % End loop over features
    
    % Display parallel ICA components
    ica_fuse_display_parallel_ica_comp('plot_data', plotData, 'title_color', 'c', 'color_map', cm, 'time_course_color', 'm', 'number_per_figure', ...
        images_per_figure, 'slice_plane', anatomical_plane);
    
    set(handles, 'userdata', displayParameters);
    
    set(handles, 'pointer', 'arrow');
    
catch
    set(handles, 'pointer', 'arrow');
    ica_fuse_displayErrorMsg;
end


%%%%%%%%%%%%%%%%%%%% FUNCTION CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%

function display_defaults_callback(hObject, event_data, handles)
% get the display defaults

displayParameters = get(handles, 'userdata');

compListH = findobj(handles, 'tag', 'selComp');
controlPara = get(compListH, 'userdata');


[controlPara, displayParameters] = ica_fuse_get_display_parameters(displayParameters, controlPara, 'on');

set(compListH, 'userdata', controlPara);

set(handles, 'userdata', displayParameters);


function featureListCallback(hObject, event_data, handles)
% Feature listbox callback

% Listbox value
listVal = get(hObject, 'value');
listString = get(hObject, 'string');

selCompH = findobj(handles, 'tag', 'selComp');
compListVal = get(selCompH, 'value');

displayParameters = get(handles, 'userdata');

if (listVal ~= displayParameters.selectedFeatureVal)
    set(selCompH, 'value', 1);
    set(selCompH, 'string', displayParameters.featureListbox(listVal).string);
end

displayParameters.selectedFeatureVal = listVal;

set(handles, 'userdata', displayParameters);



function sortPopupCallback(hObject, event_data, handles)
% Sort popup callback


ud = get(handles, 'userdata');
groupNames = ud.groupNames;

% Get selection string
getVal = get(hObject, 'value');
getString = get(hObject, 'string');

% Component Listbox
compListH = findobj(handles, 'tag', 'selComp');

% Feature Listbox
selFeatureH = findobj(handles, 'tag', 'selFeature');

selectedStr = deblank(getString(getVal, :));

set(selFeatureH, 'enable', 'on');
set(compListH, 'enable', 'on');

if ~strcmpi(selectedStr, 'none')
    
    if (strcmpi(selectedStr, 'two sample t-test on mixing coefficients'))
        if (size(groupNames, 1) < 2)
            error('No. of groups must be atleast two in order to two sample t-test on mixing coefficients');
        end
        
        if (size(groupNames, 1) > 2)
            selectedGroups = ica_fuse_listdlg('PromptString', 'Select 2 groups', 'SelectionMode', 'multiple',...
                'ListString', groupNames, 'movegui', 'center', 'windowStyle', 'modal', 'title_fig', 'Select groups');
            if (length(selectedGroups) < 2)
                error('No. of groups must be atleast two in order to two sample t-test on mixing coefficients');
            end
            
        else
            selectedGroups = [1, 2];
        end
        
        ud.selectedGroups = selectedGroups;
        set(handles, 'userdata', ud);
        set(compListH, 'enable', 'off');
        
    else
        set(selFeatureH, 'enable', 'off');
        set(compListH, 'enable', 'off');
    end
    
end


function selectLocusFile(hObject, event_data, handles, selectFile)
% Select locus file containing SNP names

if ~exist('selectFile', 'var')
    selectFile = 1;
end

dispParameters = get(handles, 'userdata');

parallel_ica_fusion_file = dispParameters.parallel_ica_fusion_file;

load(parallel_ica_fusion_file);

if selectFile
    % Select locus file name
    locusFileName = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', '*.txt', 'title', ...
        'Select SNPs locus file ...');
else
    locusFileName = paraICAInfo.run_analysis.locusFileName;
end

if isempty(locusFileName)
    error('Locus file containing SNPs is not selected');
end

if ~exist(locusFileName, 'file')
    error('Error:LocusFile', 'Locus file %s doesn''t exist', locusFileName);
end

% Save parallel ICA fusion information
paraICAInfo.run_analysis.locusFileName = locusFileName;
ica_fuse_save(parallel_ica_fusion_file, 'paraICAInfo');


function [outputVec] = uniqueVals(inputVec)

indices = zeros(1, length(inputVec));
for nV = 1:length(inputVec)
    ind = find(inputVec == inputVec(nV));
    indices(nV) = ind(1);
end

[indices] = unique(indices);
outputVec = inputVec(indices);


function tstat = r_to_t(r, n)
% correlation value to t stat

tstat = r/sqrt((1 - r^2)/(n - 2));



function helpCallback(hObject, event_data, handles)
% Help callback

msgStr = char('1. Correlation between modalities - Components are sorted based on the correlation measure between the modalities. Each component displayed consists of mixing coefficients of each modality and source of each modality.', ...
    '', ...
    '2. Two sample t-test on mixing coefficients - Two sample t-test is done on the mixing coefficients of the selected modality.');
ica_fuse_dialogBox('title', 'Sorting Components', 'textBody', msgStr, 'textType', 'large');


function resultsSummaryCallback(hObject, event_data, handles)
% Results summary callback
%

formatName = questdlg('Select results format', 'Results format', 'HTML', 'PDF', 'HTML');

if (~isempty(formatName))
    
    displayParameters = get(handles, 'userdata');
    param_file = displayParameters.parallel_ica_fusion_file;
    
    load(param_file);
    
    if (~exist('paraICAInfo', 'var'))
        error('Selected file is not a valid parallel fusion MAT file');
    end
    
    outDir = fullfile(fileparts(param_file), [paraICAInfo.setup_analysis.prefix, '_pica_results']);
    opts.outputDir = outDir;
    opts.showCode = false;
    opts.useNewFigure = false;
    opts.format = lower(formatName);
    opts.createThumbnail = true;
    if (strcmpi(opts.format, 'pdf'))
        opt.useNewFigure = false;
    end
    assignin('base', 'param_file', param_file);
    assignin('base', 'displayParameters', displayParameters);
    opts.codeToEvaluate = 'ica_fuse_pica_summary(param_file, displayParameters);';
    %publish('icatb_gica_html_report', 'outputDir', outDir, 'showCode', false, 'useNewFigure', false);
    disp('Generating reults summary. Please wait ....');
    drawnow;
    publish('ica_fuse_pica_summary', opts);
    
    close all;
    
    if (strcmpi(opts.format, 'html'))
        ica_fuse_openHTMLHelpFile(fullfile(outDir, 'ica_fuse_pica_summary.html'));
    else
        open(fullfile(outDir, 'ica_fuse_pica_summary.pdf'));
    end
    
    disp('Done');
    
end
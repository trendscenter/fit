function ica_fuse_anyway_display(param_file, displayParameters)
% Anyway ICA report generator

ica_fuse_defaults;
global ANATOMICAL_FILE;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGE_VALUES;
global ANATOMICAL_PLANE;
global FIG_FG_COLOR;


if (~isempty(param_file))
    load(param_file);
end

outputDir = fileparts(param_file);

if (isempty(outputDir))
    outputDir = pwd;
end


try
    imageValues = displayParameters.image_values;
catch
    imageValues = IMAGE_VALUES;
end

try
    threshold = abs(displayParameters.z_threshold);
catch
    threshold = Z_THRESHOLD;
end

try
    anatomical_plane = lower(displayParameters.anatomical_plane);
catch
    anatomical_plane = ANATOMICAL_PLANE;
end

try
    slices_in_mm = displayParameters.slices_in_mm;
catch
    slices_in_mm = -40:4:72;
end

try
    convertToZ = displayParameters.convert_to_z;
catch
    convertToZ = strcmpi('yes', CONVERT_TO_Z);
end

convertToZStr = 'no';
if (convertToZ)
    convertToZStr = 'yes';
end

disp_opts = {'positive and negative', 'positive', 'absolute value', 'negative'};
if (isnumeric(imageValues))
    imageValueStr = disp_opts{imageValues};
else
    imageValueStr = imageValues;
end

sz = get(0, 'ScreenSize');
defaultFigPos = get(0, 'DefaultFigurePosition');
defaultFigPos(3) = 0.9*sz(3);
defaultFigPos(2) = 0.5*sz(4) - 0.5*defaultFigPos(4) + 10;
defaultFigPos(1) = 0.5*sz(3) - 0.5*defaultFigPos(3) + 10;

% Get params from the output files
%slices_in_mm = (-40:4:72);
%threshold = threshold;
feature_info = fusionInfo.setup_analysis.feature_info;
convert_to_z = convertToZStr;
imageValues = lower(imageValueStr);
%anatomical_plane = ANATOMICAL_PLANE;
anatomical_file = ANATOMICAL_FILE;
outputFiles = fusionInfo.run_analysis.outputFiles;
modalities = cellstr(char(feature_info.modality));
featureNames = cellstr(char(feature_info.feature_name));
numComp = [feature_info.comp];
numSubjects = fusionInfo.run_analysis.numSubjects;

fusionInfo.run_analysis.numSubjects = numSubjects;
fusionInfo.run_analysis.threshold = threshold;
fusionInfo.run_analysis.convert_to_z = convert_to_z;
fusionInfo.run_analysis.anatomical_plane = anatomical_plane;
fusionInfo.run_analysis.anatomical_file = anatomical_file;
fusionInfo.run_analysis.slices_in_mm =slices_in_mm;
fusionInfo.run_analysis.image_values = imageValues;
mask_ind = fusionInfo.run_analysis.mask_ind;
[dims, voxels] = ica_fuse_getFeatureDIM(mask_ind);



feature_files = cell(1, length(featureNames));
loadingCoeff = cell(1, length(featureNames));
for nF = 1:length(featureNames)
    
    feature_files{nF} = ica_fuse_rename_4d_file(outputFiles(nF).name);
    loadingCoeff{nF} = ica_fuse_loadData(fullfile(outputDir, outputFiles(nF).loadingCoeffFile));
end


%% Anyway ICA Fusion Parameters
dispParaStr = dispPara(fusionInfo);

fprintf('\n\n');
disp(char(dispParaStr));
fprintf('\n\n');


%%

%% Anyway ICA components
% Mixing coefficients and the source maps are displayed.
%

meanEEG = cell(1, length(featureNames));
for nFea = 1:length(featureNames)
    
    if (strcmpi(modalities{nFea}, 'eeg'))
        input_data = ica_fuse_loadData(fusionInfo.setup_analysis.feature_info(nFea).files);
        meanData = mean(squeeze(input_data(fusionInfo.run_analysis.mask_ind(nFea).ind, 2, :)), 2);
        meanEEG{nFea} = meanData;
    end
end


for nFea = 1:length(feature_info)
    
    for nComp = 1:feature_info(nFea).comp
        axesTitle = ['Mixing coeff Feature ', num2str(nFea), ' Comp ', num2str(nComp)];
        gH = ica_fuse_getGraphics(axesTitle, 'timecourse', '', 'on');
        sh = get(gH, 'currentAxes');
        set(gH, 'resize', 'on');
        tmpMixingCoeff = loadingCoeff{nFea}(:, nComp);
        sh = axes('units', 'normalized', 'parent', gH, 'position', [0.1, 0.1, 0.8, 0.8]);
        plot(tmpMixingCoeff, 'parent', sh, 'color', 'm', 'linewidth', 2);
        axis(sh, 'tight');
        title(axesTitle, 'parent', sh, 'color', FIG_FG_COLOR);
        % ica_fuse_plotParaICALoading('data', tmpMixingCoeff, 'parent', gca, 'num_features', length(featureNames), 'colors', allColors, 'legend_string', featureNames, 'titlecolor', FIG_FG_COLOR);
        xlabel('Subjects', 'parent', sh);
        ylabel('Z-scores', 'parent', sh);
        set(sh, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
        
        modality1_comp = fullfile(outputDir, deblank(feature_files{nFea}(nComp, :)));
        title_str = [featureNames{nFea},  ' Comp ', ica_fuse_returnFileIndex(nComp)];
        if (strcmpi(modalities{nFea}, 'smri') || strcmpi(modalities{nFea}, 'fmri'))
            % fMRI or sMRI
            gH = ica_fuse_getGraphics(title_str, 'graphics', '', 'on');
            set(gH, 'resize', 'on');
            ica_fuse_image_viewer(modality1_comp, 'structfile', anatomical_file, 'threshold', threshold, 'convert_to_zscores', convert_to_z, 'image_values', imageValues, ...
                'slices_in_mm', slices_in_mm, 'anatomical_view', anatomical_plane, 'labels', title_str, 'axesh', gca);
        elseif strcmpi(modalities{nFea}, 'eeg')
            % EEG
            %title_str = [featureNames{nFea},  ' Comp ', ica_fuse_returnFileIndex(nComp)];
            tmp_dat = ica_fuse_loadData(modality1_comp);
            tmpMeanData = meanEEG{nFea};
            interpFactor = ceil(voxels / length(tmpMeanData));
            if (interpFactor > 0)
                % resample mean data
                tmpMeanData = ica_fuse_resample(tmpMeanData, voxels, length(tmpMeanData));
                % resample component data
                tmp_dat = ica_fuse_resampleCompData(tmp_dat, length(tmpMeanData), voxels, size(tmp_dat, 1));
            end
            
            gH = ica_fuse_getGraphics(title_str, 'timecourse', '', 'on');
            set(gH, 'resize', 'on');
            ica_fuse_plotTimecourse('parent', gca, 'data', tmp_dat, 'color', 'c', 'titleColor', FIG_FG_COLOR, 'title', title_str, 'YAxisLocation', 'left', 'meanData', ...
                tmpMeanData, 'meanDataLegend', {'Mean'},  'groupCompData',[], 'groupCompLegend', {});
            legend('Mean', 'Component', 'Location', 'best');
        else
            % Gene
            title_str = [featureNames{nFea},  ' Comp ', ica_fuse_returnFileIndex(nComp)];
            tmp_dat = ica_fuse_loadData(modality1_comp);
            tmp_dat = squeeze(tmp_dat(:, 2, :));
            gH = ica_fuse_getGraphics(title_str, 'timecourse', '', 'on');
            set(gH, 'resize', 'on');
            ica_fuse_plotSNP('data',tmp_dat, 'parent', gca, 'xlabel', 'SNPs', 'ylabel', ...
                'Z-scores', 'titleColor', FIG_FG_COLOR, 'color', 'm', 'title', title_str);
        end
        
    end
    
end



function dispParaStr = dispPara(paraICAInfo)


dispParaStr{1} = '....................................................';
dispParaStr{end + 1} = ['Feature names: ', ica_fuse_formatStr(cellstr(char(paraICAInfo.setup_analysis.feature_info.feature_name)), ',')];
dispParaStr{end + 1} = ['Modalities: ', ica_fuse_formatStr(cellstr(char(paraICAInfo.setup_analysis.feature_info.modality)), ',')];
dispParaStr{end + 1} = ['Number of components: ', num2str([paraICAInfo.setup_analysis.feature_info.comp])];
dispParaStr{end + 1} = ['Anatomical file: ', paraICAInfo.run_analysis.anatomical_file];
dispParaStr{end + 1} = ['Slice Plane: ', upper(paraICAInfo.run_analysis.anatomical_plane(1)), paraICAInfo.run_analysis.anatomical_plane(2:end)];
dispParaStr{end + 1} = ['Image values: ', upper(paraICAInfo.run_analysis.image_values(1)), paraICAInfo.run_analysis.image_values(2:end)];
dispParaStr{end + 1} = ['Convert to Z-scores: ', paraICAInfo.run_analysis.convert_to_z];
dispParaStr{end + 1} = ['Threshold: ', num2str( paraICAInfo.run_analysis.threshold)];
dispParaStr{end + 1} = '....................................................';

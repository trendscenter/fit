function ica_fuse_cmica_display(param_file, displayParameters)
%% Joint cmICA display
%

ica_fuse_defaults;
global ANATOMICAL_FILE;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGE_VALUES;
global ANATOMICAL_PLANE;


if (~exist('param_file', 'var'))
    param_file = ica_fuse_selectEntry('title', 'Select Joint cmICA fusion MAT file', 'filter', '*cmica_info.mat');
    drawnow;
end


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
convert_to_z = convertToZStr;
imageValues = lower(imageValueStr);
%anatomical_plane = ANATOMICAL_PLANE;
anatomical_file = ANATOMICAL_FILE;
outputFiles = cmICAInfo.outputFiles;
featureNames = cellstr(char(feature_info.feature_name));
numComp = 2*cmICAInfo.numComp;

cmICAInfo.run_analysis.threshold = threshold;
cmICAInfo.run_analysis.convert_to_z = convert_to_z;
cmICAInfo.run_analysis.anatomical_plane = anatomical_plane;
cmICAInfo.run_analysis.anatomical_file = anatomical_file;
cmICAInfo.run_analysis.slices_in_mm =slices_in_mm;
cmICAInfo.run_analysis.image_values = imageValues;

%% cmICA Fusion Parameters
dispParaStr = dispPara(cmICAInfo);

fprintf('\n\n');
disp(char(dispParaStr));
fprintf('\n\n');


%%

%% Joint cmICA components
for nFea = 1:length(feature_info)
    
    for nComp = 1:numComp
        
        
        modality1_comp = fullfile(outputDir, [deblank(outputFiles(nFea).name), ',', num2str(nComp)]);
        modality1_conn = fullfile(outputDir, [deblank(outputFiles(nFea).conn_name), ',', num2str(nComp)]);
        feature_name = featureNames{nFea};
        title_str = [feature_name,  ' Comp ', ica_fuse_returnFileIndex(nComp)];
        
        % fMRI or sMRI
        gH = ica_fuse_getGraphics(title_str, 'graphics', '', 'on');
        set(gH, 'resize', 'on');
        ica_fuse_image_viewer(modality1_comp, 'structfile', anatomical_file, 'threshold', threshold, 'convert_to_zscores', convert_to_z, 'image_values', imageValues, ...
            'slices_in_mm', slices_in_mm, 'anatomical_view', anatomical_plane, 'labels', title_str, 'axesh', gca);
        drawnow;
        
        title_str = [feature_name,  ' Connectivity ', ica_fuse_returnFileIndex(nComp)];
        gH = ica_fuse_getGraphics(title_str, 'graphics', '', 'on');
        set(gH, 'resize', 'on');
        ica_fuse_image_viewer(modality1_conn, 'structfile', anatomical_file, 'threshold', threshold, 'convert_to_zscores', convert_to_z, 'image_values', imageValues, ...
            'slices_in_mm', slices_in_mm, 'anatomical_view', anatomical_plane, 'labels', title_str, 'axesh', gca);
        drawnow;
        
    end
    
end



function dispParaStr = dispPara(cmICAInfo)


dispParaStr{1} = '....................................................';
dispParaStr{end + 1} = ['Feature names: ', ica_fuse_formatStr(cellstr(char(cmICAInfo.dataInfo.feature_name)), ',')];
dispParaStr{end + 1} = ['Modalities: ', ica_fuse_formatStr(cellstr(char(cmICAInfo.dataInfo.modality)), ',')];
dispParaStr{end + 1} = ['Number of components: ', num2str([cmICAInfo.numComp])];
dispParaStr{end + 1} = ['Anatomical file: ', cmICAInfo.run_analysis.anatomical_file];
dispParaStr{end + 1} = ['Slice Plane: ', upper(cmICAInfo.run_analysis.anatomical_plane(1)), cmICAInfo.run_analysis.anatomical_plane(2:end)];
dispParaStr{end + 1} = ['Image values: ', upper(cmICAInfo.run_analysis.image_values(1)), cmICAInfo.run_analysis.image_values(2:end)];
dispParaStr{end + 1} = ['Convert to Z-scores: ', cmICAInfo.run_analysis.convert_to_z];
dispParaStr{end + 1} = ['Threshold: ', num2str( cmICAInfo.run_analysis.threshold)];
dispParaStr{end + 1} = '....................................................';

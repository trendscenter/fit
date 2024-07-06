function ica_fuse_anyway_fusion_summary(param_file, displayParameters)

ica_fuse_defaults;
global IMAGE_VALUES;
global Z_THRESHOLD;
global ANATOMICAL_PLANE;
global CONVERT_TO_Z;


if (~exist('param_file', 'var'))
    param_file = ica_fuse_selectEntry('title', 'Select anyway fusion MAT file', 'filter', '*anyway_fusion.mat');
    drawnow;
end

try
    imageValues = displayParameters.image_values;
catch
    imageValues = IMAGE_VALUES;
end

try
    threshold = abs(displayParameters.threshold);
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
    convertToZ = CONVERT_TO_Z;
end

if (ischar(convertToZ))
    convertToZ = strcmpi('yes', convertToZ);
end

disp_opts = {'positive and negative', 'positive', 'absolute value', 'negative'};
if (isnumeric(imageValues))
    imageValueStr = disp_opts{imageValues};
else
    imageValueStr = imageValues;
end


displayParameters.image_values = imageValueStr;
displayParameters.z_threshold = threshold;
displayParameters.anatomical_plane = anatomical_plane;
displayParameters.slices_in_mm = slices_in_mm;
displayParameters.convert_to_z = convertToZ;

formatName = 'html';
try
    formatName = displayParameters.format;
catch
end

opts.format = formatName;

load(param_file);
outDir = fullfile(fileparts(param_file), [fusionInfo.setup_analysis.prefix, '_anyway_fusion_results']);
opts.outputDir = outDir;
opts.showCode = false;
opts.useNewFigure = false;
opts.format = lower(formatName);
opts.createThumbnail = true;
if (strcmpi(opts.format, 'pdf'))
    opts.useNewFigure = false;
end
assignin('base', 'param_file', param_file);
assignin('base', 'displayParameters', displayParameters);
opts.codeToEvaluate = 'ica_fuse_anyway_display(param_file, displayParameters);';
disp('Generating reults summary. Please wait ....');
drawnow;
publish('ica_fuse_anyway_display', opts);

close all;


if (strcmpi(opts.format, 'html'))
    ica_fuse_openHTMLHelpFile(fullfile(outDir, 'ica_fuse_anyway_display.html'));
else
    open(fullfile(outDir, 'ica_fuse_anyway_display.pdf'));
end

disp('Done');
fprintf('\n');

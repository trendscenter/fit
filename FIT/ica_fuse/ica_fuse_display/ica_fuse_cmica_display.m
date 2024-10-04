function ica_fuse_cmica_display(param_file, displayParameters)
%% Joint cmICA display
%

ica_fuse_defaults;
global ANATOMICAL_FILE;
global CONVERT_TO_Z;
global Z_THRESHOLD;
global IMAGE_VALUES;
global ANATOMICAL_PLANE;


useGUI = 0;
if (~exist('param_file', 'var'))
    param_file = ica_fuse_selectEntry('title', 'Select Joint cmICA fusion MAT file', 'filter', '*cmica_info.mat');
    useGUI = 1;
    drawnow;
end


if (~isempty(param_file))
    load(param_file);
end

outputDir = fileparts(param_file);


if (useGUI)
    displayParameters = ica_fuse_cmICA_display_defaults;
    drawnow;
end

if (isempty(outputDir))
    outputDir = pwd;
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
    convertToZ = strcmpi('yes', CONVERT_TO_Z);
end

anatomical_file = ANATOMICAL_FILE;
try
    anatomical_file = displayParameters.structFile;
catch
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
%anatomical_file = ANATOMICAL_FILE;
feature_info = cmICAInfo.dataInfo;
outputFiles = cmICAInfo.outputFiles;
featureNames = cellstr(char(feature_info.feature_name));
numComp = 2*cmICAInfo.numComp;

cmICAInfo.run_analysis.threshold = threshold;
cmICAInfo.run_analysis.convert_to_z = convert_to_z;
cmICAInfo.run_analysis.anatomical_plane = anatomical_plane;
cmICAInfo.run_analysis.anatomical_file = anatomical_file;
cmICAInfo.run_analysis.slices_in_mm =slices_in_mm;
cmICAInfo.run_analysis.image_values = imageValues;

modalities = cellstr(char(cmICAInfo.dataInfo.modality));
chkDTIFMRI = 0;

if (length(modalities) == 2)
    if any(strcmpi('fmri', modalities)) && any(strcmpi('dti', modalities))
        chkDTIFMRI = 1;
    end
end



%% cmICA Fusion Parameters
dispParaStr = dispPara(cmICAInfo);

fprintf('\n\n');
disp(char(dispParaStr));
fprintf('\n\n');


%%

%% Joint cmICA components
if (chkDTIFMRI)
    
    ica_fuse_joint_cmICA_display(cmICAInfo, imageValues, convertToZStr, threshold, anatomical_file, anatomical_plane, slices_in_mm);
    
else
    
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




function ica_fuse_joint_cmICA_display(cmICAInfo, image_values, convert_to_zscores, threshold, structFile, slice_plane, slices_in_mm)
% Lei: display joint-cmICA results!
% July, 2024
%%%%%%%%%

resultsdir = cmICAInfo.outputDir;
prefix = cmICAInfo.output_prefix;
modalities = cellstr(char(cmICAInfo.dataInfo.modality));

fmri_ind = strmatch('fmri', lower(modalities), 'exact');
dti_ind = strmatch('dti', lower(modalities), 'exact');

compFiles = ica_fuse_rename_4d_file(ica_fuse_fullFile('directory', resultsdir, 'files', [prefix, '_agg_joint_cmica_comps2.nii']));
compFilesSd = ica_fuse_rename_4d_file(ica_fuse_fullFile('directory', resultsdir, 'files', cmICAInfo.outputFiles(dti_ind).name));
compFilesSf = ica_fuse_rename_4d_file(ica_fuse_fullFile('directory', resultsdir, 'files', cmICAInfo.outputFiles(fmri_ind).name));
compFilesRd = ica_fuse_rename_4d_file(ica_fuse_fullFile('directory', resultsdir, 'files', cmICAInfo.outputFiles(dti_ind).conn_name));
compFilesRf = ica_fuse_rename_4d_file(ica_fuse_fullFile('directory', resultsdir, 'files', cmICAInfo.outputFiles(fmri_ind).conn_name));

components = [1:size(compFiles,1)]';
figVisible = 'on';

% structFile = fullfile(fileparts(which('groupica.m')), 'ica_fuse_templates', 'ch2bet_3x3x3.nii');
% convert_to_zscores = 'Yes';
% threshold = 1;
% slices_in_mm = [-40:6:72];
% slice_plane = 'Axial';


%image_values = 'positive';
load ica_fuse_colors coldhot coldhot_sensitive;
returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');

%%
for nF = 1:size(compFiles, 1)
    %%
    cn = deblank(compFiles(nF, :));
    cn1 = deblank(compFilesRd(nF, :));
    cn2 = deblank(compFilesRf(nF, :));
    
    [dd, fN, extn] = fileparts(cn);
    
    gH = ica_fuse_getGraphics([fN, extn, '(Comp ', ica_fuse_returnFileIndex(components(nF)), ')'], 'graphics', 'imviewer', figVisible);
    set(gH, 'resize', 'on');
    %
    
    xOffSet = 0.1;
    yOffSet = 0.05;
    
    width = 0.3;
    height = 0.3;
    xPos = 0.5-(width/2);
    yPos = 0.65;
    sh = axes('parent', gH, 'units', 'normalized', 'position', [xOffSet, yPos, width, height]);
    hD=plotStackedOrtho_cmica(cn, 'structfile', structFile, 'image_values', image_values,'convert_to_zscores', convert_to_zscores, 'threshold', threshold, 'set_to_max_voxel', 1, ...
        'axesh', sh, 'labels', {['Shared Source Map S (Comp ', ica_fuse_returnFileIndex(components(nF)), ')']; 'Peak Coordinates (mm)'}, 'colorbar', true, 'colorbar_label', true);
    
    %
    width2 = 0.4;
    height2 = 0.4;
    sh = axes('parent', gH, 'units', 'normalized', 'position', [0.05, 0.3-(height2/2), width2, height2]);
    ica_fuse_image_viewer(cn1, 'structfile', structFile, 'image_values', 'positive dmri', 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, ...
        'slices_in_mm', slices_in_mm, 'anatomical_view', slice_plane, 'axesh', sh, 'colorbar', 1, 'labels', ...
        {['Structural Connectivity R (Conn ', ica_fuse_returnFileIndex(components(nF)), ')']}, 'colorbar_orientation','horiz');
    
  
    %
    sh = axes('parent', gH, 'units', 'normalized', 'position', [0.5+0.05, 0.3-(height2/2), width2, height2]);
    ica_fuse_image_viewer(cn2, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, ...
        'slices_in_mm', slices_in_mm, 'anatomical_view', slice_plane, 'axesh', sh, 'colorbar', 1, 'labels', ...
        {['Functional Connectivity R (Conn ', ica_fuse_returnFileIndex(components(nF)), ')']}, 'colorbar_orientation','horiz');
    %
    sh1 = axes('parent', gH, 'units', 'normalized', 'position', [xOffSet+0.5, yPos, width, height]);
    files = [{compFilesSf(nF,:)};{compFilesSd(nF,:)}];
    ica_fuse_display_composite_cmica(files, structFile,hD.coords,sh1);
    drawnow;
end



function [hD]=plotStackedOrtho_cmica(file_name, varargin)

ica_fuse_defaults;
global FIG_FG_COLOR;

useColorbar = 1;
threshold = 1;
convert_to_zscores = 'no';
image_values = 'positive';
labels = '';
cmap = hot(256); cmap = cmap(1:4:end,:); %cmap=hot(64);
structFile = fullfile (fileparts(which('gift.m')), 'ica_fuse_templates', 'ch2bet_3x3x3.nii');

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'structfile'))
        structFile = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'image_values'))
        image_values = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'cmap'))
        cmap = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'threshold'))
        threshold = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'convert_to_zscores'))
        convert_to_zscores = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'axesh'))
        sh = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'colorbar'))
        useColorbar = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'labels'))
        labels = varargin{n + 1};
    end
end

structVol = ica_fuse_spm_vol(structFile);
structVol = structVol(1);

hD = ica_fuse_orth_views(file_name, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, 'set_to_max_voxel', 1, ...
    'get_interp_data', 1);
if (~exist('sh', 'var'))
    sh = gca;
end
[sliceXY, sliceXZ, sliceYZ] = returnSlices(hD.data, hD.maxVoxelPos);
tmp = stackData({sliceYZ, sliceXZ, sliceXY}, 101);
plotImage(sh, tmp, [1, 200], structVol);
%imagesc(tmp, [1, 200]);
colormap(hD.cmap);
%axis(sh, 'image');
%set(sh, 'Ytick', []);
%set(sh, 'Xtick', [])
realCoords = (hD.maxVoxelPos - hD.voxelOrigin).*hD.VOX;
hD.coords = realCoords;

str = [labels; ['(', num2str(realCoords(1)), ',', num2str(realCoords(2)), ',', num2str(realCoords(3)), ')'] ];
title(str, 'parent', sh, 'horizontalalignment', 'center', 'fontweight', 'bold');
if (useColorbar)
    ch = colorbar('location', 'eastoutside');
    ch.Position(1) = ch.Position(1)+0.08;
    set(ch, 'xLim', [1, 100]);
    minICAVal = hD.minICAIM;
    maxICAVal = hD.maxICAIM;
    set(ch, 'xTick', [1, 100]);
    set(ch, 'xTicklabel', cellstr(char(num2str(minICAVal, '%0.1f'), num2str(maxICAVal, '%0.1f'))));
    set(ch, 'XColor', FIG_FG_COLOR, 'YColor', FIG_FG_COLOR);
end




function [sliceXY, sliceXZ, sliceYZ] = returnSlices(data, voxelcoords)

sliceXY = rot90(reshape(data(:, :, voxelcoords(end)), size(data, 1), size(data, 2)));
sliceXZ = rot90(reshape(data(:, voxelcoords(2), :), size(data, 1),size(data, 3)));
sliceYZ = rot90(reshape(data(voxelcoords(1), :, :), size(data, 2), size(data, 3)));

function data = stackData(slices, minVal)

[m1, n1] = size(slices{1});
[m2, n2] = size(slices{2});
[m3, n3] = size(slices{3});

maxSizeX = max([m1, m2, m3]);
maxSizeY = max([n1, n2, n3]);

data = minVal*ones(maxSizeX, n1 + n2 + n3);

e = 0;
for nS = 1:length(slices)
    tmp = slices{nS};
    s = e + 1;
    e = e + size(tmp, 2);
    inda = ceil(maxSizeX/2) - ceil(size(tmp, 1)/2) + 1;
    indb = inda + size(tmp, 1) - 1;
    data(inda:indb, s:e) = tmp;
end

function cmap = getCmap(image_values)

load ica_fuse_colors coldhot coldhot_sensitive;

returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');

if (returnValue == 1)
    coldhot_sensitive = coldhot;
    cmap = coldhot_sensitive(1:4:end, :);
elseif (returnValue == 4)
    cmap = coldhot_sensitive(1:128, :);
    cmap = cmap(1:2:end, :);
else
    cmap = coldhot_sensitive(129:end, :);
    cmap = cmap(1:2:end, :);
end


function subH = plotImage(subH, data, CLIM, structVol)
% Function to plot the image at the specified position
%

%structVol = get(get(subH, 'parent'), 'userdata');
%structVol = tmpV.structVol;
VOX = double(structVol(1).private.hdr.pixdim(2:4));
imageAxisPos = get(subH, 'position');
image(data, 'parent', subH, 'CDataMapping', 'scaled');
set(subH, 'units', 'normalized');
oldWidth = imageAxisPos(3);
oldHeight = imageAxisPos(4);
yAxisRatio = imageAxisPos(4)/imageAxisPos(3);
xAxisRatio = imageAxisPos(3)/imageAxisPos(4);
if(yAxisRatio>1)
    yAxisRatio = 1;
else
    xAxisRatio = 1;
end
newWidth = imageAxisPos(3)*yAxisRatio;
newHeight = imageAxisPos(4)*xAxisRatio;

if (newHeight < oldHeight)
    imageAxisPos(2) = imageAxisPos(2) + 0.5*(oldHeight - newHeight);
end

imageAxisPos = [imageAxisPos(1), imageAxisPos(2), newWidth, newHeight];
set(subH, 'position', imageAxisPos);
%setImagePos(subH, imagePos);
set(subH, 'clim', CLIM); % set the axis positions to the specified
if (all(VOX == VOX(1)))
    axis(subH, 'image');
end
set(subH, 'XTick', []);
set(subH, 'XTickLabel', []);
set(subH, 'YTick', []);
set(subH, 'YTickLabel', []);



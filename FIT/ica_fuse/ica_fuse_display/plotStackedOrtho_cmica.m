function [hD]=plotStackedOrtho(file_name, varargin)

icatb_defaults;
global FONT_COLOR;

useColorbar = 1;
threshold = 1;
convert_to_zscores = 'no';
image_values = 'positive';
labels = '';
cmap = hot(256); cmap = cmap(1:4:end,:); %cmap=hot(64);
structFile = fullfile (fileparts(which('gift.m')), 'icatb_templates', 'ch2bet_3x3x3.nii');

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

hD = icatb_orth_views(file_name, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, 'set_to_max_voxel', 1, ...
    'get_interp_data', 1, 'cmap', cmap);
if (~exist('sh', 'var'))
    sh = gca;
end
[sliceXY, sliceXZ, sliceYZ] = returnSlices(hD.data, hD.maxVoxelPos);
tmp = stackData({sliceYZ, sliceXZ, sliceXY}, 101);
imagesc(tmp, [1, 200]);
colormap(hD.cmap);
axis(sh, 'image');
set(sh, 'Ytick', []);
set(sh, 'Xtick', [])
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
    set(ch, 'XColor', FONT_COLOR, 'YColor', FONT_COLOR);
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

load icatb_colors coldhot coldhot_sensitive;

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

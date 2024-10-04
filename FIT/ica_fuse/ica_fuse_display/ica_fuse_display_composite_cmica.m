function  varargout = ica_fuse_display_composite_cmica(files, structFile,coords,sh )
%% Make composite plots for joint cmICA.

c=ica_fuse_display_composite(files,'anatomical_file', structFile,'coords',coords);%,vh,'coords',hD.maxVoxelPos);


c.modality = [{'Functional S'},{'Structural S'}];
%
imagesc(sh,c.slices, [1, 192]);
colormap(sh,c.cmap);
axis(sh, 'image');
axis(sh, 'off');
str = ['Structural vs. Functional Source S'];
title(str, 'parent', sh, 'horizontalalignment', 'center', 'fontweight', 'bold');
%
ch = colorbar(sh,'horiz');
pos = get(ch, 'position');
offset = 0.;
maxWidth = pos(3);
cwdiths = (maxWidth-offset*length(files))/(length(files));
pos(2) = pos(2)-0.03;
pos(3) = cwdiths;
cMax = 0;
FONT_COLOR = [1 1 1];
labels = [{0},{0}];
for n = 1:length(files)
    cMin = cMax + 1;
    cMax = cMax + c.colorLen;
    ch = colorbar('horiz');
    set(ch, 'tag', ['Colorbar', num2str(n)]);
    set(ch, 'position',pos);
    set(ch, 'xlim', [1, size(c.cmap, 1)]);
    pos(1) = pos(1) + pos(3) + offset;
    set(ch, 'xlim', [cMin, cMax]);
    set(ch, 'xtick', []);
    set(ch, 'color', FONT_COLOR);
    xlabel(char(c.modality{n}), 'parent', ch);
end

function [F,A,C,I] = ica_fuse_plot_FNC(FNC, CLIM, LABEL, RSN_I, F, T, sh, MOD, MOD_NAMES,excludeDiagonals)
% [F,A,C,I] = plot_FNC(FNC, CLIM, LABEL, RSN_I, F, T)

% define the size of each "MODULE"
%MOD = [4, 2, 8, 10, 14, 9];
if (~exist('MOD', 'var'))
    MOD = ones(1, length(RSN_I));
end
cMOD = cumsum(MOD);

if (~exist('MOD_NAMES', 'var'))
    MOD_NAMES = [];
end

if (~isempty(MOD_NAMES))
    tmpNames = repmat({''}, sum(MOD), 1);
    
    for nM = 1:length(MOD_NAMES)
        indx = ceil(MOD(nM)/2);
        if (nM > 1)
            indx = indx + sum(MOD(1:nM - 1));
        end
        tmpNames{indx} = MOD_NAMES{nM};
    end
    MOD_NAMES = tmpNames;
end

if (~exist('excludeDiagonals', 'var'))
    excludeDiagonals = 1;
end



if isvector(FNC)
    FNC = ica_fuse_vec2mat(FNC, 1);
else
    % make sure Nans are along the diagonal
    if (excludeDiagonals)
        FNC = ica_fuse_vec2mat(ica_fuse_mat2vec(FNC),1);
    end
end


if ~exist('CLIM', 'var') || isempty(CLIM)
    if any(FNC(:)<0)
        CLIM = [-max(abs(FNC(:))), max(abs(FNC(:)))];
    else
        temp = ica_fuse_mat2vec(FNC);
        CLIM = [min(temp(:)), max(temp(:))];
    end
end

if ~exist('T', 'var') || isempty(T)
    T = 'correlation';
end

if ~exist('F', 'var') || isempty(F)
    F = figure('Color', 'w');
end

set(F, 'invertHardcopy', 'off');

if (~exist('sh', 'var') || isempty(sh))
    sh = gca;
end

bgColor = get(F, 'color');
if (all(bgColor == 0))
    foregroundcolor = 'w';
else
    foregroundcolor = 'k';
end


%imagesc(FNC, CLIM);
ica_fuse_simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), FNC');
set(sh, 'clim', CLIM); axis ij
c = get(sh, 'Children');
set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', foregroundcolor);
set(c(length(RSN_I)-cMOD+1), 'Color', bgColor)
set(c(2*(length(RSN_I)+1)-cMOD), 'Color', bgColor)
uistack(c(2*(length(RSN_I)+1)-cMOD), 'top')
A = sh;
I = c(find(strcmp(get(c, 'Type'),'image')));
C = colorbar;
set(get(C, 'YLabel'), 'String', T)
axis square

try
    set(C, 'color', foregroundcolor);
    set(get(C, 'label'), 'color', foregroundcolor);
catch
end

%% labeling
if exist('LABEL', 'var') && ~isempty(LABEL)
    
    set(A, 'YTick', 1:length(RSN_I), 'YTickLabel', LABEL(RSN_I), 'TickDir', 'out', 'XTick', 1:length(RSN_I),'XTickLabel',[]);
    b = get(A, 'XTick');
    c=get(A,'YTick');
    ypos = c(end) +1.25;
    th=text(1:length(b),repmat(ypos,length(b),1), LABEL(RSN_I),'HorizontalAlignment','right','rotation',90, ...
        'color', foregroundcolor);
end

%% labeling
if exist('MOD_NAMES', 'var') && ~isempty(MOD_NAMES)
    
    %set(A, 'YTick', 1:length(RSN_I), 'YTickLabel', MOD_NAMES(RSN_I), 'TickDir', 'out', 'XTick', 1:length(RSN_I),'XTickLabel',[]);
    b = get(A, 'XTick');
    c=get(A,'YTick');
    ypos = c(end) + 4;
    xpos = b(end) + 4;
    th=text(repmat(-4,length(b),1), 1:length(b), MOD_NAMES(RSN_I),'HorizontalAlignment','right','rotation',90, ...
        'color', foregroundcolor, 'fontweight', 'bold');
    th=text(1:length(c),repmat(xpos,length(c),1), MOD_NAMES(RSN_I),'HorizontalAlignment','right','rotation',90, ...,
        'color', foregroundcolor, 'fontweight', 'bold');
end


if exist('RSN_I', 'var') && ~isempty(RSN_I)
    hold on;
    
    for ii = 1:length(RSN_I);
        T =text(ii,ii,LABEL{ii});
        set(T, 'Color', foregroundcolor, 'HorizontalAlignment', 'Center');
    end
end
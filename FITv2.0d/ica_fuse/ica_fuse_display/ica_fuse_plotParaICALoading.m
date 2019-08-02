function ica_fuse_plotParaICALoading(varargin)
% Plot Parallel ICA Loading parameters

ica_fuse_defaults;
global FIG_FG_COLOR;
global FIG_FG_COLOR;


titleColor = 'c';
line_width = 1.5;

% Loop over arguments
for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'data')
        dataStruct = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'parent')
        axesH = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'num_features')
        numFeatures = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'colors')
        colors = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'legend_string')
        legendString = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'titlecolor')
        titleColor = varargin{ii + 1};
    end
end
% End loop over arguments


% Initialise loading coefficient length
lengthLoadingCoeff = zeros(1, numFeatures);
for nF = 1:numFeatures
    temp = dataStruct.feature(nF).loadingCoeff;
    if (isnumeric(colors{nF}))
        pH = plot(temp, 'parent', axesH, 'linewidth', line_width);
        set(pH, 'color', colors{nF});
    else
        pH = plot(temp, colors{nF}, 'parent', axesH, 'linewidth', line_width);
    end
    lengthLoadingCoeff(nF) = length(temp);
    hold on;
end


hold off;
%ica_fuse_legend(legendString);

numSubjects = dataStruct.numSubjects;
groupNames = dataStruct.groupNames;

plotGroupsText = 1;
if length(find(lengthLoadingCoeff == sum(numSubjects))) ~= numFeatures
    plotGroupsText = 0;
end

if plotGroupsText
    if numSubjects > 1
        axisYlim = get(axesH, 'ylim');
        ySpacing = (axisYlim(2) - axisYlim(1))/4;
        plotYY = (axisYlim(1):ySpacing:axisYlim(2));
        % Loop over group names
        endTp = 0;
        for nS = 1:length(numSubjects)
            valX =  sum(numSubjects(1:nS));
            if nS ~= length(numSubjects)
                hold on;
                plotXX(1:length(plotYY)) = valX;
                plot(plotXX, plotYY, '--y', 'parent', axesH, 'linewidth', line_width);
            end
            xPos = (endTp + 0.025*numSubjects(nS)) / sum(numSubjects);
            yPos = 0.1;
            text(xPos, yPos, deblank(groupNames(nS, :)), 'parent', axesH, 'units', 'normalized');
            endTp = endTp + numSubjects(nS);
        end
        % End loop over group names
    end
    hold off;
end

%xlabel('Subjects', 'parent', axesH);
title(dataStruct.axesTitle, 'parent', axesH, 'color', titleColor);
xlabel('Subjects', 'parent', axesH);
ylabel('Normalized Units', 'parent', axesH);
axis(axesH, 'tight');
set(axesH, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
ica_fuse_legend(legendString);
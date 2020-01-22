function ica_fuse_histogram_plot2(fusionFile, bestComp, selGroups, selFeatures, threshold, histogram_criteria)
% Plot cross-task histograms of individual subjects of the selected groups
%
% Input:
% 1. fusionFile - full file path of fusion MAT file

try

    if isempty(which('ttest2.m'))
        error('Need statistics toolbox to plot histograms for subjects');
    end

    % Load defaults
    ica_fuse_defaults;

    global FUSION_INFO_MAT_FILE;
    global UI_FONT_NAME;
    global AX_COLOR;
    global COLORMAPFILE;
    global FIG_FG_COLOR;
    global Z_THRESHOLD;
    global HISTOGRAM_THRESHOLD;
    global NUM_BINS;

    fontSize = 0.01;
    % Select fusion file
    if ~exist('fusionFile', 'var')
        fusionFile = ica_fuse_selectEntry('title', 'Select fusion information file for plotting histograms', ...
            'typeEntity', 'file', 'typeSelection', 'Single', 'filter', ['*', FUSION_INFO_MAT_FILE, '*.mat']);
        drawnow;
    end

    load(fusionFile);

    if ~exist('fusionInfo', 'var')
        error(['Selected file: ', fusionFile, ' is not a valid fusion parameter file']);
    end

    if ~isfield(fusionInfo, 'run_analysis')
        error('Please run the analysis inorder to plot histograms of individual subjects');
    end

    [outputDir, fileN, extn] = fileparts(fusionFile);
    if isempty(outputDir)
        outputDir = pwd;
    end

    cd(outputDir);

    %%%%%%% Get the required parameters from run_analysis field %%%%%%%%%%%

    % Get the information from the fusion file
    numComp = fusionInfo.run_analysis.numComp;

    % Get the feature names
    featureNames = str2mat(fusionInfo.run_analysis.dataInfo(1).feature.name);

    % Get the groups names
    groupNames = str2mat(fusionInfo.run_analysis.dataInfo.name);

    % Get the ouput files
    outputFiles = fusionInfo.run_analysis.outputFiles;

    % Get the dataInfo
    dataInfo = fusionInfo.run_analysis.dataInfo;

    % Back reconstruct files
    backReconstructFiles = fusionInfo.run_analysis.backReconstructFiles;

    % Number of groups
    numGroups = fusionInfo.run_analysis.numGroups;

    % Number of features
    numFeatures = fusionInfo.run_analysis.numFeatures;

    % voxel indices
    maskIndices = fusionInfo.run_analysis.mask_ind;
    
    [maskIndices] = ica_fuse_form_maskInd(maskIndices, dataInfo);

    % Number of subjects
    numSubjects = fusionInfo.run_analysis.numSubjects;
    
    if length(fusionInfo.run_analysis.numSubjects) ~= length(dataInfo)
        numSubjects = fusionInfo.run_analysis.numSubjects;
        numSubjects = repmat(numSubjects, 1, length(dataInfo));
    else
        numSubjects = fusionInfo.run_analysis.numSubjects;
    end

    % Feature data length
    featureDataLength = fusionInfo.run_analysis.featureDataLength;
    dataLength = featureDataLength(1).Length;

    normalize = fusionInfo.run_analysis.normalize;

    % Normalization parameters
    featureNormPara = fusionInfo.run_analysis.featureNormPara;
    stdParameters = fusionInfo.run_analysis.stdParameters;


    % Fusion prefix
    fusionPrefix = fusionInfo.run_analysis.prefix;

    clear fusionInfo;

    listSize = [300 400];

    % Form component numbers string
    for nn = 1:numComp
        compStr(nn).name =  ['Component ', num2str(nn)];
    end


    % Selected groups
    if ~exist('selGroups', 'var')

        if numGroups > 2

            for nGroup = 1:numGroups
                answerString(nGroup).name = [deblank(groupNames(nGroup, :))];
            end

            % Select atmost 2 groups
            [selGroups, name_button] = ica_fuse_listdlg('PromptString','Select Groups', 'SelectionMode','multiple',...
                'ListString', str2mat(answerString.name), 'movegui', 'east', 'windowStyle', 'modal', 'maxselection', 2, 'title_fig', 'Select Groups');

            if name_button == 0
                error('Groups must be selected in order to do generate histograms');
            end

        elseif numGroups == 1
            selGroups = 1;
        else
            selGroups = [1 2];
        end

    end
    % End for selected groups


    % Selected Features
    if ~exist('selFeatures', 'var')

        % Select features
        if numFeatures > 2
            clear answerString;
            for nFea = 1:numFeatures
                answerString(nFea).name = ['Feature ', deblank(featureNames(nFea, :))];
            end

            % Select atmost 2 groups
            [selFeatures, name_button] = ica_fuse_listdlg('PromptString','Select Features', 'SelectionMode','multiple',...
                'ListString', str2mat(answerString.name), 'movegui', 'east', 'windowStyle', 'modal', 'maxselection', 2, 'title_fig', 'Select Features');

            if name_button == 0
                error('Features must be selected in order to generate cross task histograms');
            end

        elseif numFeatures == 1
            selFeatures = 1;
        else
            selFeatures = [1 2];
        end

    end
    % End for selected features

    if ~exist('bestComp', 'var')
        [bestComp, name_button] = ica_fuse_listdlg('PromptString','Select a component', 'SelectionMode', 'single',...
            'ListString', str2mat(compStr.name), 'movegui', 'east', 'windowStyle', 'modal', ...
            'title_fig', 'Select a component');

        if name_button == 0
            error('Component must be selected in order to generate cross task histograms');
        end

    end

    fprintf('\n');
    disp(['Using component ', num2str(bestComp), ' indices to plot the histograms of subjects ']);
    disp(['and the threshold applied will be first ', num2str(threshold), ' maximum values']);

    fprintf('\n');

    helpMsg = 'Calculating histograms. Please wait...';
    disp(helpMsg);
    helpHandle = helpdlg(helpMsg, helpMsg);

    load(COLORMAPFILE); % colormap file


    % Load back Reconstruction parameters
    backReconstructFile = fullfile(outputDir, backReconstructFiles(1).name);
    load(backReconstructFile);

    if ~exist('threshold', 'var')
        threshold = HISTOGRAM_THRESHOLD;
    end

    if ~exist('histogram_criteria', 'var')
        histogram_criteria = 'features';
    end


    fprintf('\n');

    [feature1StartInd, feature1EndInd] = ica_fuse_get_featureInd(selFeatures(1), dataLength);

    % Determine start and end indices for feature 2
    if length(selFeatures) > 1
        [feature2StartInd, feature2EndInd] = ica_fuse_get_featureInd(selFeatures(2), dataLength);
    end

    groupIndices = repmat(struct('ind', []), 1, length(selGroups));
    % Selected group indices
    for nGroup = 1:length(selGroups)
        groupIndices(nGroup).ind = ica_fuse_get_groupInd(selGroups(nGroup), numSubjects);
    end


    if strcmpi(histogram_criteria, 'features')
        % Plot histograms of the features

        % Load data and normalize
        % stack the data
        [stackInfo] = ica_fuse_stack_data('data_info', dataInfo, 'mask_ind', maskIndices, 'normalize_scheme', ...
            normalize, 'num_subjects', numSubjects);
        drawnow;
        dat = stackInfo.data; clear stackInfo;
        dataN = dat(1).data; clear dat;

        dataN = dataN([groupIndices.ind], :);

        %%%%%% Sort the component voxel values for their respective features
        group1_task1 = icasig(bestComp, feature1StartInd:feature1EndInd);
        %group1_task1 = group1_task1(mask_ind);

        [g1, xIndex] = sort(abs(group1_task1(:)));
        xIndex = xIndex(end:-1:1);


        % Calculate the default range
        temp1 = dataN(:, feature1StartInd:feature1EndInd);
        temp1 = temp1(:, xIndex);
        %temp1 = temp1(:, mask_ind);

        rng1 = [ica_fuse_minN(temp1), ica_fuse_maxN(temp1)];
        clear temp1;

        if length(selFeatures) > 1

            group1_task2 = icasig(bestComp, feature2StartInd:feature2EndInd);
            %group1_task2 = group1_task2(mask_ind);
            [g2, yIndex] = sort(abs(group1_task2(:)));
            yIndex = yIndex(end:-1:1);

            if length(xIndex) > length(yIndex)
                xIndex = xIndex(1:length(yIndex));
            elseif length(yIndex) > length(xIndex)
                yIndex = yIndex(1:length(xIndex));
            end

            if threshold < length(xIndex)
                xIndex = xIndex(1:threshold);
                yIndex = yIndex(1:threshold);
            end

            temp2 = dataN(:, feature2StartInd:feature2EndInd);
            temp2 = temp2(:, yIndex);
            %temp2 = temp2(:, mask_ind);

            rng2 = [ica_fuse_minN(temp2), ica_fuse_maxN(temp2)];
            clear temp2;

            temp1 = dataN(:, feature1StartInd:feature1EndInd);
            temp1 = temp1(:, xIndex);
            rng1 = [ica_fuse_minN(temp1), ica_fuse_maxN(temp1)];
            clear temp1;

        else
            
            if threshold < length(xIndex)
                xIndex = xIndex(1:threshold);
            end
            
             temp1 = dataN(:, feature1StartInd:feature1EndInd);
            temp1 = temp1(:, xIndex);
            rng1 = [ica_fuse_minN(temp1), ica_fuse_maxN(temp1)];
            clear temp1;

            yIndex = xIndex;
            rng2 = rng1;

        end

        clear group1_task1;
        clear group1_task2;
        %%%%%%%%%% End for determining the voxel values %%%%%%%%%


        % Bin size
        xsize = (max(rng1) - min(rng1)) / NUM_BINS;
        ysize = (max(rng2) - min(rng2)) / NUM_BINS;

        % Generate indices
        %xind = (min(rng1) + xsize):xsize:max(rng1);
        %yind = (min(rng2) + ysize):ysize:max(rng2);

        xind = linspace(min(rng1), max(rng1), NUM_BINS);
        yind = linspace(min(rng2), max(rng2), NUM_BINS);


        selSubjects = numSubjects(selGroups);

        % Loop over number of subjects
        %for nSub = 1:length(selSubjects)
        % Loop over selected groups
        for nGroup = 1:length(selGroups)
            %groupInd = 1:sum(groupIndices(1:nGroup).ind);

            [groupInd] = ica_fuse_get_groupInd(nGroup, selSubjects);
            countN = 0;
            for nSub = groupInd

                countN = countN + 1;

                % task1 data
                %task1 = dataN(selFeatures(1)).data(groupInd(nSub), :);
                task1 = dataN(nSub, feature1StartInd:feature1EndInd);
                task1 = task1(xIndex);
                %task1 = task1(mask_ind);
                task1 = task1(:)';

                if length(selFeatures) > 1
                    %task2 = dataN(selFeatures(2)).data(groupInd(nSub), :);
                    task2 = dataN(nSub, feature2StartInd:feature2EndInd);
                    task2 = task2(yIndex);
                    %task2 = task2(mask_ind);
                    task2 = task2(:)';
                else
                    task2 = task1;
                end

                disp(['Calculating cross task histograms for group ', deblank(groupNames(selGroups(nGroup), :)), ' subject ', ...
                    num2str(countN)]);

                % generate histogram
                hData{1} = task1;
                hData{2} = task2;
                hIndex{1} = xind;
                hIndex{2} = yind;
                [out] = ica_fuse_histnd(hData, hIndex);
                clear hIndex; clear hData;

                graphData.group(nGroup).data(:, :, countN) = log(out + 1);
                clear out;
            end
        end
        %end

        clear dataN;


        %%%%%%%% Display histograms %%%%%%%
        % Plot Histogram of groups
        for nn = 1:length(graphData.group)
            [graphicsH(nn).H] = plotIm(graphData.group(nn).data, ...
                ['Group ', deblank(groupNames(selGroups(nn), :))], hot);
        end

        % Plot difference histograms
        % Plot First group followed by second group
        m1 = mean(graphData.group(1).data, 3);

        % Mean for group 1
        graphCount = length(graphicsH);
        axesTitle = ['Mean of Group ', deblank(groupNames(selGroups(1), :))];
        [graphicsH(graphCount + 1).H] = plotIm(m1, axesTitle, hot, axesTitle);

        graphCount = length(graphicsH);
        drawArrows(graphicsH(graphCount).H, featureNames, selFeatures);

        if length(graphData.group) == 2

            m2 = mean(graphData.group(2).data, 3);

            graphCount = length(graphicsH);
            axesTitle = ['Mean of Group ', deblank(groupNames(selGroups(2), :))];
            % Mean for group 2
            [graphicsH(graphCount + 1).H] = plotIm(m2, axesTitle, hot, axesTitle);
            graphCount = length(graphicsH);
            drawArrows(graphicsH(graphCount).H, featureNames, selFeatures);

            out = m1 - m2;
            group1Name = deblank(groupNames(selGroups(1), :));
            group2Name = deblank(groupNames(selGroups(2), :));

            graphCount = length(graphicsH);
            axesTitle = ['Group ', group1Name , ' - Group ', group2Name];
            [graphicsH(graphCount + 1).H] = plotIm(out, ['Group ', group1Name , ' - Group ', group2Name], coldhot, axesTitle);

            % Number of graphs
            graphCount = length(graphicsH);

            axesHandle = get(graphicsH(graphCount).H, 'currentAxes');
            axesPos = get(axesHandle, 'position');

            drawArrows(graphicsH(graphCount).H, featureNames, selFeatures);

            set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
            colorbarPos = [axesPos(1) + axesPos(3) + 0.03, axesPos(2), 0.05, axesPos(4)];
            colorbarHandle = colorbar('peer', axesHandle);
            set(colorbarHandle, 'position', colorbarPos);
            % Maximum
            xPos = 0.1; yPos = 1.05;
            text(xPos, yPos, group1Name, 'units', 'normalized', 'parent', colorbarHandle, 'color', FIG_FG_COLOR,  ...
                'fontsize', fontSize, 'HorizontalAlignment', 'left', ...
                'FontName', UI_FONT_NAME, 'fontunits', 'normalized');

            % Minimum
            yPos = -0.05;
            text(xPos, yPos, group2Name, 'units', 'normalized', 'parent', colorbarHandle, 'color', FIG_FG_COLOR,  ...
                'fontsize', fontSize, 'HorizontalAlignment', 'left', 'FontName', UI_FONT_NAME, 'fontunits', ...
                'normalized');
            set(colorbarHandle, 'YTickLabel', []);
            set(colorbarHandle, 'YTick', []);

        end

        delete(helpHandle);

        drawnow;

        % Plot Next previous and exit buttons
        ica_fuse_plotNextPreviousExitButtons(graphicsH);

        disp('Done calculating histograms');

    else

        data = A*icasig;
        % Initialise groups icasig variable
        groups_icasig = zeros(numGroups, size(icasig, 2));

        % Loop over groups
        for nGroups = 1:numGroups
            groupInd = ((nGroups - 1)*numSubjects + 1 : nGroups*numSubjects);
            % get the back reconstruted component
            groups_icasig(nGroups, :) = pinv(A(groupInd, bestComp))*data(groupInd, :);
        end

        %%%%%%%%%%%%%% End for groups back reconstruction %%%%%

        index1 = [feature1StartInd:feature1EndInd];
        if length(selFeatures) > 1
            index2 = [feature2StartInd:feature2EndInd];
        else
            index2 = [];
        end

        [histData(1).task1, histData(1).task2, histData(1).rng1, histData(1).rng2, ...
            histData(1).xIndex, histData(1).yIndex] = getDataRange(groups_icasig, selGroups(1), ...
            selFeatures, index1, index2);

        if length(selGroups) > 1
            [histData(2).task1, histData(2).task2, histData(2).rng1, histData(2).rng2, ...
                histData(2).xIndex, histData(2).yIndex] = getDataRange(groups_icasig, selGroups(2), ...
                selFeatures, index1, index2);
        end


        for nn = 1:length(histData)
            % Bin size
            xsize = (max(histData(nn).rng1) - min(histData(nn).rng1)) / NUM_BINS;
            ysize = (max(histData(nn).rng2) - min(histData(nn).rng2)) / NUM_BINS;

            % Generate indices
            xind = (min(histData(nn).rng1) + xsize):xsize:max(histData(nn).rng1);
            yind = (min(histData(nn).rng2) + ysize):ysize:max(histData(nn).rng2);

            disp(['Calculating cross task histograms for group ', deblank(groupNames(selGroups(nn), :))]);

            hData{1} = histData(nn).task1;
            hData{2} = histData(nn).task2;
            hIndex{1} = xind;
            hIndex{2} = yind;
            % generate histogram
            [out] = ica_fuse_histnd(hData, hIndex);
            clear hIndex; clear hData;

            graphData(nn).data = log(out + 1);
            graphData(nn).groupName = deblank(groupNames(selGroups(nn), :));

        end

        %%%%%%%% Display histograms %%%%%%%
        % Plot Histogram of groups
        for nn = 1:length(graphData)
            axesTitle = ['Group ', graphData(nn).groupName];
            [graphicsH(nn).H] = plotIm(graphData(nn).data, axesTitle, hot, axesTitle);
            axesHandle = get(graphicsH(nn).H, 'currentAxes');
            axesPos = get(axesHandle, 'position');
            drawArrows(graphicsH(nn).H, featureNames, selFeatures);
            set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);

        end

        if length(graphData) == 2
            out = graphData(1).data - graphData(2).data;
            group1Name = graphData(1).groupName;
            group2Name = graphData(2).groupName;

            graphCount = length(graphicsH);
            axesTitle = ['Group ', group1Name , ' - Group ', group2Name];
            [graphicsH(graphCount + 1).H] = plotIm(out, ['Group ', group1Name , ' - Group ', group2Name], ...
                coldhot, axesTitle);

            % Number of graphs
            graphCount = length(graphicsH);

            axesHandle = get(graphicsH(graphCount).H, 'currentAxes');
            axesPos = get(axesHandle, 'position');

            drawArrows(graphicsH(graphCount).H, featureNames, selFeatures);

            set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
            colorbarPos = [axesPos(1) + axesPos(3) + 0.03, axesPos(2), 0.05, axesPos(4)];
            colorbarHandle = colorbar('peer', axesHandle);
            set(colorbarHandle, 'position', colorbarPos);
            % Maximum
            xPos = 0.1; yPos = 1.05;
            text(xPos, yPos, group1Name, 'units', 'normalized', 'parent', colorbarHandle, 'color', FIG_FG_COLOR,  ...
                'fontsize', fontSize, 'HorizontalAlignment', 'left', ...
                'FontName', UI_FONT_NAME, 'fontunits', 'normalized');

            % Minimum
            yPos = -0.05;
            text(xPos, yPos, group2Name, 'units', 'normalized', 'parent', colorbarHandle, 'color', FIG_FG_COLOR,  ...
                'fontsize', fontSize, 'HorizontalAlignment', 'left', 'FontName', UI_FONT_NAME, 'fontunits', ...
                'normalized');
            set(colorbarHandle, 'YTickLabel', []);
            set(colorbarHandle, 'YTick', []);
        end

        delete(helpHandle);

        drawnow;

        % Plot Next previous and exit buttons
        ica_fuse_plotNextPreviousExitButtons(graphicsH);

    end


catch

    if exist('helpHandle', 'var')
        if ishandle(helpHandle)
            delete(helpHandle);
        end
    end

    %rethrow(lasterror);
    ica_fuse_displayErrorMsg;

end


function [figHandle] = plotIm(data, titleStr, cmap, axesTitle)
% Sub function to plot images in a sub plot

ica_fuse_defaults;
global AX_COLOR; % axes color
global LEGEND_COLOR; % LEGEND_COLOR for Matlab 7.0
global FIG_FG_COLOR;
global UI_FONT_NAME;
global UI_FONT_UNITS;
global UI_FONT_SIZE;

figHandle = ica_fuse_getGraphics(titleStr, 'normal', titleStr, 'on');

if length(size(data)) == 3
    numPlots = size(data, 3);
else
    numPlots = 1;
end

userdata = zeros(size(data));

for nn = 1:numPlots
    userdata(:, :, nn) = (squeeze(data(:, :, nn)))';
end

% Display GUI Options Menu
figOptions = uimenu('parent', figHandle, 'tag', 'extract_data', 'label', 'Extract Data', 'callback', ...
    {@extractDataCallback, figHandle}, 'userdata', userdata);

set(figHandle, 'colormap', cmap);

numRows = ceil(sqrt(numPlots));
numCols = ceil(numPlots/numRows);

fontSize = (0.05);

if numRows*numCols > 1
    fontSize = 0.1;
end

countN = 0;
for nRow = 1:numRows

    for nCol = 1:numCols
        countN = countN + 1;

        if countN <= numPlots
            if numPlots > 1
                figure(figHandle);
                axesHandle = subplot(numRows, numCols, countN);
            else
                axesHandle = axes('parent', figHandle, 'units', 'normalized', 'position', [0.1 0.1 0.7 0.7]);
            end

            set(axesHandle, 'units', 'normalized', 'color', AX_COLOR, 'FontName', UI_FONT_NAME, 'fontunits', ...
                'normalized', 'fontsize', fontSize);
            set(axesHandle, 'YColor', FIG_FG_COLOR, 'XColor', FIG_FG_COLOR);
            set(axesHandle, 'YTickLabel', []);
            set(axesHandle, 'YTick', []);
            set(axesHandle, 'XTickLabel', []);
            set(axesHandle, 'XTick', []);
            subject_str = ['Sub ', num2str(countN)];
            temp = squeeze(data(:, :, countN));
            temp = temp';
            ImageAxis = image(temp, 'parent', axesHandle, 'CDataMapping', 'scaled');

            set(axesHandle, 'clim', [min(data(:)), max(data(:))]);
            clear temp;
            axis(axesHandle, 'off');
            if exist('axesTitle', 'var')
                title(axesTitle, 'parent', axesHandle);
            else
                title(['Sub ', num2str(countN)], 'parent', axesHandle);
            end
        else
            break;
        end

    end

end



function [group1_task1, group1_task2, rng1, rng2, xIndex, yIndex] = getDataRange(groups_icasig, groupNum, selFeatures, ...
    index1, index2);


%%%%%% Sort the component voxel values for their respective features
group1_task1 = groups_icasig(groupNum, index1);
%group1_task1 = group1_task1(mask_ind);

[g1, xIndex] = sort(abs(group1_task1(:)));
xIndex = xIndex(end:-1:1);

if length(selFeatures) > 1
    %%%%%% Sort the component voxel values for their respective features
    group1_task2 = groups_icasig(groupNum, index2);
    %group1_task1 = group1_task1(mask_ind);

    [g2, yIndex] = sort(abs(group1_task2(:)));
    yIndex = yIndex(end:-1:1);

    if length(xIndex) > length(yIndex)
        xIndex = xIndex(1:length(yIndex));
    elseif length(yIndex) > length(xIndex)
        yIndex = yIndex(1:length(xIndex));
    end

    [xIndex] = sort(xIndex);
    [yIndex] = sort(yIndex);

    group1_task1 = group1_task1(xIndex);
    group1_task2 = group1_task2(yIndex);

    rng1 = [ica_fuse_minN(group1_task1), ica_fuse_maxN(group1_task1)];
    rng2 = [ica_fuse_minN(group1_task2), ica_fuse_maxN(group1_task2)];

else
    [xIndex] = sort(xIndex);
    group1_task1 = group1_task1(xIndex);
    rng1 = [ica_fuse_minN(group1_task1), ica_fuse_maxN(group1_task1)];
    group1_task2 = group1_task1;
    yIndex = xIndex;
    rng2 = rng1;
end


function drawArrows(graphicsHandle, featureNames, selFeatures)
% Draw arrows

ica_fuse_defaults;
global FIG_FG_COLOR;
global UI_FONT_NAME;

axesHandle = get(graphicsHandle, 'currentAxes');
axesPos = get(axesHandle, 'position');
xPos = 0.45; yPos = -0.02;

if length(selFeatures) == 1
    selFeatures = repmat(selFeatures, 1, 2);
end

fontSize = 0.01;

% Feature 1
text(xPos, yPos, deblank(featureNames(selFeatures(1), :)), 'units', 'normalized', 'parent', ...
    axesHandle, 'color', FIG_FG_COLOR, 'fontsize', fontSize, 'HorizontalAlignment', ...
    'center', 'FontName', UI_FONT_NAME, 'fontunits', 'normalized');

%% Draw a horizontal arrow
try

    % Annotation
    annoX = [axesPos(1), axesPos(1) + axesPos(3)];

    annoY = [axesPos(2) + yPos - 0.03, axesPos(2) + yPos - 0.03];

    % annotation object
    annoH = annotation('arrow', annoX, annoY, 'color', FIG_FG_COLOR);

catch

end


xPos = 0.05; yPos = 0.45;
% Feature 2
text(xPos, yPos, deblank(featureNames(selFeatures(2), :)), 'units', 'normalized', 'parent', ...
    axesHandle, 'color', FIG_FG_COLOR, 'fontsize', fontSize, 'HorizontalAlignment', ...
    'center', 'FontName', UI_FONT_NAME, 'fontunits', 'normalized', ...
    'rotation', 90);

%% Draw a vertical arrow
try

    % Annotation
    annoX = [axesPos(1) - 0.03, axesPos(1) - 0.03];

    annoY = [axesPos(2), axesPos(2) + axesPos(4)];

    % annotation object
    annoH = annotation('arrow', annoX, annoY, 'color', FIG_FG_COLOR);

catch

end


function extractDataCallback(hObject, event_data, handles) 
% Extract histogram data

data = get(hObject, 'userdata');

% open input dialog box
prompt = {'Enter a valid file name'};
dlg_title = 'Save Figure Data As';
num_lines = 1;
def = {'histogram_data'};
% save the file with the file name specified
fileName = ica_fuse_inputdlg2(prompt, dlg_title, num_lines, def);

if ~isempty(fileName)
    [dd, fileN, extn] = fileparts(fileName{1});
    if ~isempty(dd)
        outputDir = dd;
    else
        outputDir = pwd;
    end
    clear fileName;
    fileName = [fileN, '.mat'];
    fileName = fullfile(outputDir, fileName);
    ica_fuse_save(fileName, 'data');
    disp(['Figure data is saved in file: ', fileName]);
end

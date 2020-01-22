function ica_fuse_stats_neural_net_fusion(multi_param_file, opts)

ica_fuse_defaults;
global IMAGE_VALUES;
global Z_THRESHOLD;
global ANATOMICAL_PLANE;
global CONVERT_TO_Z;
global ANATOMICAL_FILE;
global SCREENSIZE;


if (~exist('multi_param_file', 'var'))
    multi_param_file = ica_fuse_selectEntry('typeEntity', 'file', 'filter', '*neural*net*fusion*mat', 'title', 'Select SBM ICA Parameter file/Neural net fusion file');
    drawnow;
end

if (isempty(multi_param_file))
    error('Select neural net fusion file/sbm ica parameter file');
end

load(multi_param_file);


if (~exist('neuralNetFusionInfo', 'var'))
    error('Selected file is not a valid neural network fusion file');
end


thresh = 0.05;
try
    thresh = opts.thresh;
catch
    
end

thresh_criteria = 'none';
try
    thresh_criteria = opts.thresh_criteria;
catch
end



try
    imageValues = opts.image.image_values;
catch
    imageValues = IMAGE_VALUES;
end

try
    threshold_maps = abs(opts.image.z_threshold);
catch
    threshold_maps = Z_THRESHOLD;
end

try
    anatomical_plane = lower(opts.anatomical_plane);
catch
    anatomical_plane = ANATOMICAL_PLANE;
end

try
    slices_in_mm = opts.image.slices_in_mm;
catch
    slices_in_mm = -40:4:72;
end

try
    convertToZ = opts.image.convert_to_z;
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

opts.image.disp_opts=disp_opts;
opts.image.image_values = imageValueStr;
opts.image.convertToZ = convertToZStr;
opts.image.slices_in_mm = slices_in_mm;
opts.image.threshold_maps=Z_THRESHOLD;
opts.thresh_criteria = thresh_criteria;
opts.thresh = thresh;


outputDir = fileparts(multi_param_file);
if (isempty(outputDir))
    outputDir = pwd;
end

fname = fullfile(outputDir, neuralNetFusionInfo.alignment_file);
load(fname);
num_states = size(states_all, 1);
num_comps = size(states_all, 2);
num_subjects = size(subject_alignments, 1);

stateLabels = cellstr([repmat('State ', num_states, 1), num2str((1:num_states)')]);

subjectLabels = cellstr([repmat('Subject ', num_subjects, 1), num2str((1:num_subjects)')]);


R = [1, 1];
screenSize = SCREENSIZE;
SZ = min([screenSize(3), screenSize(4)]);
extendLeft = round(SZ*R(1)*.85);
extendUp = round(SZ*R(2)*.85);
x0= (screenSize(3)/2)-(extendLeft/2);
y0 = (screenSize(4)/2)-(extendUp/2);
Rect = [x0 y0 extendLeft extendUp];

if (~exist('opts', 'var') || ~isfield(opts, 'groupsInfo') || isempty(opts.groupsInfo))
    
    [groupName1, groupVal1] = ica_fuse_select_groups_gui(subjectLabels, 'Group 1 Info','groups_gui');
    if (isempty(groupVal1))
        error('Figure window was quit');
    end
    groupsInfo(1).name = groupName1;
    groupsInfo(1).value = groupVal1;
    [groupName2, groupVal2] = ica_fuse_select_groups_gui(subjectLabels, 'Group 2 Info','groups_gui');
    if (isempty(groupVal2))
        error('Figure window was quit');
    end
    groupsInfo(2).name = groupName2;
    groupsInfo(2).value = groupVal2;
    %group_cov_file = ica_fuse_selectEntry('typeEntity', 'file', 'filter', '*txt', 'title', 'Select groups covariate file ...');
    drawnow;
    opts.groupsInfo = groupsInfo;
    
    opts = getDefs(opts);
    imageValues = opts.image.image_values;
    convertToZ = opts.image.convertToZ;
    slices_in_mm = opts.image.slices_in_mm;
    threshold_maps = opts.image.threshold_maps;
    thresh_criteria = opts.thresh_criteria;
    thresh = opts.thresh;
    
    
    ica_fuse_neural_net_fusion_report(multi_param_file, opts);
    return;
    
else
    
    groupsInfo = opts.groupsInfo;
    
end

subject_alignments(subject_alignments == 0) = NaN;

[pval, tval]= ica_fuse_ttest(subject_alignments);

pval = squeeze(pval);
tval = squeeze(tval);

if (strcmpi(thresh_criteria, 'fdr'))
    p_fdr = ica_fuse_fdr(pval, thresh);
else
    p_fdr = thresh;
end

tval(pval > p_fdr) = NaN;
pval(pval > p_fdr) = NaN;

good_inds = find(isnan(pval) == 0);


one_sample_ttest.pval = pval;
one_sample_ttest.tval = tval;

if (~isempty(good_inds))
    
    count = 0;
    for nState = 1:size(pval, 1)
        for nComp = 1:size(pval, 2)
            if ~isnan(pval(nState, nComp))
                
                re_values = squeeze(subject_alignments(:, nState, nComp));
                
                tmpx = re_values(groupsInfo(1).value);
                tmpy = re_values(groupsInfo(2).value);
                
                tmpx = tmpx (tmpx > 0);
                tmpy = tmpy(tmpy > 0);
                
                if (length(tmpx) > 1 && length(tmpy) > 1)
                    count = count + 1;
                    [p_ttest2, t_ttest2] = icatb_ttest2(tmpx ,tmpy);
                    
                    [h, p_kstest,k_stat] = kstest2(tmpx, tmpy, 'Alpha', thresh);
                    groupStats(count).k_stat = k_stat;
                    groupStats(count).p_kstest = p_kstest;
                    
                    groupStats(count).value = re_values;
                    groupStats(count).pval = p_ttest2;
                    groupStats(count).tval = t_ttest2;
                    groupStats(count).state_no = nState;
                    groupStats(count).comp_no = nComp;
                end
                
            end
        end
    end
end

drawnow;


if (exist('groupStats', 'var'))
    
    if (strcmpi(thresh_criteria, 'fdr'))
        p_fdr_ttest2 = ica_fuse_fdr([groupStats.pval], thresh);
    else
        p_fdr_ttest2 = thresh;
    end
    
    groupStats([groupStats.pval] > p_fdr_ttest2) = [];
    
    load(neuralNetFusionInfo.smri_param_file);
    smri_comp_files = ica_fuse_rename_4d_file(ica_fuse_fullFile('directory', fileparts(neuralNetFusionInfo.smri_param_file), 'files', sesInfo.icaOutputFiles(1).ses(1).name));
    
    load(fullfile(outputDir, [neuralNetFusionInfo.prefix, '_neural_net_fusion_dfnc_centroids.mat']));
    load(neuralNetFusionInfo.dfnc_param_file);
    network_values = zeros(1, length(dfncInfo.userInput.comp));
    for nV = 1:length(network_values)
        network_values(nV) = length(dfncInfo.userInput.comp(nV).value);
    end
    network_names =  cellstr(char(dfncInfo.userInput.comp.name));
    
    if (length(network_names) == 1)
        network_names = '';
    end
    
    dfnc_comps = dfncInfo.comps(:);
    
    
    CLIM = max(abs(C_all(:)));
    CLIM = [-CLIM, CLIM];
    
    
    
    %% dFNC centroids and linked sMRI components are shown. Two sample t-statistic is also shown between the groups.
    for nState = 1:num_states
        
        inds = find([groupStats.state_no] == nState);
        
        if (~isempty(inds))
            
            H = figure('position', Rect);
            
            FNC = ica_fuse_vec2mat(C_all(nState, :));
            ica_fuse_plot_FNC(FNC, CLIM, cellstr(num2str(dfnc_comps)), (1:length(dfnc_comps)), H, ['State # ', num2str(nState), ' corr(z)'], ...
                gca, network_values, network_names);
            colormap(jet);
            
            title(['State # ', num2str(nState)], 'parent', gca);
            
            
            % Plot States FNC
            
            for nI = 1:length(inds)
                
                ica_fuse_image_viewer(deblank(smri_comp_files(groupStats(nI).comp_no, :)), 'structfile', ANATOMICAL_FILE, ...
                    'anatomical_view', anatomical_plane, 'slices_in_mm', slices_in_mm, 'image_values', imageValues, 'threshold', threshold_maps, 'convert_to_zscores', convertToZ,...
                    'labels', ['State # ', num2str(nState), ' vs Comp No: ', num2str(groupStats(nI).comp_no), ' T-stat: ', num2str(groupStats(nI).tval,'%0.3f')]);
                
            end
            
        end
        
    end
    
    drawnow;
    
    %% Density Plots: Density plot is computed for each group and Two-sample Kolmogorov-Smirnov test statistic is also reported.
    %
    for n = 1:length(groupStats)
        
        re_values = groupStats(n).value;
        tmpx = re_values(groupsInfo(1).value);
        tmpy = re_values(groupsInfo(2).value);
        
        tmpx = tmpx (tmpx > 0);
        tmpy = tmpy(tmpy > 0);
        [tmp1, xi] = ksdensity(re_values(re_values > 0));
        [f1, xx] = ksdensity(tmpx, xi);
        [f2, xx] = ksdensity(tmpy, xi);
        
        figure;
        color_val = [0.5, 0, 0.5];
        fill(xi, f1, color_val);
        hold on;
        color_val= [0, 0.5, 0.5];
        fill(xi, f2, color_val);
        title(['State # ', num2str(groupStats(n).state_no), ' comp # ', num2str(groupStats(n).comp_no), ' K-stat and p-value: k = ', num2str(groupStats(n).k_stat,'%0.3f'), ...
            '  p = ', num2str(groupStats(n).p_kstest,'%0.3f')]);
        legend(groupsInfo(1).name, groupsInfo(2).name, 'location', 'best');
        axis tight;
        
    end
    
    
    
end




drawnow;

%% Alignment scores between dFNC states and sMRI components in stacked bar plots.Alignment scores between dFNC states and sMRI components shown as an image.
H = figure('color', 'w', 'position', Rect, 'tag', 'alignment_scores');
subplot(2, 1, 1);
ah = get(H, 'currentAxes');
bh = bar(states_all,'stacked');
ylabel('sMRI components');
xlabel('dFNC states');
set(gca, 'XtickLabel', stateLabels);
title('Alignment scores (dFNC states vs sMRI)');
axis tight;


subplot(2, 1, 2);
imagesc(states_all);
set(gca,'Ytick',(1:num_states));
xlabel('sMRI components');
ylabel('dFNC states');
set(gca, 'YtickLabel', stateLabels);
title('Alignment scores (dFNC states vs sMRI)');
colormap(hot);
colorbar;



pause(0.3);

drawnow;

% Save files
fname = fullfile(outputDir, [neuralNetFusionInfo.prefix, '_stats.mat']);
statsInfo.one_sample_ttest = one_sample_ttest;
if (exist('groupStats', 'var'))
    statsInfo.groupStats = groupStats;
end

save(fname, 'statsInfo', 'opts');


function opts = getDefs(opts)

numParameters = 1;

inputText(numParameters).promptString = 'Enter p-value threshold';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(opts.thresh);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'thresh';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;
threshCriteriaOpts = {'FDR', 'None'};
chkImgInd = strmatch(lower(opts.thresh_criteria), lower(threshCriteriaOpts), 'exact');
inputText(numParameters).promptString = 'Select threshold criteria (t-tests)';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = threshCriteriaOpts;
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'thresh_criteria';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;


numParameters = numParameters + 1;

chkImgInd = strmatch(lower(opts.image.image_values), lower(opts.image.disp_opts), 'exact');
inputText(numParameters).promptString = 'Select image values (for display)';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = opts.image.disp_opts;
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'image_values';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = chkImgInd;


numParameters = numParameters + 1;
z_opts = {'Yes', 'No'};
chkZInd = strmatch(lower(opts.image.convertToZ), lower(z_opts), 'exact');
inputText(numParameters).promptString = 'Do you want to scale image values to z-scores (for display)?';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = z_opts;
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'convert_to_z';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = chkZInd;

numParameters = numParameters + 1;
inputText(numParameters).promptString = 'Enter threshold for maps (for display)';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(opts.image.threshold_maps);
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'thresh_maps';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;


answer = ica_fuse_inputDialog('inputtext', inputText, 'Title', 'Select nnet statistics params', 'handle_visibility',  'on', 'windowStyle', 'modal');

drawnow;

if (~isempty(answer))
    opts.thresh = answer{1};
    opts.thresh_criteria = answer{2};
    opts.image.image_values = answer{3};
    opts.image.convertToZ = answer{4};
    opts.image.threshold_maps = answer{5};
else
    error('Figure window was quit');
    
end
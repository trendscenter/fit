function fusionInfo = ica_fuse_copyfields_to_run(fusionInfo)
%% Copy fields from fusionInfo.setup_analysis to fusionInfo.run_analysis
%

%% Fields required
fieldsRequired = {'prefix', 'numGroups', 'numFeatures', 'dataInfo', 'normalize', 'numComp', 'algorithm'};

%% Add the fields required from setup_analysis to run_analysis
for nF = 1:length(fieldsRequired)
    currentField = fieldsRequired{nF};
    temp = getfield(fusionInfo.setup_analysis, currentField);
    fusionInfo.run_analysis = setfield(fusionInfo.run_analysis, currentField, temp);
    clear temp;
end

%% Mask indices
mask_ind = fusionInfo.setup_analysis.mask_ind;
mask_ind = ica_fuse_form_maskInd(mask_ind, fusionInfo.setup_analysis.dataInfo);
fusionInfo.setup_analysis.mask_ind = mask_ind;
fusionInfo.run_analysis.mask_ind = mask_ind;

%% Scaling to z-scores
if ~isfield(fusionInfo.setup_analysis, 'z_scores')
    fusionInfo.setup_analysis.z_scores = 0;
end

% Set the z_scores field in run_analysis
fusionInfo.run_analysis.z_scores = fusionInfo.setup_analysis.z_scores;

%%  Number of subjects per group
if ~isfield(fusionInfo.setup_analysis, 'numSubjects')
    numSubjects = ica_fuse_countSubjects(fusionInfo.setup_analysis.dataInfo);
else
    if length(fusionInfo.setup_analysis.numSubjects) ~= length(fusionInfo.setup_analysis.dataInfo)
        numSubjects = fusionInfo.setup_analysis.numSubjects;
        numSubjects = repmat(numSubjects, 1, length(fusionInfo.setup_analysis.dataInfo));
    else
        numSubjects = fusionInfo.setup_analysis.numSubjects;
    end
end

if (fusionInfo.run_analysis.numComp > sum(numSubjects))
    error(['Number of components (', num2str(fusionInfo.run_analysis.numComp), ...
        ') is greater than the number of data-sets (', num2str(sum(numSubjects)), ')']);
end


fusionInfo.run_analysis.numSubjects = numSubjects;

%% Type of PCA
if ~isfield(fusionInfo.setup_analysis, 'type_pca')
    fusionInfo.setup_analysis.type_pca = 'standard';
end

fusionInfo.run_analysis.type_pca = fusionInfo.setup_analysis.type_pca;

%% Check type of PCA
% Print to command window
if ~strcmpi(fusionInfo.setup_analysis.type_pca, 'standard')
    % Initialise reference vector
    reference = ones(sum(numSubjects), 1);
    if isfield(fusionInfo.setup_analysis, 'reference') && ~isempty(fusionInfo.setup_analysis.reference)
        reference = fusionInfo.setup_analysis.reference;
    else

        if (length(numSubjects) == 2)
            reference(1:numSubjects(1)) = -1;
        else
            for nS = 1:length(numSubjects)
                groupInd = ica_fuse_get_groupInd(nS, numSubjects);
                reference(groupInd) = 1/numSubjects(nS);
            end
        end
    end
else
    reference = [];
end

%% ICA Options
ICA_Options = {};
if isfield(fusionInfo.setup_analysis, 'ICA_Options')
    ICA_Options = fusionInfo.setup_analysis.ICA_Options;
end
fusionInfo.run_analysis.ICA_Options = ICA_Options;

% Add reference field to fusionInfo
fusionInfo.setup_analysis.reference = reference;
fusionInfo.run_analysis.reference = reference;

%% DIMS and voxels
if (~isfield(fusionInfo.setup_analysis, 'dims'))
    [dims, voxels] = ica_fuse_getFeatureDIM(mask_ind);
else
    dims = fusionInfo.setup_analysis.dims;
    voxels = fusionInfo.setup_analysis.voxels;
end

fusionInfo.run_analysis.dims = dims;
fusionInfo.run_analysis.voxels = voxels;


function sortResults = ica_fuse_sort_components(sortParameters, comb_number, saveResults)
%% Sort the components based on the selected sorting criteria
%
% Inputs:
% 1. sortParameters - structure containing the necessary fields for sorting
% the components
% 2. comb_number - MAT file number
% 3. saveResults - Options are 0 and 1
%
% Outputs:
% sortResults - Sorting results data structure

try
    
    if ~exist('comb_number', 'var')
        comb_number = 1;
    end
    
    if ~exist('saveResults', 'var')
        saveResults = 1;
    end
    
    % Get the necessary parameters
    outputDir = sortParameters.outputDir;
    dataInfo = sortParameters.dataInfo;
    output_prefix = sortParameters.output_prefix;
    
    %% Number of Subjects, components
    numSubjects = sortParameters.numSubjects;
    numComp = sortParameters.numComp;
    
    sortingCriteria = sortParameters.sortingCriteria;
    threshold = sortParameters.z_threshold;
    back_reconstruct_files = sortParameters.backReconstructFiles;
    
    % Histogram criteria
    if isfield(sortParameters, 'histogramCriteria')
        histogram_criteria = sortParameters.histogramCriteria;
    else
        histogram_criteria = 'feature';
    end
    
    
    %% Selected group names and indices
    selected_group_names = sortParameters.selGroupNames;
    selGroupsVal = sortParameters.selGroupsVal;
    
    %% Selected feature names and indices
    selectedFeature = sortParameters.selectedFeature;
    selectedFeatureVal = sortParameters.selectedFeatureVal;
    
    %selectedFeature = deblank(sortParameters.featureNames(selectedFeatureVal, :));
    
    % Feature data length
    %featureDataLength = sortParameters.featureDataLength;
    %dataLength = featureDataLength(comb_number).Length;
    mask_ind = sortParameters.mask_ind;
    normalizeVal = sortParameters.normalize;
    sortingCriteria = lower(sortingCriteria);
    all_comb = sortParameters.all_comb;
    
    % Load back Reconstruction parameters
    fileName = fullfile(outputDir, back_reconstruct_files(comb_number).name);
    load(fileName);
    
    if (comb_number == 1)
        if isfield(sortParameters, 'scaleCompFiles')
            scaleCompFile = sortParameters.scaleCompFiles(1).name;
            scaleCompFile = fullfile(outputDir, scaleCompFile);
            load(scaleCompFile, 'A');
        end
        
        icaFile = fullfile(outputDir, sortParameters.icaFiles(1).name);
        varsInICAFile = whos('-file', icaFile);
        checkCorrModalities = strmatch('corr_modalities', char(varsInICAFile.name), 'exact');
        if (~isempty(checkCorrModalities))
            load(icaFile, 'corr_modalities');
        end
        
    end
    
    if (~iscell(A))
        A = {A};
    else
        A = A(selectedFeatureVal);
    end
    
    group1_indices = ica_fuse_get_groupInd(selGroupsVal(1), numSubjects);
    
    group2_indices = ica_fuse_get_groupInd(selGroupsVal(2), numSubjects);
    
    % Print group names and feature names
    group1Name = deblank(selected_group_names(1, :));
    group2Name = deblank(selected_group_names(2, :));
    disp(['Selected groups are ', group1Name, ' and ', group2Name]);
    
    if length(selectedFeatureVal) == 1
        initStr = 'Selected feature is';
    else
        initStr = 'Selected features are';
    end
    
    for nn = 1:length(selectedFeatureVal)
        if nn == 1
            msgStr = [initStr, ' ', deblank(selectedFeature(nn, :))];
        else
            msgStr = [msgStr, ', ', deblank(selectedFeature(nn, :))];
        end
    end
    
    disp(msgStr);
    
    % Initialise ttest2 and kstest2 results
    %ttest2_results = zeros(1, numComp);
    sortCompData = repmat(struct('title', [], 'data', [], 'colorbarText', [], 'colorbarLim', []), 1, numComp);
    
    fprintf('\n');
    
    switch sortingCriteria
        
        case 'ttest2'
            % Calculate the ttest2 and kstest2
            %string1 = ['Calculating two sample t-test and  Kolmogorov-Smirnov test for the mixing coefficient between '];
            string1 = 'Calculating two sample t-test for the mixing coefficient between ';
            string2 = [deblank(selected_group_names(1, :)), ' and ', deblank(selected_group_names(2, :)), ' ...'];
            disp(str2mat(string1, string2));
            helpMsg = [string1, string2];
            
            if saveResults
                clear string1; clear string2;
                helpTitle = 'Ranking joint components by ttest2';
                helpHandle = helpdlg(helpMsg, helpTitle);
            end
            
            % Store the data information
            varStruct = repmat(struct('tag', '', 'value', []), 1, 2*length(A));
            ttest2_results = zeros(length(A), numComp);
            sorted_comp = zeros(length(A), numComp);
            
            for nA = 1:length(A)
                
                tmpA = A{nA};
                % Loop over components
                for nComp = 1:numComp
                    % calculate ttest2 and kstest2 for each component
                    pSignificance = ica_fuse_ttest2(tmpA(group1_indices, nComp), tmpA(group2_indices, nComp));
                    ttest2_results(nA, nComp) = pSignificance;
                end
                % end loop over components
                
                % Sort the components
                [tmpRes, sorted_comp(nA, :)] = sort(ttest2_results(nA, :));
                
                outFileName = [output_prefix, '_ttest2_mixing_coeff.txt'];
                varStruct(2*nA - 1).tag = 'Component Number';
                varStruct(2*nA - 1).value =  sorted_comp(nA, :);
                
                if (length(A) > 1)
                    varStruct(2*nA).tag = ['Two sample t-test (p-value) ', deblank(selectedFeature(nA, :))];
                else
                    varStruct(2*nA).tag = 'Two sample t-test (p-value)';
                end
                
                varStruct(2*nA).value = tmpRes;
                
            end
            
            % Form strings to print to a file
            string1 = 'Two sample t-test results for the mixing coefficient between ';
            string2 = [deblank(selected_group_names(1, :)), ' and ', deblank(selected_group_names(2, :))];
            titlePrint = [string1, string2];
            
            newA = cat(2, A{:});
            
            sorted_comp = sorted_comp(1, :);
            sorted_values = ttest2_results(1, sorted_comp);
            
            % Loop over components
            for nComp = 1:numComp
                
                compIndex = ica_fuse_returnFileIndex(sorted_comp(nComp));
                %                 sortCompData(nComp).title = ['Comp ', compIndex, ' ttest2 = ', ...
                %                         num2str(ttest2_results(nComp)), ' kstest2 = ', num2str(kstest2_results(nComp))];
                
                if (length(A) > 1)
                    str = [];
                    for nA = 1:length(A)
                        if (nA == 1)
                            str = ['Comp ', compIndex, ' ', deblank(selectedFeature(nA, :)), ' p = ', num2str(ttest2_results(nA, sorted_comp(nComp)), '%0.3f')];
                        else
                            str = [str, ' ', deblank(selectedFeature(nA, :)), ' p = ',  num2str(ttest2_results(nA, sorted_comp(nComp)), '%0.3f')];
                        end
                    end
                    
                    if (exist('corr_modalities', 'var'))
                        str = [str, ' r = ', num2str(corr_modalities(sorted_comp(nComp)), '%0.3f')];
                    end
                    
                    sortCompData(nComp).title = str;
                else
                    sortCompData(nComp).title = ['Comp ', compIndex, ' p = ', num2str(ttest2_results(sorted_comp(nComp)), '%0.3f')];
                end
                
                %mixingData = [group1_mixing_parameters(:, nComp), group2_mixing_parameters(:, nComp)];
                
                sortCompData(nComp).data(1).dat = newA(group1_indices, sorted_comp(nComp):numComp:end);
                sortCompData(nComp).data(2).dat = newA(group2_indices, sorted_comp(nComp):numComp:end);
                sortCompData(nComp).colorbarText = ['0'; '0'];
                sortCompData(nComp).colorbarLim = [0 0];
                clear mixingData;
                
            end
            
            % end loop over components
            plotType = 'ICA Loading';
            
        case 'divergence'
            
            if isfield(sortParameters, 'div_name')
                div_name = sortParameters.div_name;
                div_value = sortParameters.div_value;
            else
                [div_name, div_value] = ica_fuse_getDivergencePara;
            end
            
            helpMsg = ['Ranking joint components by ', div_name, ' divergence. Please wait...'];
            helpTitle = ['Ranking joint components by ', div_name, ' divergence'];
            disp(helpTitle);
            
            if saveResults
                helpHandle = helpdlg(helpMsg, helpTitle);
            end
            
            divergence_value = zeros(1, numComp);
            
            if ica_fuse_findstr(histogram_criteria, 'feature')
                if ~isfield(sortParameters, 'dataN')
                    currentCombNums = all_comb{comb_number};
                    stackInfo = ica_fuse_prepare_data('dataInfo', dataInfo, 'mask_ind', mask_ind, 'normalize_scheme', ...
                        normalizeVal, 'voxels', sortParameters.voxels, 'sel_groups', ...
                        selGroupsVal, 'sel_features', currentCombNums(selectedFeatureVal));
                    
                    sortParameters.dataN = stackInfo.data;
                    clear stackInfo;
                end
            end
            
            for nComp = 1:numComp
                
                if ica_fuse_findstr(histogram_criteria, 'feature')
                    %disp(['Using component ', num2str(nComp), ' voxels that are ranked from max to min as a mask to generate histograms ...']);
                    [hData, meanHist] = ica_fuse_calculate_histogram(sortParameters, histogram_criteria, ...
                        threshold, comb_number, selectedFeatureVal, nComp);
                    out1 = meanHist.data{1};
                    out2 = meanHist.data{2};
                    
                    % Get in
                    ind = find((out1(:) ~= 0) & (out2(:) ~=0));
                    
                    if isempty(ind)
                        error('Error:Histogram', 'Histograms of groups are disjoint with a Z-threshold of %s.', ...
                            num2str(threshold));
                    end
                    
                else
                    disp(['Using component ', num2str(nComp), ' voxels that are ranked from max to min to generate histograms ...']);
                    [hData] = ica_fuse_calculate_histogram(sortParameters, histogram_criteria, ...
                        threshold, comb_number, selectedFeatureVal, nComp);
                    out1 = hData.group(1).data;
                    out2 = hData.group(2).data;
                end
                
                % Get the KL divergence
                [divNames, divergence_value(nComp)] = ica_fuse_divergence(out1, out2, div_name, div_value);
                
                clear out1; clear out2;
            end
            
            sortParameters.dataN = [];
            clear data;
            plotType = 'ICA Loading';
            clear groups_icasig;
            
            % Sort the KL divergence
            [distances, sorted_comp] = sort(divergence_value);
            
            sorted_comp = sorted_comp(end:-1:1);
            
            %histData = histData(sorted_comp(1));
            
            newA = cat(2, A{:});
            
            % Sort the optimization data
            divergence_value = divergence_value(sorted_comp);
            
            sorted_values = divergence_value;
            
            if ~isempty(div_value)
                div_value_str = ['(', num2str(div_value), ')'];
            else
                div_value_str = '';
            end
            
            div_string = [div_name, div_value_str, ' Divergence'];
            
            div_string(1) = upper(div_string(1));
            
            % Form title here
            for nComp = 1:numComp
                compIndex = ica_fuse_returnFileIndex(sorted_comp(nComp));
                sortCompData(nComp).title = ['Component ', compIndex, ' ', div_string, ...
                    ' = ', num2str(divergence_value(nComp));];
                sortCompData(nComp).data(1).dat =  newA(group1_indices, sorted_comp(nComp):numComp:end);
                sortCompData(nComp).data(2).dat = newA(group2_indices, sorted_comp(nComp):numComp:end);
                sortCompData(nComp).colorbarText = ['0'; '0'];
                sortCompData(nComp).colorbarLim = [0 0];
                clear mixingData;
            end
            
            outFileName = [output_prefix, '_divergence_groups.txt'];
            % Store the data information
            numPara = 1;
            varStruct(numPara).tag = 'Component Number';
            varStruct(numPara).value = sorted_comp;
            
            numPara = numPara + 1;
            varStruct(numPara).tag = div_string;
            varStruct(numPara).value = divergence_value;
            
            % Form strings to print to a file
            string1 = ['Divergence results for groups ', deblank(selected_group_names(1, :)), ' and ', deblank(selected_group_names(2, :))];
            titlePrint = string1;
            
    end
    
    if exist('helpHandle', 'var')
        if ishandle(helpHandle)
            delete(helpHandle);
        end
    end
    
    sortResults.values = sorted_values;
    sortResults.sorted_comp = sorted_comp;
    sortResults.plotType = plotType;
    sortResults.sortCompData = sortCompData;
    sortResults.selGroupNames = selected_group_names;
    sortResults.selFeatureNames = selectedFeature;
    sortResults.selGroupsVal = selGroupsVal;
    sortResults.selFeaturesVal = selectedFeatureVal;
    
    %     if exist('histData', 'var')
    %         sortResults.histData = histData;
    %     end
    
    
    if saveResults
        % Form full file
        outFileName = fullfile(outputDir, outFileName);
        
        fprintf('\n');
        
        % Print to a file
        ica_fuse_printToFile(outFileName, varStruct, titlePrint, 'row_wise', 'append');
        
        disp(['Sorting information is stored in ', outFileName]);
        fprintf('\n');
    end
    
catch
    
    if exist('helpHandle', 'var')
        if ishandle(helpHandle)
            delete(helpHandle);
        end
    end
    ica_fuse_displayErrorMsg;
    
end
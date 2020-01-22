function pca_output = ica_fuse_parallelICA_dataReduction(featureData, modalities, numComp, reference, pca_type, type_parallel_ica)
% Data reduction step

pca_output = repmat(struct('dewhiteM', [], 'whiteM', [], 'whitesig', []), 1, length(modalities));

if ~exist('pca_type', 'var')
    pca_type = 'standard';
end


if ~exist('type_parallel_ica', 'var')
    type_parallel_ica = 'aa';
end

% Remove the mean only for fmri modality
for nModality = 1:length(modalities)
    fprintf('\n');
    currentData = featureData(nModality).data;
    if strcmpi(modalities{nModality}, 'fmri')
        currentData = currentData';
        disp(['Removing mean for feature ', featureData(nModality).feature_name, '...']);
        % Remove mean of fmri data for each subject
        try
            tempVar = ones(size(currentData, 1), 1)*mean(currentData);
            currentData = currentData - tempVar;
            clear tempVar;
        catch
            clear tempVar;
            % --slower way to zero mean data
            disp('Using slower but less memory constrained way to zero mean data');

            for j = 1:size(currentData, 2)
                currentData(:, j) = currentData(:,  j) - mean(currentData(:, j));
            end
        end
        currentData = currentData';
    elseif strcmpi(modalities{nModality}, 'gene')
        if strcmpi(type_parallel_ica, 'as')
            % Arrange SNP data as SNPS by subjects
            currentData = currentData';
        end
    end

    featureData(nModality).data = currentData;
    clear currentData;

    disp(['Doing data reduction for feature ', featureData(nModality).feature_name]);
    %%%%%%%%%%%%%%% PCA Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(pca_type, 'standard')

        % Calculate PCA
        [V, Lambda] = ica_fuse_v_pca(featureData(nModality).data, 1, numComp(nModality), 0, 'transpose', 'no');

    else
        % Calculate PCA with reference
        [V, Lambda] = ica_fuse_pca_reference(featureData(nModality).data, 1, numComp(nModality), reference, 0, 'transpose', 'no');
    end

    % Whitening
    [whitesig, whiteM, dewhiteM] = ica_fuse_v_whiten(featureData(nModality).data, V, Lambda, 'transpose');
    disp('Done');


    fprintf('\n');

    dataFiles(nModality).name = featureData(nModality).files;
    dataFiles(nModality).feature_name = featureData(nModality).feature_name;

    % Store PCA output
    pca_output(nModality).dewhiteM = dewhiteM;
    pca_output(nModality).whiteM = whiteM;
    pca_output(nModality).whitesig = whitesig;
    pca_output(nModality).modality = modalities{nModality};
    pca_output(nModality).feature_name = featureData(nModality).feature_name;

end
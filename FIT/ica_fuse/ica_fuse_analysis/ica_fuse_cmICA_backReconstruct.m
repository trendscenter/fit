function ica_fuse_cmICA_backReconstruct(cmICAInfo)
%% Run cmICA

ica_fuse_defaults;


%% Get params
outputDir = cmICAInfo.outputDir;
output_prefix = cmICAInfo.output_prefix;
numPC1 = cmICAInfo.numPC1;
numComp = cmICAInfo.numComp;
dataInfo = cmICAInfo.dataInfo;
sub_size = size(cmICAInfo.dataInfo(1).files, 1);

icaFile = fullfile(outputDir, [output_prefix, '_joint_cmica_ica.mat']);
load(icaFile, 'A', 'W', 'icasig');

disp('Back-reconstructing components');

%% Loop over modalities
endC = 0;
for nModality = 1:length(dataInfo)
    startC = endC + 1;
    endC = endC + numComp;
    AM = A(startC:endC, :);
    
    featureName = dataInfo(nModality).feature_name;
    modality  = dataInfo(nModality).modality;
    disp(['Loading components of ', featureName]);
    in_file = fullfile(outputDir, [output_prefix, '_', featureName, '_joint_cmica_pca_r2-1.mat']);
    load(in_file, 'whiteM', 'dewhiteM');
    
    whiteM2 = whiteM;
    whiteM2 = reshape(whiteM2, numComp, numPC1, sub_size);
    dewhiteM2 = dewhiteM;
    dewhiteM2 = reshape(dewhiteM2, numPC1, sub_size, numComp);
    dewhiteM2 = permute(dewhiteM2, [1 3 2]);
    clear dewhiteM whiteM pcasig;
    %% Loop over subjects
    for ii = 1:sub_size
        
        in_file = fullfile(outputDir, [output_prefix, '_', featureName, '_joint_cmica_pca_r1-', num2str(ii), '.mat']);
        load(in_file, 'dewhiteM', 'pcasig', 'Lambda');
        
        % gica1
        ICdinfo = pinv(dewhiteM2(:, :, ii));
        TCdinfo = squeeze(dewhiteM2(:,:,ii));
        
        ic = pinv(AM)*ICdinfo*pcasig';
        if strcmpi(modality, 'fmri')
            rc = pcasig*Lambda*TCdinfo*AM;
        else
            rc = dewhiteM*TCdinfo*AM;
        end
        output_file = fullfile(outputDir, [output_prefix, '_', featureName, '_joint_cmica_br_r1-', num2str(ii), '.mat']);
        save(output_file, 'ic', 'rc');
        
    end
    %% End loop over subjects
    
    
end
%% End loop over modalities

disp('Done');
fprintf('\n');
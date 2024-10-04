function ica_fuse_cmICA_scaleComponents(cmICAInfo)
%% Run cmICA

ica_fuse_defaults;


%% Get params
outputDir = cmICAInfo.outputDir;
output_prefix = cmICAInfo.output_prefix;
dataInfo = cmICAInfo.dataInfo;
sub_size = size(cmICAInfo.dataInfo(1).files, 1);

disp('Scaling components');
for nModality = 1:length(dataInfo)
    featureName = dataInfo(nModality).feature_name;
    modality  = dataInfo(nModality).modality;
    disp(['Loading components of ', featureName]);
    pca_file = fullfile(outputDir, [output_prefix, '_', featureName, '_joint_cmica_pca_r2-1.mat']);
    pcinfo = load(pca_file);
    gm_mask_inds = pcinfo.gm_mask_inds;
    if (strcmpi(modality, 'dti'))
        wm_mask_inds = pcinfo.wm_mask_inds;
    end
    
    for ii = 1:sub_size
        br_file = fullfile(outputDir, [output_prefix, '_', featureName, '_joint_cmica_br_r1-', num2str(ii), '.mat']);
        load(br_file, 'ic', 'rc');
        if (ii == 1)
            icm = zeros(size(ic));
            rcm = zeros(size(rc));
        end
        icm = icm + ic;
        rcm = rcm + rc;
    end
    icm = icm/sub_size;
    rcm = rcm/sub_size;
    
    
    if (strcmpi(modality, 'fmri'))
        writeV = ica_fuse_spm_vol(deblank(cmICAInfo.dataInfo(nModality).files(1, :)));
        writeV = writeV(1);
        ic_4d = zeros(prod(writeV(1).dim(1:3)), size(icm, 1));
        ic_4d(gm_mask_inds, :) = icm';
        rc_4d = zeros(prod(writeV(1).dim(1:3)), size(icm, 1));
        rc_4d(gm_mask_inds, :) = rcm;
    else
        % dti
        writeV = ica_fuse_spm_vol(cmICAInfo.dataInfo(nModality).mask);
        writeV = writeV(1);
        ic_4d = zeros(prod(writeV(1).dim(1:3)), size(icm, 1));
        ic_4d(gm_mask_inds, :) = icm';
        rc_4d = zeros(prod(writeV(1).dim(1:3)), size(icm, 1));
        rc_4d(wm_mask_inds, :) = rcm;
    end
    
    % reshape data
    ic_4d = reshape(ic_4d, [writeV(1).dim(1:3), size(icm, 1)]);
    rc_4d = reshape(rc_4d, [writeV(1).dim(1:3), size(icm, 1)]);
    
    disp(['Writing components for modality ', featureName]);
    out_fname1 = [output_prefix, '_', featureName, '_joint_cmica_comps.nii'];
    out_fname2 = [output_prefix, '_', featureName, '_joint_cmica_conn.nii'];
    
    writeV = repmat(writeV, size(icm, 1), 1);
    ica_fuse_write_nifti_data(fullfile(outputDir, out_fname1), writeV, ic_4d, 'IC');
    
    ica_fuse_write_nifti_data(fullfile(outputDir, out_fname2), writeV, rc_4d, 'CONN');
    fprintf('Done\n');
    clear icm rcm;
    
    outputFiles(nModality).name = out_fname1;
    outputFiles(nModality).conn_name = out_fname2;
    outputFiles(nModality).feature_name = featureName;
    outputFiles(nModality).modality = modality;
    
    
end

cmICAInfo.outputFiles = outputFiles;
fusionFile = fullfile(outputDir, [output_prefix, '_joint_cmica_info.mat']);
save(fusionFile, 'cmICAInfo');
fprintf('\n');
disp(['Fusion information is saved in file: ', fusionFile]);


fprintf('Done scaling components \n');

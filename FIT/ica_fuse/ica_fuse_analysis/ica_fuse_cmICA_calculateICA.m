function ica_fuse_cmICA_calculateICA(cmICAInfo)
%% Run cmICA

ica_fuse_defaults;


%% Get params
outputDir = cmICAInfo.outputDir;
output_prefix = cmICAInfo.output_prefix;
numPC1 = cmICAInfo.numPC1;
numComp = cmICAInfo.numComp;
dataInfo = cmICAInfo.dataInfo;
algorithm = lower(cmICAInfo.algorithm);

disp('Calculating joint cmICA');

%% Stack data temporally
pcasigM = cell(1, length(dataInfo));
for nModality = 1:length(dataInfo)
    featureName = dataInfo(nModality).feature_name;
    in_file = fullfile(outputDir, [output_prefix, '_', featureName, '_joint_cmica_pca_r2-1.mat']);
    load(in_file, 'pcasig');
    pcasigM{nModality} = pcasig;
end

pcasig = cat(2, pcasigM{:});
pcasig = pcasig';

[icaAlgo, W, A, icasig] = ica_fuse_icaAlgorithm(algorithm, pcasig);

numOfIC = size(icasig, 1);

% or  icasso
%[A, W, icasig, sR] = icasso_tmp(pcasig,numOfIC); % 'posact' is on, works better than the codes below!!

for compNum = 1:size(icasig,1)
    v = ica_fuse_recenter_image(icasig(compNum, :));
    skew = ica_fuse_skewness(v);
    clear v;
    if (sign(skew) == -1)
        disp(['Changing sign of component ',num2str(compNum)]);
        icasig(compNum, :) = icasig(compNum, :)*-1;
        A(:,compNum) = A(:,compNum)*-1;
        W(compNum,:) = W(compNum,:)*-1;
    end
end

%% Sort components by mean variance
meanvar = zeros(1, numOfIC);
for i=1:numOfIC
    compproj=A(:,i)*icasig(i,:);
    meanvar(i)=var(compproj(:));
end

[sortvar, windex] = sort(meanvar,'descend');
icasig = icasig(windex,:);
W = W(windex,:);% reorder the weight matrix
A = A(:,windex);

% outputFile = fullfile(outputDir, [output_prefix, 'test_joint_cmica_ica.mat']); %lei
outputFile = fullfile(outputDir, [output_prefix, '_joint_cmica_ica.mat']);
save(outputFile, 'A', 'W', 'icasig');

%add by lei
gm_mask_inds = find(ica_fuse_read_data(cmICAInfo.dataInfo(1).mask)~=0);
writeV = ica_fuse_spm_vol(deblank(cmICAInfo.dataInfo(1).mask));
ic_4d = zeros(prod(writeV(1).dim(1:3)), size(icasig, 1));
ic_4d(gm_mask_inds, :) = icasig';
ic_4d = reshape(ic_4d, [writeV(1).dim(1:3), size(icasig, 1)]);
disp(['Writing components for aggregated maps']);
out_fname = [output_prefix, '_agg_joint_cmica_comps2.nii'];
writeV = repmat(writeV, size(icasig, 1), 1);
ica_fuse_write_nifti_data(fullfile(outputDir, out_fname), writeV, ic_4d, 'IC');
%%%

disp('Done');
fprintf('\n');

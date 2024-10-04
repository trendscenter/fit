function  ica_fuse_cmICA_dataReduction(cmICAInfo)
%% Data reduction cmICA

ica_fuse_defaults;


%% Get params
outputDir = cmICAInfo.outputDir;
output_prefix = cmICAInfo.output_prefix;
numPC1 = cmICAInfo.numPC1;
numComp = cmICAInfo.numComp;
dataInfo = cmICAInfo.dataInfo;

disp('Computing data reduction step');

for nModality = 1:length(dataInfo)
    featureName = dataInfo(nModality).feature_name;
    dataInfo(nModality).outputDir = outputDir;
    dataInfo(nModality).prefix = output_prefix;
    disp(['Applying two step pca on feature ', featureName]);
    
    apply_gpca(dataInfo(nModality), numPC1, numComp);
    
end

disp('Done');
fprintf('\n');


function apply_gpca(dataInfo, numPC1, numComp)

ica_fuse_defaults;
global GREY_MATTER;
global WHITE_MATTER;

prefix = dataInfo.prefix;
files = cellstr(dataInfo.files);
sub_size = length(files);
modalityName = dataInfo.modality;
outputDir = dataInfo.outputDir;
feature_name = dataInfo.feature_name;

if (strcmpi(modalityName, 'fmri'))
    refImage = deblank(dataInfo.files(1, :));
else
    refImage = deblank(dataInfo.mask);
end
[inDir, fN, extn] = fileparts(refImage);
refImage = ica_fuse_fullFile('directory', inDir, 'files', ica_fuse_listFiles_inDir(inDir, [fN, extn]));
refImage = deblank(refImage(1, :));

rdata = ica_fuse_spm_vol(refImage);
rdata = rdata(1);
gdata = ica_fuse_spm_vol(GREY_MATTER);
gdata = gdata(1);
wdata = ica_fuse_spm_vol(WHITE_MATTER);
wdata = wdata(1);

% Load gray and white matter files

if all(rdata(1).dim == gdata.dim)
    resliced_gm_files = GREY_MATTER;
else
    resliced_gm_files = ica_fuse_spm_corgeister(refImage, GREY_MATTER, GREY_MATTER, outputDir);
end

if all(rdata(1).dim == wdata.dim)
    resliced_wm_files = WHITE_MATTER;
else
    resliced_wm_files = ica_fuse_spm_corgeister(refImage, WHITE_MATTER, WHITE_MATTER, outputDir);
end



gmmask = ica_fuse_loadData(resliced_gm_files);
wmmask = ica_fuse_loadData(resliced_wm_files);
gmmask = gmmask(:, :, :, 1);
wmmask = wmmask(:, :, :, 1);
gmmask_ind = find(gmmask~=0);
gmmask_len = length(gmmask_ind);
wmmask_ind = find(wmmask~=0);
wmmask_len = length(wmmask_ind);

%if (strcmpi(modalityName, 'fmri'))

pc1_data = cell(1, sub_size);
disp(['Doing First Step PCA on feature ', feature_name]);
for ii = 1:sub_size
    disp(['Loading Subject', num2str(ii), '...'])
    currentFiles = files{ii};
    if (size(currentFiles, 1) == 1)
        [inDir, fN, extn] = fileparts(files{ii});
        sub_files = ica_fuse_fullFile('directory', inDir, 'files', ica_fuse_listFiles_inDir(inDir, [fN, extn]));
    else
        sub_files = currentFiles;
    end
    
    if (strcmpi(modalityName, 'fmri'))
        %% FMRI
        if (~isempty(dataInfo.mask))
            mask_inds = ica_fuse_loadData(dataInfo.mask);
            mask_inds = find(abs(mask_inds) > eps);
            [~, gmmask_loc] = intersect(mask_inds, gmmask_ind);
            gm_mask_inds = mask_inds(gmmask_loc);
        else
            gm_mask_inds = gmmask_ind;
        end
        
        
        data = ica_fuse_read_data(sub_files, [], gm_mask_inds);
        data = icatb_preproc_data4cmica(data);
        data(isnan(data)) = 0;
        
        disp('Doing standard PCA!');
        [V, Lambda] = ica_fuse_svd(data, numPC1);
        
    else
        
        %% DTI
        mask_inds = ica_fuse_loadData(dataInfo.mask);
        mask_inds = find(abs(mask_inds) > eps);
        [~,gmmask_loc] = intersect(mask_inds, gmmask_ind);
        [~, wmmask_loc] = intersect(mask_inds, wmmask_ind);
        
        fdt_matrix3 = load(currentFiles);  % taking a while to load, or load fdt_matrix2 depending on masking
        disp('create connectivity matrix');
        data = sparse(fdt_matrix3(:, 2), fdt_matrix3(:, 1), fdt_matrix3(:, 3));
        data = data + data';
        data = data(gmmask_loc, wmmask_loc);
        data = spfun(@(x) log10(x+1), data);  %sparse compute
        
        disp('1st data reduction - sparse PCA');
        pw = 2;  % # power
        [U, Lambda, V] = fsvd(data, numPC1, pw);
        Lambda = Lambda.*diag(std(U)); % adjusted, the original Lambda ~= sqrtm(Lambda_eig), un order to make std(pcasig)==1 (sphering)
        
        gm_mask_inds = mask_inds(gmmask_loc);
        wm_mask_inds = mask_inds(wmmask_loc);
        
    end
    
    %     else
    %         %% DTI
    %         if (~isempty(dataInfo.mask))
    %             mask_inds = ica_fuse_loadData(dataInfo.mask);
    %             mask_inds = find(abs(mask_inds) > eps);
    %         else
    %             mask_inds = (1:size(data, 1)*size(data, 2)*size(data, 3));
    %         end
    %
    %         data = ica_fuse_read_data(sub_files, [], mask_inds);
    %         data(isnan(data)) = 0;
    %         [~, gmmask_loc]=intersect(mask_inds, gmmask_ind);
    %         [~, wmmask_loc]=intersect(mask_inds, wmmask_ind);
    %
    %         [V, Lambda] = compute_pca_conn(data, gmmask_loc, wmmask_loc, numComp);
    %
    %         gmmask = gmmask_ind(gmmask_loc);
    %         wmmask = wmmask_ind(wmmask_loc);
    
    % end
    
    whiteM = sqrtm(abs(Lambda)) \ V';
    dewhiteM = V * sqrtm(abs(Lambda));
    pcasig = data*whiteM';
    pc1_data{ii} = pcasig;
    out_file = fullfile(outputDir, [prefix, '_', feature_name, '_joint_cmica_pca_r1-', num2str(ii), '.mat']);
    disp(['Saving file ', out_file, ' ...']);
    
    save(out_file, 'Lambda', 'V', 'whiteM', 'dewhiteM', 'pcasig');
    
end
disp('Done First Step PCA');

disp('Loading PCA from fMRI');
pcasig2 = cat(2, pc1_data{:});
clear pc1_data;

disp(['Doing Second Step PCA on ', feature_name, ' ...']);
[V, Lambda] = ica_fuse_svd(pcasig2, numComp);
% Whitening and de-whitening matrix
whiteM = sqrtm(abs(Lambda)) \ V';
dewhiteM = V * sqrtm(abs(Lambda));
pcasig = pcasig2*whiteM';  %dfmri_matrix rotated already!!!!!!!!!!!!!!!!! column zero means

out_file = fullfile(outputDir, [prefix, '_', feature_name,  '_joint_cmica_pca_r2-1.mat']);
disp(['Saving file ', out_file, ' ...']);
save(out_file, 'Lambda', 'V', 'whiteM', 'dewhiteM', 'pcasig', 'gm_mask_inds');

if (exist('wm_mask_inds', 'var'))
    save(out_file,  'wm_mask_inds', '-append');
end


function out = icatb_preproc_data4cmica(out)
%% Use 50 voxels at a time for variance or intensity normalization
blockSize = 50;
voxels = size(out, 1);
nLoops = ceil(voxels / blockSize);
endT = 0;
for n = 1:nLoops
    startT = endT + 1;
    endT = endT + blockSize;
    if (endT > voxels)
        endT = voxels;
    end
    
    tmp = out(startT:endT, :)';
    tmp = detrend(tmp,'constant');  %lei modified from original 'vn', w.o. global signal removal.
    tmp = tmp.*repmat(1./(std(tmp)+eps), size(tmp, 1), 1);
    out(startT:endT, :) = tmp';
    clear tmp;
end


function [V2, D2, U2]  = compute_pca_conn(data, grey_mask, white_mask, numComp)
%% Compute pca on connectivity matrix

max_iter = 1000;
tol = 1e-4;
residual_err = 1;
verbose = 1;

C = randn(length(white_mask), numComp);
C = norm2(C);

count = 0;

while ((residual_err > tol) && (count <= max_iter))
    
    count = count + 1;
    tmpS = data(grey_mask, :)*(data(white_mask, :)'*C);
    dd = (pinv(tmpS)*data(grey_mask, :))*data(white_mask, :)';
    
    C_old = C;
    
    C = dd';
    
    C = norm2(C);
    
    residual_err = norm_resid(C, C_old);
    
    if (mod(count, 5) == 0)
        if (verbose)
            disp(['Step No: ', num2str(count), ' Norm of residual error: ', num2str(residual_err, '%0.6f')]);
        end
    end
    
end

C = orth(C);

cov_m = data(grey_mask, :)*(data(white_mask, :)'*C);

cov_m = (cov_m'*cov_m)/(size(cov_m, 1) - 1 );

[V2, D2] = eig(cov_m, 'nobalance');

V2 = C*V2;
U2 = data(grey_mask, :)*(data(white_mask, :)'*V2);

%V2 = V2(:, end:-1:1);
%U2 = U2(:, end:-1:1);

U2 = U2*diag(1./sqrt(sum(U2.^2)));


function C = norm2(C)
%% Normalize transformation matrix using slower way
%

for nC = 1:size(C, 1)
    C(nC, :) =  C(nC, :) ./ norm(C(nC, :), 2);
end

function residual_err = norm_resid(C, C_old)
%% Use norm2 of residual error
%

residual_err = 0;
for nC = 1:size(C, 1)
    res = C(nC, :) - C_old(nC, :);
    residual_err = residual_err + sum(res.^2);
end

residual_err = sqrt(residual_err);


function [U,S,V] = fsvd(A, k, i, usePowerMethod)
% FSVD Fast Singular Value Decomposition
%
%   [U,S,V] = FSVD(A,k,i,usePowerMethod) computes the truncated singular
%   value decomposition of the input matrix A upto rank k using i levels of
%   Krylov method as given in [1], p. 3.
%
%   If usePowerMethod is given as true, then only exponent i is used (i.e.
%   as power method). See [2] p.9, Randomized PCA algorithm for details.
%
%   [1] Halko, N., Martinsson, P. G., Shkolnisky, Y., & Tygert, M. (2010).
%   An algorithm for the principal component analysis of large data sets.
%   Arxiv preprint arXiv:1007.5510, 0526. Retrieved April 1, 2011, from
%   http://arxiv.org/abs/1007.5510.
%
%   [2] Halko, N., Martinsson, P. G., & Tropp, J. A. (2009). Finding
%   structure with randomness: Probabilistic algorithms for constructing
%   approximate matrix decompositions. Arxiv preprint arXiv:0909.4061.
%   Retrieved April 1, 2011, from http://arxiv.org/abs/0909.4061.
%
%   See also SVD.
%
%   Copyright 2011 Ismail Ari, http://ismailari.com.

if nargin < 3
    i = 1;
end

% Take (conjugate) transpose if necessary. It makes H smaller thus
% leading the computations to be faster
if size(A,1) < size(A,2)
    A = A';
    isTransposed = true;
else
    isTransposed = false;
end

n = size(A,2);
l = k + 2;

% Form a real n×l matrix G whose entries are iid Gaussian r.v.s of zero
% mean and unit variance
G = randn(n,l);


if nargin >= 4 && usePowerMethod
    % Use only the given exponent
    H = A*G;
    for j = 2:i+1
        H = A * (A'*H);
    end
else
    % Compute the m×l matrices H^{(0)}, ..., H^{(i)}
    % Note that this is done implicitly in each iteration below.
    H = cell(1,i+1);
    H{1} = A*G;
    for j = 2:i+1
        H{j} = A * (A'*H{j-1});
    end
    
    % Form the m×((i+1)l) matrix H
    H = cell2mat(H);
end

% Using the pivoted QR-decomposiion, form a real m×((i+1)l) matrix Q
% whose columns are orthonormal, s.t. there exists a real
% ((i+1)l)×((i+1)l) matrix R for which H = QR.
% XXX: Buradaki column pivoting ile yapılmayan hali.
[Q,~] = qr(H,0);

% Compute the n×((i+1)l) product matrix T = A^T Q
T = A'*Q;

% Form an SVD of T
[Vt, St, W] = svd(T,'econ');

% Compute the m×((i+1)l) product matrix
Ut = Q*W;

% Retrieve the leftmost m×k block U of Ut, the leftmost n×k block V of
% Vt, and the leftmost uppermost k×k block S of St. The product U S V^T
% then approxiamtes A.

if isTransposed
    V = Ut(:,1:k);
    U = Vt(:,1:k);
else
    U = Ut(:,1:k);
    V = Vt(:,1:k);
end
S = St(1:k,1:k);

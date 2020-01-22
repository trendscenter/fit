function [comp_est, mdl, aic, kic] = ica_fuse_estim_dim(data, varargin)
% function [comp_est, mdl, aic, kic] = ica_fuse_estim_dim(files, maskvec)
% Select the order of the multivariate data using Information theoretic
% criteria with the option of sample dependence correction;
% Inputs:
%
% 1. data - Masked data of dimensions time points by voxels.
% 2. varargin - Variable no. of arguments passed in pairs.
%
% Output:
% 1. comp_est:  estimated order using MDL
% 2. mdl - MDL vector
% 3. aic - AIC vector
% 4. kic - KIC vector

% Please cite the following paper if you use this code for publication
% Y.-O. Li, T. Adali and V. D. Calhoun, "Estimating the number of independent
% components for fMRI data," Human Brain Mapping, vol. 28, pp. 1251--66, 2007


%COPYRIGHT NOTICE
%This function is a part of Fusion software library
%Copyright (C) 2003-2009
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


doSampling = 1;
scaleVal = 1;

%% Select files
if (~exist('data', 'var'))
    data = ica_fuse_selectEntry('typeSelection', 'multiple', 'filter', '*.img;*.nii;*.dat;*.asc', 'title', 'Select files ...');
end

drawnow;

%% Check arguments
for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'maskvec'))
        maskvec = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'dosampling'))
        doSampling = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'scalefactor'))
        scaleVal = varargin{n + 1};
    end
end


%% Load data
if (ischar(data))
    
    files = data;
    clear data;
    [dd, pp, extn] = fileparts(deblank(files(1, :)));
    data = ica_fuse_loadData(files);
    
    if (~strcmpi(extn, '.img') &&  ~strcmpi(extn, '.nii'))
        % Treat data as ascii
        data = squeeze(data(:, 2, :))';
        dimInfo = size(data, 2);
    else
        % Reshape .img or .nii files as Time by voxels
        dimInfo = [size(data, 1), size(data, 2), size(data, 3)];
        data = reshape(data, prod(dimInfo), size(data, 4))';
    end
    
end

data(isnan(data)) = 0;

%% Compute mask
if (~exist('maskvec', 'var'))
    maskvec = (prod(data) ~= 0);
end

if (ischar(maskvec))
    maskvec = deblank(maskvec(1, :));
    [dd, pp, extn] = fileparts(maskvec);
    maskvec = ica_fuse_loadData(maskvec);
    if (~strcmpi(extn, '.img') &&  ~strcmpi(extn, '.nii'))
        maskvec = squeeze(maskvec(:, 2, :));
    end
end

% Convert mask to logical vector
maskvec = (maskvec ~= 0);

if ~exist('dimInfo', 'var')
    dimInfo = size(maskvec);
end

maskvec = maskvec(:);
verbose = 1;
tdim = size(data, 1);

%% Mask data if possible
if (size(data, 2) ~= length(find(maskvec ~= 0)))
    data = data(:, maskvec);
end

% %remove mean
for ii = 1:size(data, 1)
    data(ii, :) = data(ii, :) - mean(data(ii, :));
end

if (doSampling)
    
    xdim = dimInfo(1);
    
    ydim = 1;
    if (length(dimInfo) >= 2)
        ydim = dimInfo(2);
    end
    
    zdim = 1;
    if (length(dimInfo) >= 3)
        zdim = dimInfo(3);
    end
    
    %   (data, 1, maskvec, dimv, 1)
    indp_eff = 1;
    
    %% Perform variance normalization
    if verbose
        fprintf('\n Performing variance normalization ...');
    end
    
    std_x = std(data');
    data = diag(1./std_x)*data;
    %% (completed)
    
    fprintf('\n P1:');
    [EigenVectors, EigenValues, dataN] = ica_fuse_pcsquash(data, tdim, 'no'); % x: (T x XYZ)
    
    
    %% Select Gaussian principal components for effectively i.i.d. sample
    %% estimation
    if verbose
        fprintf('\n Selecting Gaussian principal components ...');
    end
    for ik = 1:tdim
        if EigenValues(ik) > mean(EigenValues);
            kurtv1(ik) = 1000;
        else
            kurtv1(ik) = kurtn(dataN(ik, :));
        end
        
    end
    idx_gauss = find(kurtv1<0.3); % 0.3 is heuristic threshold to decide
    % (approx.) Gaussian distribution
    num_gauss_pc = length(idx_gauss);
    if num_gauss_pc > 2
        idx = idx_gauss([num_gauss_pc:-round((num_gauss_pc)/2)+1:1]);
    elseif num_gauss_pc == 0
        idx = [tdim-2:tdim];
    else
        idx = idx_gauss;
    end
    %% (completed)
    
    %% Estimate the subsampling depth for effectively i.i.d. samples
    if verbose
        fprintf('\n\n Estimate effectively i.i.d. samples: ');
    end
    
    mask_ND = reshape(maskvec, xdim, ydim, zdim);
    ms = length(idx);
    s = zeros(ms,1);
    for i = 1:ms
        x_single = zeros(xdim*ydim*zdim, 1);
        x_single(maskvec) = dataN(idx(i), :);
        x_single = reshape(x_single, xdim, ydim, zdim);
        [s(i)] = est_indp_sp(x_single);
        dim_n = prod(size(size(x_single)));
        clear x_single;
    end
    clear dataN;
    fprintf('\n');
    s1 = round(mean(s));
    if floor((sum(maskvec)/(1*tdim))^(1/dim_n)) < s1
        s1 = floor((sum(maskvec)/(1*tdim))^(1/dim_n));
    end
    N = round(sum(maskvec)/s1^dim_n);
    %% (completed)
    
    %% Use the subsampled dataset to calculate eigen values
    if verbose
        fprintf('\n Perform EVD on the effectively i.i.d. samples ...');
    end
    if s1 ~= 1
        mask_s = subsampling(mask_ND, s1, [1,1,1]);
        mask_s_1d = reshape(mask_s, 1, prod(size(mask_s)));
        dat = zeros(tdim, length(find(mask_s_1d == 1)));
        for i = 1:tdim
            x_single = zeros(xdim*ydim*zdim, 1);
            x_single(maskvec) = data(i, :);
            x_single = reshape(x_single, xdim, ydim, zdim);
            dat0 = subsampling(x_single, s1, [1, 1, 1]);
            dat0 = reshape(dat0,1,prod(size(dat0)));
            dat(i,:) = dat0(mask_s_1d);
            clear x_single;
        end
        clear data;
        %% Perform variance normalization
        std_dat = std(dat');
        dat = diag(1./std_dat)*dat;
        %% (completed)
        fprintf('\n P2:');
        [EigenVectors, EigenValues] = ica_fuse_pcsquash(dat, tdim);
        %lam = EigenValues;
    end
    lam = EigenValues;
    clear dat;
    
    if verbose
        fprintf('\n Effective number of i.i.d. samples: %d \n',N);
    end
    
    if size(lam,1)==1
        lam = fliplr(sort(lam));
    else
        lam = flipud(sort(lam))';
    end
    
    %% Make eigen spectrum adjustment
    if verbose
        fprintf('\n Perform eigen spectrum adjustment ...');
    end
    lam = eigensp_adj(lam,N,length(lam));
    %% (completed)
    
    if sum(imag(lam))
        error('Invalid eigen value found for the subsampled data.');
    end
    
    %% Correction on the ill-conditioned results (when tdim is large, some
    %% least significant eigenvalues become small negative numbers)
    lam(real(lam) <= 0) = min(lam(real(lam)>=0));
    
else
    
    N = round(size(data, 2)*scaleVal);
    lam = eig(data*data');
    
end

if verbose
    fprintf('\n Estimating the dimension ...');
end
p = tdim;
aic = zeros(1, p - 1);
kic = zeros(1, p - 1);
mdl = zeros(1, p - 1);
for k = 1:p-1
    LH(k) = log(prod( lam(k+1:end).^(1/(p-k)) )/mean(lam(k+1:end)));
    mlh(k) = 0.5*N*(p-k)*LH(k);
    df(k) = 1 + 0.5*k*(2*p-k+1);
    aic(k) =  -2*mlh(k) + 2*df(k);
    kic(k) =  -2*mlh(k) + 3*df(k);
    mdl(k) =  -mlh(k) + 0.5*df(k)*log(N);
end

% Find the first local minimum of each ITC
itc = zeros(3, length(mdl));
itc(1,:) = aic;
itc(2,:) = kic;
itc(3,:) = mdl;

for i = 1:size(itc, 1)
    dlap = squeeze(itc(i,2:end)-itc(i,1:end-1));
    a = find(dlap > 0);
    if isempty(a)
        rst_dim(i) = length(squeeze(itc(i,:)));
    else
        rst_dim(i) = a(1);
    end
end

% estimated components using MDL
comp_est = rst_dim(3);


%%%%%%% Sub functions for estimate dimension %%%%%%%%

function out = subsampling(x,s,x0)
% Subsampling the data evenly with space 's'

n = size(x);

if max(size(n)) == 2 & min(n) == 1  % 1D
    out = x([x0(1):s:max(n)]);
    
elseif max(size(n)) == 2 & min(n) ~= 1  % 2D
    out = x([x0(1):s:n(1)],[x0(2):s:n(2)]);
    
elseif max(size(n)) == 3 & min(n) ~= 1  % 3D
    out = x([x0(1):s:n(1)],[x0(2):s:n(2)],[x0(3):s:n(3)]);
    
else
    error('Unrecognized matrix dimension!(subsampling)');
    
end


function [s] = est_indp_sp(x)
% estimate the effective number of independent samples based on the maximum
% entropy rate principle of stationary random process


dimv = size(x);
s0 = 0;
fprintf('\n');
for j = 1:min(dimv)-1
    x_sb = subsampling(x,j,[1,1,1]);
    if j == 1
        fprintf('\n Estimating the entropy rate of the Gaussian component with subsampling depth %d,',j);
    else
        fprintf(' %d,',j);
    end
    entrate_m = entrate_sp(x_sb,1);
    
    ent_ref = 1.41;
    if entrate_m > ent_ref
        s0 = j; break;
    end
end
fprintf(' Done;');
if s0 == 0
    error('Ill conditioned data, can not estimate independent samples.(est_indp_sp)');
else
    s = s0;
end


function out = entrate_sp(x, sm_window)
% Calculate the entropy rate of a stationary Gaussian random process using
% spectrum estimation with smoothing window

n = size(x);

% Normalize x_sb to be unit variance
x_std = std(reshape(x,prod(n),1));
if x_std < 1e-10; x_std = 1e-10; end;
x = x/x_std;

if( sm_window == 1)
    
    M = ceil(n/10);
    
    if(max(size(n) >= 3))
        parzen_w_3 = zeros(2*n(3)-1,1);
        parzen_w_3(n(3)-M(3):n(3)+M(3)) = parzen_win(2*M(3)+1);
    end
    
    if(max(size(n) >= 2))
        parzen_w_2 = zeros(2*n(2)-1,1);
        parzen_w_2(n(2)-M(2):n(2)+M(2)) = parzen_win(2*M(2)+1);
    end
    
    if(max(size(n) >= 1))
        parzen_w_1 = zeros(2*n(1)-1,1);
        parzen_w_1(n(1)-M(1):n(1)+M(1)) = parzen_win(2*M(1)+1);
    end
    
end

if max(size(n)) == 2 & min(n) == 1  % 1D
    xc = xcorr(x,'unbiased');
    xc = xc.*parzen_w;
    xf = fftshift(fft(xc));
    
elseif max(size(n)) == 2 & min(n) ~= 1  % 2D
    xc = xcorr2(x); % default option: computes raw correlations with NO
    % normalization -- Matlab help on xcorr
    
    % Bias correction
    v1 = [1:n(1),n(1)-1:-1:1];
    v2 = [1:n(2),n(2)-1:-1:1];
    
    vd = v1'*v2;
    xc = xc./vd;
    parzen_window_2D = parzen_w_1*parzen_w_2';
    xc = xc.*parzen_window_2D;
    xf = fftshift(fft2(xc));
    
elseif max(size(n)) == 3 & min(n) ~= 1  % 3D
    xc = zeros(2*n-1);
    for m3 = 0:n(3)-1
        temp = zeros(2*n(1:2)-1);
        for k = 1:n(3)-m3
            temp = temp + xcorr2(x(:,:,k+m3),x(:,:,k)); % default option:
            % computes raw correlations with NO normalization
            % -- Matlab help on xcorr
        end
        xc(:,:,n(3)-m3) = temp;
        xc(:,:,n(3)+m3) = temp;
    end
    
    % Bias correction
    v1 = [1:n(1),n(1)-1:-1:1];
    v2 = [1:n(2),n(2)-1:-1:1];
    v3 = [n(3):-1:1];
    
    vd = v1'*v2;
    vcu = zeros(2*n-1);
    for m3 = 0:n(3)-1
        vcu(:,:,n(3)-m3) = vd*v3(m3+1);
        vcu(:,:,n(3)+m3) = vd*v3(m3+1);
    end
    
    xc = xc./vcu;
    
    parzen_window_2D = parzen_w_1*parzen_w_2';
    for m3 = 0:n(3)-1
        parzen_window_3D(:,:,n(3)-m3) = parzen_window_2D*parzen_w_3(n(3)-m3);
        parzen_window_3D(:,:,n(3)+m3) = parzen_window_2D*parzen_w_3(n(3)+m3);
    end
    xc = xc.*parzen_window_3D;
    
    xf = fftshift(fftn(xc));
    
else
    error('Unrecognized matrix dimension.');
    
end

xf = abs(xf);
xf(xf<1e-4) = 1e-4;
out = 0.5*log(2*pi*exp(1)) + sumN(log(abs((xf))))/2/sumN(abs(xf));


function w = parzen_win(n)
% PARZENWIN Parzen window.
%   PARZENWIN(N) returns the N-point Parzen (de la Valle-Poussin) window in
%   a column vector.

% Check for valid window length (i.e., n < 0)
[n,w,trivialwin] = checkOrder(n);
if trivialwin, return, end;

% Index vectors
k = -(n-1)/2:(n-1)/2;
k1 = k(k<-(n-1)/4);
k2 = k(abs(k)<=(n-1)/4);

% Equation 37 of [1]: window defined in three sections
w1 = 2 * (1-abs(k1)/(n/2)).^3;
w2 = 1 - 6*(abs(k2)/(n/2)).^2 + 6*(abs(k2)/(n/2)).^3;
w = [w1 w2 w1(end:-1:1)]';


function [n_out, w, trivalwin] = checkOrder(n_in)
% CHECKORDER Checks the order passed to the window functions.

w = [];
trivalwin = 0;

% Special case of negative orders:
if n_in < 0,
    error('Order cannot be less than zero.');
end

% Check if order is already an integer or empty
% If not, round to nearest integer.
if isempty(n_in) | n_in == floor(n_in),
    n_out = n_in;
else
    n_out = round(n_in);
    warning('Rounding order to nearest integer.');
end

% Special cases:
if isempty(n_out) | n_out == 0,
    w = zeros(0,1);               % Empty matrix: 0-by-1
    trivalwin = 1;
elseif n_out == 1,
    w = 1;
    trivalwin = 1;
end


function lam_adj = eigensp_adj(lam,n,p)
% Eigen spectrum adjustment for EVD on finite samples

r = p/n;

bp = (1+sqrt(r))^2;
bm = (1-sqrt(r))^2;

vv = [bm:(bp-bm)/(5*p-1):bp];

gv = 1./(2*pi*r*vv).*sqrt((vv-bm).*(bp-vv));

for i = 1:length(gv)
    gvd(i) = sum(gv(1:i));
end
gvd = gvd/max(gvd);

for i = 1:p
    i_norm = i/p;
    [minv,minx] = min(abs(i_norm-gvd));
    lam_emp(i) = vv(minx);
end
lam_emp = rot90(lam_emp,2);

lam_adj = lam./lam_emp;


function [sum_dat] = sumN(dat)
% sum of the all elements of the dat matrix

sum_dat = sum(dat(:));


function c = xcorr2(a,b)
% two dimensional cross correlation

if nargin == 1
    b = a;
end

c = conv2(a, rot90(conj(b),2));

function kurt = kurtn(x)
% Normilized kurtosis funtion so that for a Gaussian r.v. the kurtn(g) = 0
% Input: x 1:N vector


x = x-mean(x);

x = x/std(x);

kurt = mean(x.^4)-3;

function [d_mdl_maxmin, r1_mdl_maxmin] = ica_fuse_model_est_cca(X,Y)
% Implementation of order-selection using PCA-CCA. Though this code can be
% used in both the sample-rich and sample-poor regime, it was developed to
% perform order-selection in the sample-poor regime.
%
% Input: 
%       X, Y - Datasets. Note that the number of columns (samples) must be
%              the same for the two datasets.
%
% Output:
%       d_mdl_maxmin - Number of common components shared between the two
%                      datasets
%       r1_mdl_maxmin - Estimated order
%
% If using this technique, please cite the following papers:
%
% [1] Y. Levin-Schwartz, Y. Song, P. J. Schreier, V. D. Calhoun, & Adali, T. "Sample-poor estimation of order and common signal subspace with application to fusion of medical imaging data," NeuroImage, 2016, 134, 486-493
% [2] Y. Song, P. J. Schreier, D. Ramirez, & T. Hasija, "Canonical correlation analysis of high-dimensional data with very small sample support," Signal Processing, Elsevier, 2016, 128, 449-458

M = size(Y,2);
f = [1 1]; % Make the same if the order should be the same for both datasets (default)

if min(size(X,1),size(Y,1))>=M % sample-poor (multimodal fusion)

    rmax = floor(M/2);
    Wxx = 1/M*(X')*X;
    Wyy = 1/M*(Y')*Y;
    [Vx, Lx] = eig(Wxx);
    [Vy, Ly] = eig(Wyy);
    [~, ix] = sort(diag(Lx),'descend');
    [~, iy] = sort(diag(Ly),'descend');
    Vx_sort = Vx(:,ix);
    Vy_sort = Vy(:,iy);
    [d_mdl_maxmin, r1_mdl_maxmin, ~] = IC_maxmin_2(Vx_sort,Vy_sort,M,'mdl','real',f,rmax); %old detector (detector #3)

else % sample-rich (multiset)

    M = size(Y,1);
    rmax = floor(0.5*M);
    Rxx = 1/M*(X)*X';
    Ryy = 1/M*(Y)*Y';
    [Ux, Lx] = eig(Rxx);
    [Uy, Ly] = eig(Ryy);
    [~, xid] = sort(diag(Lx),'descend');
    [~, yid] = sort(diag(Ly),'descend');
    Ux_sort = Ux(:,xid);
    Uy_sort = Uy(:,yid);
    [d_mdl_maxmin, r1_mdl_maxmin, ~] = IC_maxmin_rich(Ux_sort,Uy_sort,X,Y,M,'mdl','real',f,rmax);
end



function [d, r1Est, r2Est] = IC_maxmin_2(Vx,Vy,M,Thre,RealComp,f,rmax)
% IC_maxmin.m - MDL based max-min method
%
% Thre                'MDL' or 'AIC' ('AIC' is used only for test)
% RealComp            'real' or 'comp'
% rmax                {r1,r2} = 1,...,rmax
% 
% 
%
% Input:
%
% Vx:        contains right singular vectors of X
% Vy:        contains right singular vectors of Y
% M:         number of samples
% Thre:     'MDL' ('AIC' is used only for test)
% RealComp: 'real' for real-valued data; 
%           'comp' for complex-valued data
% f:         f(1) is the number of independent signals in X
%            f(2) is the number of independent signals in Y
% rmax       {r1,r2} = 1,...,rmax
%             r1: the rank PCA keeps in X
%             r2: the rank PCA keeps in Y
%
% Output:
%
% d:         number of correlated signal between X and Y
% r1Est:     the optimum rank that PCA should keep for X
% r2Est:     the optimum rank that PCA should keep for Y

for r1 = 1:rmax
    % If f1=f2, we let r1=r2.
    if f(1)==f(2), r2range  = r1; else r2range = 1:rmax; end
    for r2 = r2range
        ga = sort(svd(Vx(:,1:r1)'*Vy(:,1:r2)),'descend');
        for r3 = 0:min(r1,r2)-1
            switch lower(RealComp)
                case 'real'
                    dof = r1*r2 - (r1-r3)*(r2-r3);
                    Loglike = M/2*log(prod(1-ga(1:r3).^2));
                case 'comp'
                    dof = 2*r3*(r1+r2-r3);
                    Loglike = M*log(prod(1-ga(1:r3).^2));
            end
            switch lower(Thre)
                case 'mdl'
                    pen = 1/2*log(M)*dof;
                case 'aic'
                    pen = dof;
            end
            IC(r3+1) = Loglike + pen;
        end
        if min(IC)==-inf, 
            idx(r1,r2)=0;
        else
            idx(r1,r2) = find(IC==min(IC))-1; 
        end
        clear IC
    end
end
d = max(max(idx));

[~, location] = max(idx(:));
[r1Est,r2Est] = ind2sub(size(idx),location);

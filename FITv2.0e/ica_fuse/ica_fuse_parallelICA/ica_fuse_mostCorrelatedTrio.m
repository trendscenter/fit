function [ ind, m12, m23, m13 ] = ica_fuse_mostCorrelatedTrio( mod1load, mod2load, mod3load )
% this function was created to find correlated triplets from ICA or 3wayICA
% modalities
%% Get the correlations between modalities and the indexes
rho12 = ica_fuse_corr(mod1load,mod2load);
rho23 = ica_fuse_corr(mod2load,mod3load);
rho13 = ica_fuse_corr(mod1load,mod3load);
%% size(modXload,2) must be the number of IC on each modality
Ac = ica_fuse_combvec(1:size(mod1load,2),  1:size(mod2load,2),  1:size(mod3load,2))';
Ncomb = size(Ac,1);
c12 = zeros(Ncomb,1);
c13 = zeros(Ncomb,1);
c23 = zeros(Ncomb,1);
for jj = 1:Ncomb
    c12(jj) = (rho12(Ac(jj,1),Ac(jj,2)));
    c23(jj) = (rho23(Ac(jj,2),Ac(jj,3)));
    c13(jj) = (rho13(Ac(jj,1),Ac(jj,3)));
end
Mean_Correlation = mean(abs([c12 c13 c23]),2);
[ddd, t] = sort(Mean_Correlation,'descend');
m12 = c12(t);
m13 = c13(t);
m23 = c23(t);
ind =  Ac(t,:);



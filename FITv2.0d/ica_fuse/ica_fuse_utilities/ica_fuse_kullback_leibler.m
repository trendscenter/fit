function [kl_metric] = ica_fuse_kullback_leibler(d1, d2)
% Calculate KL metric

d1 = d1(:);
d2 = d2(:);

% Get in
ind = find((d1~=0) & (d2~=0));
d1 = d1(ind);
d2 = d2(ind);

d1 = d1/sum(d1);
d2 = d2/sum(d2);

%kl_metric = sumN(d1.*log(d1./(d2)));
k = d1.*log(d1./(d2));

kl_metric = sum(k(:));

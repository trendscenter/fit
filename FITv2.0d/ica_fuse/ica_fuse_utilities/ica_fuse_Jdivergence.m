function j_div = ica_fuse_Jdivergence(d1, d2)
% Calculate KL metric

d1 = d1(:);
d2 = d2(:);

% Get in
ind = find((d1~=0) & (d2~=0));
d1 = d1(ind);
d2 = d2(ind);

d1 = d1/sum(d1);
d2 = d2/sum(d2);

k1 = d1.*log(d1./(d2));
k2 = d2.*log(d2./(d1));

k = 0.5*(k1 + k2);

% J divergence
j_div = sum(k(:));

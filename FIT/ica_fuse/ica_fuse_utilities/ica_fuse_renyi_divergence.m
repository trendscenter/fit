function renyi_divergence = ica_fuse_renyi_divergence(d1, d2, alpha)
% Calculate KL metric

d1 = d1(:);
d2 = d2(:);

% Get in
ind = find((d1~=0) & (d2~=0));

d1 = d1(ind);
d2 = d2(ind);

d1 = d1/sum(d1);
d2 = d2/sum(d2);

%renyi_divergence = (1/alpha)*log(sum((d1.^alpha).*(d2.^(1-alpha))));
renyi_divergence = (1/(alpha-1))*log(sum((d1.^alpha).*(d2.^(1-alpha))));
function alpha_div = ica_fuse_alpha_divergence(d1, d2, alpha)
% Calculate KL metric

d1 = d1(:);
d2 = d2(:);

% Get in
ind = find((d1~=0) & (d2~=0));
d1 = d1(ind);
d2 = d2(ind);

d1 = d1/sum(d1);
d2 = d2/sum(d2);

%k = (alpha.*d1) + ((1 - alpha).*d2) - ((d1.^alpha) .* (d2.^(1 - alpha)));

k = (alpha.*d1) + ((1 - alpha).*d2) - (((d1 ./ d2).^alpha).*d2);

k = k / (alpha*(1 - alpha));

alpha_div = sum(k(:));

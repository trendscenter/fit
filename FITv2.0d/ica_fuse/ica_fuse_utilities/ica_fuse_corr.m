function c = ica_fuse_corr(x, y)
% Linear or rank correlation

if ~exist('y', 'var')
    y = x;
    [ny1, ny2] = size(y);
else
    ny2 = size(y, 2);
end

% Check dimensions of x and y
if (size(x, 1) ~= size(y, 1))
   error('X and Y must have the same number of rows.'); 
end

[nx1, nx2] = size(x);

% Initialise pairwise correlation
c = zeros(nx2, ny2);

% Loop over rows
for ii = 1:nx2
    % Loop over cols
    for jj = 1:ny2
        c(ii, jj) = ica_fuse_corr2(x(:, ii), y(:, jj));
    end
    % End loop over cols
end
% End loop over rows
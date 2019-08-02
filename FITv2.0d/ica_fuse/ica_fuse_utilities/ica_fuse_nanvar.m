function y = ica_fuse_nanvar(x)
% Ignore Nan's while calculating variance. Currently works for 2 dimensions

y = ica_fuse_nanstd(x);
y = y.^2;

% 
% if (numel(x) == length(x))
%     % If x is a vector
%     x(isnan(x)) = [];
%     if isempty(x)
%         y = NaN;
%         return;
%     end
%     y = var(x);
% else
%     % Initialise variance to number of columns
%     y = repmat(NaN, 1, size(x, 2));
%     % Loop over columns
%     for yDim = 1:size(x, 2)
%         temp = x(:, yDim);
%         temp(isnan(temp)) = [];
%         if ~isempty(temp)
%             y(yDim) = var(temp);
%         end
%     end
%     % End loop over columns
% end
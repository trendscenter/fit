function ind = ica_fuse_good_cells( mycell )
%% Find good cells

if ~iscell(mycell)
    mycell = {mycell};
end

ind = cellfun('isempty', mycell);

% Good cells
ind = (ind == 0);

ind = ind(:)';
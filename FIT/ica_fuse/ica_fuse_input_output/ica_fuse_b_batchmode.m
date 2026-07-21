%ICA_FUSE_B_BATCHMODE Determine whether GIFT/FIT is executing in batch mode.
%
% Returns:
%   b_batch_true - Logical true if execution originated from one of the
%                  known batch-mode entry-point functions; otherwise false (e.g., run fom GUI).
%
% Method:
%   Batch mode is detected by searching the MATLAB call stack for known
%   batch entry-point functions. This helper is useful when the caller does
%   not have direct access to the analysis configuration.

% Cyrus Erik Eierud, TReNDS 7/21/26

b_batch_true = false;

struStack = dbstack;

for i = 1:numel(struStack)

    if strcmpi(struStack(i).file, 'ica_fuse_batch_file.m') || ...
       strcmpi(struStack(i).file, 'ica_fuse_anyway_fusion_batch') || ...
       strcmpi(struStack(i).file, 'ica_fuse_batch_file_dfuse_final.m')

        b_batch_true = true;
        return;
    end

end


function groups_icasig = ica_fuse_backReconstruct_groups(A, data, numGroups, numSubjects)
%% Back reconstruct components for each group
%
% Inputs:
% 1. A - Mixing matrix
% 2. data - Original data
% 3. numGroups - Number of groups
% 4. numSubjects - Number of subjects per group
%
% Outputs:
% groups_icasig - Back-reconstructed group components
%

%% Do group back reconstruction if number of groups is greater than 1
if numGroups > 1
    fprintf('\n');

    disp('Doing group back reconstruction on components ...');

    % Number of voxels & components
    voxels = size(data, 2);
    numComp = size(A, 2);

    %% Initialise groups icasig variable
    groups_icasig = zeros(voxels, numComp, numGroups);

    %% Loop over groups
    for nGroups = 1:numGroups
        groupInd = ica_fuse_get_groupInd(nGroups, numSubjects);
        % get the back reconstruted component
        groups_icasig(:, :, nGroups) = (pinv(A(groupInd, :))*data(groupInd, :))';
    end
    %% End loop over groups


    disp('Done with group back reconstruction on components');

else
    groups_icasig = [];
end
function ica_fuse_batch_file_dfuse_final(inputFile)
% Cyrus Eierud 070126
% Final run after all states have run batch file for dynamic data fusion
% by submitting the first batch file (inputFile) again, but to this function
% After all states are processed this function
% will match all the components across the states
%
% Dynamic FNC Fusion Workflow
% ---------------------------
%
% This example performs fusion ICA between dynamic functional network
% connectivity (dFNC) derived from fMRI and a structural modality (e.g.,
% gray matter volume, GMV) from the same subjects.
%
% Prerequisites
% -------------
% 1. Run the fMRI analysis in GIFT using a NeuroMark template
%    (e.g., Neuromark_fMRI_2.2_modelorder-multi.nii).
%
% 2. Run the dynamic FNC analysis in GIFT. This will generate:
%      - n_comps independent components (e.g., 105)
%      - n_states dynamic states (e.g., 3)
%      - a GIFT output directory (e.g., /location/I/saved/GIFT)
%      - a dFNC post-processing file:
%          {prefix}_dfnc_post_process.mat
%          ...
%
% 3. Prepare a second imaging modality for the same subjects. E.g.,
%    the second modality may be gray matter volume (GMV), obtained by
%    normalizing and smoothing each subject's T1-weighted MRI. Each subject
%    should have one corresponding structural data file.
%
% Running the Dynamic FNC Fusion
% ------------------------------
%
% 1. Copy the template batch file
%
%      FIT/ica_fuse_batch_files/input_fuse_dfnc_state1.m
%
%    to your own working directory, for example
%
%      /users/myfiles/dfnc_code/
%
% 2. Edit input_fuse_dfnc_state1.m and update at least the following
%    variables:
%
%      outputDir='/location/to/save/fusionresults/';
%
%      dynfitFileGiftDfncPostProcessMat = ...
%          '/location/I/saved/GIFT/{prefix}_dfnc_post_process.mat';
%
%      prefix='prefix_dFNC_state1_GMV_5comp';
%
%      featureNames={'dFNC','GMV'};
%
%      modality={'fmri','smri'};
%
%      files={ ...
%          [outputDir '/dfit/*state-01_avgConn.mat', ...
%          '/location/I/saved/GMV/*_sMRI.mat'};
%      % Asterisks (*) above denotes subject IDs, which has to match across modalities.
%      % Make sure the dir before /dfit/ is the outputDir (from above).
%
%      numComp=5;
%
% 3. Copy input_fuse_dfnc_state1.m to
%
%      input_fuse_dfnc_state2.m
%
% 4. Edit input_fuse_dfnc_state2.m as follows:
%
%      - Remove the line containing
%            dynfitFileGiftDfncPostProcessMat
%
%      - Change
%            prefix='prefix_dFNC_state1_GMV_5comp'
%        to
%            prefix='prefix_dFNC_state2_GMV_5comp'
%
%      - Change the dFNC file (related with files param) from
%            *state-01_avgConn.mat
%        to
%            *state-02_avgConn.mat
%
% 5. Repeat the previous step for all remaining dynamic states produced by
%    GIFT. For each state:
%
%      - Copy the previous batch file (input_fuse_dfnc_state2.m).
%      - Update the prefix to match the state number.
%      - Update the dFNC input file to the corresponding
%        *state-XX_avgConn.mat file.
%
% 6. Run the FIT toolbox for the first state:
%
%      ica_fuse_batch_file( ...
%          '/users/myfiles/dfnc_code/input_fuse_dfnc_state1.m');
%
%    During this first run, FIT generates processed files (under dfit
%    folder), which are reused by all subsequent states.
%
% 7. Run FIT for each remaining state:
%
%      ica_fuse_batch_file( ...
%          '/users/myfiles/dfnc_code/input_fuse_dfnc_state2.m');
%
%      ...
%
%      ica_fuse_batch_file( ...
%          '/users/myfiles/dfnc_code/input_fuse_dfnc_stateN.m');
%
% 8. After all states have completed, perform the final post-processing
%    step to sort the fusion components consistently across states:
%
%      ica_fuse_batch_file_dfuse_final( ...
%          '/users/myfiles/dfnc_code/input_fuse_dfnc_state1.m');
%
% Results
% -------
%
% The fusion results are written to files of the form
%
%      {prefix}_stateX_*_Ycomp_joint_comp_ica_feature_1_00Z.asc
%
% and
%
%      {prefix}_stateX_*_Ycomp_joint_comp_ica_feature_2_00Z.asc
%
% where
%
%      X = dynamic state number (1...N)
%      Y = number of fusion components (numComp)
%      Z = component number (1...Y)
%
% The final component matching and state ordering are stored in
%
%      dynamicFusion_postprocessing_component_matches_*.mat

[pathstr, fName, extn] = fileparts(inputFile);
if isempty(pathstr)
    pathstr = fileparts(which(inputFile));
end

addpath(pathstr); %%%%%%%%%%%% needed?

inputFile = fullfile(pathstr, [fName, extn]);

eval(fName);

% Option for dynamic FIT preprocessing
oc_dfit = ica_fuse_cls_dfit(inputFile, outputDir); % Needed to check if used wants dfit
b_dfit_selected_in_batch_or_gui = oc_dfit.get_b_dfit_selected_in_batch_or_gui();
if b_dfit_selected_in_batch_or_gui
    n_ret = oc_dfit.dyn_matching_components_across_states(); %dfit selected and engaged
else
    error('Missing information for ica_fuse_batch_file_dfuse_final');
end

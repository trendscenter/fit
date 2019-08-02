function standardize_sub = ica_fuse_flagStandardizeSub(selectedModalities)
% Return flag for standarizing subjects to z-scores or not

ica_fuse_defaults;
global STANDARDIZE_SUBJECTS;

global FLAG_STANDARDIZE_SUB;

if ~isstruct(STANDARDIZE_SUBJECTS)
    error('STANDARDIZE_SUBJECTS variable in ica_fuse_defaults.m must be defined as a structure');
end

% find if specific modality exists or not
eeg_match = strmatch('eeg', lower(selectedModalities), 'exact');
fmri_match = strmatch('fmri', lower(selectedModalities), 'exact');
smri_match = strmatch('smri', lower(selectedModalities), 'exact');

% Initialise standardize subjects variable
standardize_sub = 'no';

if (length(eeg_match) > 0) && (length(fmri_match) == 0) && (length(smri_match) == 0)
    % eeg-eeg fusion
    standardize_sub = STANDARDIZE_SUBJECTS.eeg_eeg;
elseif (length(eeg_match) == 0) && (length(fmri_match) > 0) && (length(smri_match) == 0)
    % fmri-fmri fusion
    standardize_sub = STANDARDIZE_SUBJECTS.fmri_fmri;
elseif (length(eeg_match) == 0) && (length(fmri_match) == 0) && (length(smri_match) > 0)
    % smri-smri fusion
    standardize_sub = STANDARDIZE_SUBJECTS.smri_smri;
elseif (length(eeg_match) > 0) && (length(fmri_match) > 0)
    % eeg-fmri fusion
    standardize_sub = STANDARDIZE_SUBJECTS.eeg_fmri;
elseif (length(eeg_match) > 0) && (length(smri_match) > 0) && (length(fmri_match) == 0)
    % eeg-smri fusion
    standardize_sub = STANDARDIZE_SUBJECTS.eeg_smri;
elseif (length(eeg_match) == 0) && (length(smri_match) > 0) && (length(fmri_match) > 0)
    % fmri-smri fusion
    standardize_sub = STANDARDIZE_SUBJECTS.fmri_smri;
end

FLAG_STANDARDIZE_SUB = standardize_sub;

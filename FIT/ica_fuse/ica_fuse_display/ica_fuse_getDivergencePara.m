function [div_name, div_value] = ica_fuse_getDivergencePara

ica_fuse_defaults;
global DIVERGENCE_PARAMETERS;

if isempty(DIVERGENCE_PARAMETERS)
    error(['Check the global variable DIVERGENCE_PARAMETERS in ica_fuse_defaults.m']);
end


% Check divergence name from defaults and if no match is found by default select first name 
% Get all the divergence names
[divNames] = ica_fuse_divergence;

matchIndex = strmatch(lower(DIVERGENCE_PARAMETERS{1}), divNames, 'exact');
if isempty(matchIndex)
    div_name = divNames{1};
else
    div_name = DIVERGENCE_PARAMETERS{1};
end

if strcmpi(div_name, 'renyi') | strcmpi(div_name, 'alpha')
    div_value = DIVERGENCE_PARAMETERS{2};
else
    div_value = [];
end


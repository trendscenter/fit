function [A, jicasig, wM] = ica_fuse_mcca_reference(data, numPC, reference, opts)
%% MCCA with reference
%
% Inputs:
% 1. data - Cell array of length equal to number of modalities
% 2. numPC - Number of principal components
% 3. reference - Reference vector
% 4. opts - Options
%
% Outputs:
% 1. A - Mixing coefficients
% 2. jicasig - Canonical scores
%

global MCCAR_LAMBDA;

if (~iscell(data))
    % assuming last dimension is modalities
    data = squeeze(num2cell(data, [1, 2]));
end

lam = MCCAR_LAMBDA;

if (isempty(lam))
    lam = 0.98;
end


numIC = min(numPC);
try
    numIC = opts.numc;
catch
end

V = cell(1, length(data));
Lambda = V;
whitesig = V;
whiteM = V;
dewhiteM = V;
A = V;
ccacomp  = V;

for nD = 1:length(data)
    [V{nD}, Lambda{nD}, whitesig{nD}, whiteM{nD}, dewhiteM{nD}] = ica_fuse_calculate_pca(data{nD}, numPC(1));
end

if (numel(reference) == length(reference))
    reference = reference(:);
end

whitesig = cat(3, whitesig{:});

% MCCA with reference
W = ica_fuse_mcca_ssqcor_R(whitesig, numIC, reference, lam);
wM = cell(1, length(data));
for nM = 1:length(data)
    A{nM} = (W(:, :, nM)*whitesig(:, :, nM))';
    ccacomp{nM} = pinv(A{nM})*data{nM}';
    wM{nM} = pinv(A{nM});
end

jicasig = cat(2, ccacomp{:});
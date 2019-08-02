function [divNames, div_value] = ica_fuse_divergence(out1, out2, divCriteria, value)
% Divergence

div_value = [];
divNames = {'kl', 'j', 'renyi', 'alpha'};    
if nargin ==0
    return;
end

switch lower(divCriteria)
    case 'kl'   
        % Kullback and Leibler divergence
        [div_value] = ica_fuse_kullback_leibler(out1, out2);
    case 'j'
        % J divergence
        [div_value] = ica_fuse_Jdivergence(out1, out2);
    case 'renyi'
        % Renyi divergence
        [div_value] = ica_fuse_renyi_divergence(out1, out2, value);
    case 'alpha'
        % Alpha divergence
        [div_value] = ica_fuse_alpha_divergence(out1, out2, value);
    otherwise
        error('Unknown divergence criteria passed');
end
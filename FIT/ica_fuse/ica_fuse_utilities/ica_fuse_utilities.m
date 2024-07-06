function ica_fuse_utilities(selectedStr)
%% Utilities
% Current utilities are optimal features and histogram plot

switch lower(selectedStr)
    
    case 'optimal features' 
        % Optimal features (Rank the features)
        ica_fuse_optimal_features;        
        
    case 'histogram plot'
        % Histogram plot
        ica_fuse_plotHistogram;
        
    case 'write talairach table'
        % Write talairach table
        ica_fuse_talairach;
        
    otherwise 
        
        
end
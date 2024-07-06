function ica_fuse_toolboxes(selectedStr)
% List of toolboxes

switch lower(selectedStr)
    
    case 'parallel ica' 
        
        ica_fuse_parallel_ICAToolbox;
        
    otherwise 
        error('Unrecognized toolbox');
        
end
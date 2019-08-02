% estimate intersect between two columns pairwisely

function y = ica_fuse_intersect_pairwise(x,range)

clear y
for j1 = range
    for j2 = range
        y(j1,j2) = length(intersect(x(:,j1),x(:,j2)))/size(x,1);
    end
end

    


        
    
    

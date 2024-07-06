%% cluster candidate order


function [s_cluster,summary] = ica_fuse_snp_est_cluster(sind_candid,step)

x1 = diff(sind_candid);
sind1 = find(x1 > step);

clear s_cluster summary;
if isempty(sind1)
    s_cluster(1).ind = sind_candid;
else
    s_cluster(1).ind = sind_candid(1:sind1(1));
    if length(sind1) == 1
        s_cluster(2).ind = sind_candid(sind1(end)+1:end);
    elseif length(sind1) > 1
        for j = 2:length(sind1)
            s_cluster(j).ind = sind_candid(sind1(j-1)+1:sind1(j));
        end
        s_cluster(j+1).ind = sind_candid(sind1(end)+1:end);
    end
end


for j = 1:length(s_cluster)
    sind1 = s_cluster(j).ind;
    if length(sind1)==1
        s_cluster(j).mq = 1e-3;
    else
        s_cluster(j).mq = length(sind1)/(sind1(end)-sind1(1)+1);
    end
    summary(j,:) = [length(sind1),s_cluster(j).mq];
end










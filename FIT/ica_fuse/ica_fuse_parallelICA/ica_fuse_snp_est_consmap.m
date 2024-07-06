%% consistency map construction


function [consmap_s,consmap_a,indmap_s,indmap_a] = ica_fuse_snp_est_consmap(filepath,s_ica,ncomp_snp_min,ncomp_snp_max)


% consistency map
clear consmap_s consmap_a indmap_s indmap_a;
for j1 = 1:ncomp_snp_min
    indmap_s(1:j1,j1) = 1:j1;
    indmap_a(1:j1,j1) = 1:j1;
end
for j1 = ncomp_snp_min:ncomp_snp_max-1
    %=== component ====================================
    icasig1 = s_ica(j1).icasig(indmap_s(1:j1,j1),:)';
    icasig2 = s_ica(j1+1).icasig';
    [r1,p1] = corr(icasig1,icasig2);
    numic_test = 1:j1;
    numic_match = 1:j1+1;
    while ~isempty(numic_test)
        [b1,b2] = max(max(abs(r1)));   % column index
        [b3,b4] = max(abs(r1(:,b2)));  % row index
        consmap_s(b4,j1) = b3;
        indmap_s(b4,j1+1) = b2;
        r1(:,b2) = 0;
        r1(b4,:) = 0;
        numic_test = setdiff(numic_test,b4);
        numic_match = setdiff(numic_match,b2);
    end
    indmap_s(j1+1,j1+1) = numic_match;
    %=== loading ===================================
    icasig1 = s_ica(j1).loading(:,indmap_a(1:j1,j1));
    icasig2 = s_ica(j1+1).loading;
    [r1,p1] = corr(icasig1,icasig2);
    numic_test = 1:j1;
    numic_match = 1:j1+1;
    while ~isempty(numic_test)
        [b1,b2] = max(max(abs(r1)));   % column index
        [b3,b4] = max(abs(r1(:,b2)));  % row index
        consmap_a(b4,j1) = b3;
        indmap_a(b4,j1+1) = b2;
        r1(:,b2) = 0;
        r1(b4,:) = 0;
        numic_test = setdiff(numic_test,b4);
        numic_match = setdiff(numic_match,b2);
    end
    indmap_a(j1+1,j1+1) = numic_match;
end


% save data
%save([filepath, 'consmap.mat'], 'consmap_s','indmap_s', 'consmap_a', 'indmap_a');

return;








%% reference selection


function [icasig_r,ind_testcomp] = ica_fuse_snp_est_ref(filepath,s_ica,consmap_s,indmap_s,consmap_a,indmap_a,winsize,testN,ncomp_snp_min,ncomp_snp_max,type,pheno)
switch type
    case 'pheno'
        %=== get association map =============
        gd = [];
        for j = ncomp_snp_min:ncomp_snp_max
            [r1,p1] = corr(s_ica(j).loading,pheno);
            [gd(j,1),gd(j,2)] = max(max(abs(r1),[],2));
        end
        figure; plot(gd(:,1));
        saveas(gcf,[filepath, 'gd.fig']); title('max group difference (correlation) as a function of component number');
        clear icasig_r;
        [b1,b2] = sort(gd,'descend');
        x1 = b2(1:winsize);
        med1 = median(gd(x1,1));
        med2 = median(gd(:,1));
        sigma1 = std(gd(x1,1));
        lambda = 0.9+0.1*(med1-med2)/med1;
        gd_th = lambda*(med1);
        sind1 = find(gd(:,1) > gd_th,testN,'first');
        for j = 1:testN
            ind_testcomp(j,1) = sind1(j);
            ind_testcomp(j,2) = gd(ind_testcomp(j,1),2);
            icasig_r(:,j) = s_ica(ind_testcomp(j,1)).icasig(ind_testcomp(j,2),:)';
        end
        save([filepath,'ref_pheno.mat'],'gd','ind_testcomp','icasig_r');
    case 'blind'
        clear max_avg_cons ind_testcomp icasig_r;
        x1 = cumsum(consmap_s,2);
        x2 = [];
        x2(:,ncomp_snp_min:ncomp_snp_max-1-winsize) = x1(:,ncomp_snp_min+winsize:ncomp_snp_max-1) - x1(:,ncomp_snp_min:ncomp_snp_max-1-winsize);    
        x2 = x2/(winsize+1);
        for j1 = 1:size(x2,2)
            x2(j1+1:size(x2,1),j1) = 0;
        end            
        [max_avg_cons,sind1] = max(x2,[],2);
        [b1,b2] = sort(max_avg_cons,'descend');
        for j1 = 1:testN
            ind_testcomp(j1,1) = sind1(b2(j1)) + floor(winsize/2);
            ind_testcomp(j1,2) = indmap_s(b2(j1),ind_testcomp(j1,1));
        end
        for j1 = 1:testN
            icasig_r(:,j1) = s_ica(ind_testcomp(j1,1)).icasig(ind_testcomp(j1,2),:)';
        end
        save([filepath,'ref_blind.mat'],'max_avg_cons','ind_testcomp','icasig_r');
end

close all;
return;
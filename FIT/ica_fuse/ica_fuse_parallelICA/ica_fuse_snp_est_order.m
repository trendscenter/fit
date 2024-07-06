%% stability evaluation based on selected reference snp component and component number estimation


function [ncomp_snp_est,ncomp_snp_candid] = ica_fuse_snp_est_order(filepath,s_ica,consmap_s,consmap_a,icasig_r,ind_testcomp,winsize,testN,ncomp_snp_min,ncomp_snp_max,type1,type2)

clear ncomp_snp_candid_s ncomp_snp_candid_a;
step = floor(winsize/2);
range_comN = ncomp_snp_min+step:ncomp_snp_max-winsize+step;
if strcmp(type1,'row')
    clear loading_r;
    for j1 = 1:testN
        loading_r(:,j1) = s_ica(ind_testcomp(j1,1)).loading(:,ind_testcomp(j1,2))';
    end
end

%% order selection
switch type1
    case 'column'
        th = 0.95;
        %=== component ==========================================
        avg_cons_s = sum(consmap_s)./(1:ncomp_snp_max-1);
        avg_cons_s(1) = avg_cons_s(2);
        x1 = cumsum(avg_cons_s);
        x2 = [];
        x2(1:ncomp_snp_max-1-winsize) = x1(1+winsize:ncomp_snp_max-1)-x1(1:ncomp_snp_max-1-winsize);
        [b1,b2] = sort(x2,'descend');
        med1 = median(b1(1:winsize));
        med2 = median(b1);
        sigma1 = std(b1(1:winsize));
        lambda = th-(1-th)*(med1-med2)/med1;
        cons_th = lambda*(med1);
        ncomp_snp_candid_s = ceil(median(find(x2 < cons_th,testN))) + floor(winsize/2) - 1;
        %=== loading =========================================
        avg_cons_a = sum(consmap_a)./(1:ncomp_snp_max-1);
        avg_cons_a(1) = avg_cons_a(2);
        x1 = cumsum(avg_cons_a);
        x2 = [];
        x2(1:ncomp_snp_max-1-winsize) = x1(1+winsize:ncomp_snp_max-1)-x1(1:ncomp_snp_max-1-winsize);
        [b1,b2] = sort(x2,'descend');
        med1 = median(b1(1:winsize));
        med2 = median(b1);
        sigma1 = std(b1(1:winsize));
        lambda = th-(1-th)*(med1-med2)/med1;
        cons_th = lambda*(med1);
        ncomp_snp_candid_a = ceil(median(find(x2 < cons_th,testN))) + floor(winsize/2) - 1;
        %=== final estimation ===============================
        avg_cons = [avg_cons_s',avg_cons_a'];
        ncomp_snp_candid = [ncomp_snp_candid_s',ncomp_snp_candid_a'];
        ncomp_snp_est = ceil(median([ncomp_snp_candid_s,ncomp_snp_candid_a]));
        %save([filepath,'compest_column.mat'], 'avg_cons','ncomp_snp_candid', 'ncomp_snp_est');
    case 'row'
        th = 0.90;
        %=== component ===========================================
        clear icasig_s indcomp_s loading_a indcomp_a;
        for j = ncomp_snp_min:ncomp_snp_max
            icasig1 = s_ica(j).icasig;
            [r1,p1] = corr(icasig1',icasig_r);
            [b1,b2] = max(abs(r1));
            for k = 1:testN
                icasig_s{k}(:,j) = icasig1(b2(k),:)';
                indcomp_s(j,k) = b2(k);
            end
            loading1 = s_ica(j).loading;
            [r1,p1] = corr(loading1,loading_r);
            [b1,b2] = max(abs(r1));
            for k = 1:testN
                loading_a{k}(:,j) = loading1(:,b2(k));
                indcomp_a(j,k) = b2(k);
            end
        end
        %=== check cross correlation =============================
        %         figure;
        %         for k = 1:testN
        %             [r1,p1] = corr(icasig_s{k});
        %             subplot(testN,1,k); imagesc(abs(r1)); colorbar; title(['correlation map of component ', num2str(k)]);
        %         end
        %         saveas(gcf,[filepath,'crosscorr_', type1, '_', type2, '.fig']);
        %=== component ===============================================
        clear cons_s s_overlap;
        cons_s(1,:,1:testN) = 0;
        for k = 1:testN
            icasig1 = icasig_s{k};
            [r1,p1] = corr(icasig1);
            % consistency of most contributing SNPs
            N = ceil(0.05*size(icasig1,1));
            for j = ncomp_snp_min:ncomp_snp_max-winsize+1
                clear sind1;
                clear summary;
                for j1 = j:j+winsize-1
                    [b1,b2] = sort(abs(zscore(icasig1(:,j1))),'descend');
                    sind1(:,j1) = b2(1:N);
                end
                y = intersect_pairwise(sind1,j:j+winsize-1);
                summary(:,:,1) = y(j:j+winsize-1,j:j+winsize-1);
                summary(:,:,2) = abs(corr(icasig1(:,j:j+winsize-1)));
                s_overlap(j,k).summary = summary;
                y1 = summary(:,:,1);
                y1 = y1(:);
                y1(y1 == 1) = [];
                y2 = summary(:,:,2);
                y2 = y2(:);
                y2(y2 == 1) = [];
                cons_s(j,1,k) = mean(abs(y1));
                cons_s(j,2,k) = mean(abs(y2));
            end
        end
        %=== loading ==========================================
        clear cons_a;
        cons_a(1,1:testN) = 0;
        for k = 1:testN
            loading1 = loading_a{k};
            for j = ncomp_snp_min:ncomp_snp_max-winsize+1
                r1 = abs(corr(loading1(:,j:j+winsize-1)));
                r1 = r1(:);
                r1(r1 == 1) = [];
                cons_a(j,k) = mean(abs(r1));
            end
        end
        %=== estimation =======================================
        %=== component ===========================================
        ind_feature = 1;
        for j1 = 1:testN
            [b1,b2] = sort(cons_s(:,ind_feature,j1),1,'descend');
            x1 = b2(1:winsize);
            med1 = median(cons_s(x1,ind_feature,j1));
            med2 = median(cons_s(:,ind_feature,j1));
            sigma1 = std(cons_s(x1,ind_feature,j1));
            lambda = th+(1-th)*(med1-med2)/med1;
            cons_th = lambda*(med1);
            sind_l = find(cons_s(:,ind_feature,j1) > cons_th,1,'first');
            sind_r = find(cons_s(:,ind_feature,j1) > cons_th,1,'last');
            sind1 = sind_l:sind_r;
            sind_candid = sind1(find(cons_s(sind1,ind_feature,j1) > cons_th));
            [s_cluster,summary] = ica_fuse_snp_est_cluster(sind_candid,3);
            sind1 = find(summary(:,1) >= winsize);
            if ~isempty(sind1)
                [cons_max,sind_max] = max(summary(sind1,2));
                for j2 = 1:sind_max
                    if summary(sind1(j2),2) > lambda*cons_max
                        cluster_ind = sind1(j2);
                        break;
                    end
                end
                ncomp_snp_candid_s(j1,1) = s_cluster(cluster_ind).ind(1) + winsize - 1;
            else
                [l_max,sind_max] = max(summary(:,1));
                cluster_ind = sind_max;
                ncomp_snp_candid_s(j1,1) = ceil(median(s_cluster(cluster_ind).ind)) + floor(winsize/2) - 1;
            end
        end
        %=== loading ==========================================
        for j1 = 1:testN
            [b1,b2] = sort(cons_a(:,j1),1,'descend');
            x1 = b2(1:winsize);
            med1 = median(cons_a(x1,j1));
            med2 = median(cons_a(:,j1));
            sigma1 = std(cons_a(x1,j1));
            lambda = th+(1-th)*(med1-med2)/med1;
            cons_th = lambda*(med1);
            sind_l = find(cons_a(:,j1) > cons_th,1,'first');
            sind_r = find(cons_a(:,j1) > cons_th,1,'last');
            sind1 = sind_l:sind_r;
            sind_candid = sind1(find(cons_a(sind1,j1) > cons_th));
            [s_cluster,summary] = ica_fuse_snp_est_cluster(sind_candid,3);
            sind1 = find(summary(:,1) >= winsize);
            if ~isempty(sind1)
                [cons_max,sind_max] = max(summary(sind1,2));
                for j2 = 1:sind_max
                    if summary(sind1(j2),2) > lambda*cons_max
                        cluster_ind = sind1(j2);
                        break;
                    end
                end
                ncomp_snp_candid_a(j1,1) = s_cluster(cluster_ind).ind(1) + winsize - 1;
            else
                [l_max,sind_max] = max(summary(:,1));
                cluster_ind = sind_max;
                ncomp_snp_candid_a(j1,1) = ceil(median(s_cluster(cluster_ind).ind)) + floor(winsize/2) - 1;
            end
        end
        %=== final estimation =================================
        ncomp_snp_candid = [ncomp_snp_candid_s,ncomp_snp_candid_a];
        x1 = ceil(median([min(ncomp_snp_candid_s),max(ncomp_snp_candid_s)]));
        x2 = ceil(median([min(ncomp_snp_candid_a),max(ncomp_snp_candid_a)]));
        ncomp_snp_est = ceil(median([x1,x2]));
       % save([filepath, 'compest_', type1, '_', type2, '.mat'], 's_overlap','cons_s','cons_a', 'ncomp_snp_candid','ncomp_snp_est', 'icasig_s','indcomp_s','loading_a','indcomp_a');
end
return;




function visualize_spatial_loading_afterFDR(s_batch)

    %%Parallel Multilink Group Joint ICA
    
    %%%% Parameter setup
    stru_in = ica_fuse_eval_script(s_batch);

    if ~exist( stru_in.source, 'dir' )
        mkdir(stru_in.source);
    end
    
    num_subjects = stru_in.num_subjects;
    hc_sub = stru_in.hc_sub;
    ad_sub = stru_in.ad_sub;
    num_fmri = stru_in.num_fmri;

    fig_position = [50 50 1000 800];
    savenifti=1;
    mont_indx2=1;
    sub_l=[];


    xtck=1;
    ld=1;

    stru=load(stru_in.param_pmlg);
    mask_indx=stru.sesInfo.mask_ind;
    V=stru.sesInfo.HInfo.V;
    dim=V.dim(1:3);
    maskindxgm=stru.sesInfo.userInput.maskmod1;
    maskindxicn=stru.sesInfo.userInput.maskmod2;

    %%%%%%%%%%%%%% Finding most significant p-values to compare results  %%%%%%
    %%%%%%%%%%%%%% FDR 

    file_ttest=[stru_in.source, '/','spatial_Groupttest_',stru_in.procedure,'_fdr.pdf'];
    for j=1:stru_in.num_fmri    
            %%%% 
            sub_file=[stru.sesInfo.outputDir, [filesep stru.sesInfo.userInput.prefix '_ica_br'], num2str(j), '.mat'];

            [ic_gmcm, ic_gm, ic_icn, tc_gmcm, tc_gm, tc_icn]  = fun_inner(sub_file);
            %stru_iter = load(sub_file);
            recon_maps=ic_gmcm;  %%%20x50485
            recon_maps_gm=recon_maps(:, 1:length(maskindxgm));
            recon_maps_icn=recon_maps(:, length(maskindxgm)+1:end);
            tc_=tc_gmcm;
            clear recon_maps;
            clear ic_gmcm tc_gmcm;

            sub_hc=tc_(1:stru_in.hc_sub, :);
            sub_ad=tc_(stru_in.hc_sub+1:end, :);

            mn_sub_hc=mean(sub_hc);
            mn_sub_ad=mean(sub_ad);

            [h,p,ci,stats] = ttest2(sub_hc, sub_ad);
            %p
            %%p_fdr=mafdr(p)
            %%p_masked = p_fdr < 0.05
            [pID, ~, p_masked] = ica_fuse_fdr(p, 0.05); %ce042525 0.99
            %p_masked            
            log_p=-log10(p).*sign(stats.tstat);
            log_p(~p_masked)=0;
            p(~p_masked)=0;

            %%sum(p_masked)
            comp_num=find(abs(log_p)>3.5);  %%% p-value> 0.0003 (>3.5) %ce042525 0.01

            if comp_num
                sel_logp=log_p(abs(log_p)>3.0);  %%% p-value> 0.0008 (>3.0)%ce042525 0.01
                p_val1=p(abs(log_p)>3.0);
                p_val=[stru_in.p_val, p_val1];
                sub_l=[sub_l, sel_logp]
                %%%%%%%%%%%%%%%%%%%%%%%%%
                for k=comp_num
                    mn_hc_l(ld)=mn_sub_hc(k);
                    mn_ad_l(ld)=mn_sub_ad(k);
                    ld=ld+1;

                    xtckl{xtck}=['Com-', num2str(j), '_C-', num2str(k)];
                    xtck=xtck+1;

                    recon_mapgm=recon_maps_gm(k, :);
                    recon_mapgm_img=zeros(dim);
                    recon_mapgm_img(maskindxgm)=recon_mapgm;        

                    recon_mapicn=recon_maps_icn(k, :);
                    recon_mapicn_img=zeros(dim); %ce042425
                    recon_mapicn_img(maskindxicn)=recon_mapicn; 

                    if (exist(stru_in.source, 'dir') ~= 7)
                        mkdir(stru_in.source)
                    end

                    V.fname=[stru_in.source, '/', 'subgmtmp22', '.nii'];
                    ica_fuse_spm_write_vol(V, recon_mapgm_img);                
                    V.fname=[stru_in.source, '/', 'subicntmp22', '.nii'];
                    ica_fuse_spm_write_vol(V, recon_mapicn_img);  


                    if savenifti
                        % apply display parameters
                        StrucDim=[]; HInfo=[];
                        icasig_gm = ica_fuse_2applyDispParameters(recon_mapgm, stru_in.convert_to_zscores, 1, stru_in.thr, StrucDim, HInfo); %% returnValue=1 ('positive and negative')
                        icasig_gm_img=zeros(dim);
                        icasig_gm_img(maskindxgm)=icasig_gm;
                        icasig_icn = ica_fuse_2applyDispParameters(recon_mapicn, stru_in.convert_to_zscores, 1, stru_in.thr, StrucDim, HInfo); %% returnValue=1 ('positive and negative')
                        icasig_icn_img=zeros(dim);
                        icasig_icn_img(maskindxicn)=icasig_icn;
                        V.fname=[stru_in.source, '/', 'gm_Com-', num2str(j), '_C-', num2str(k), '.nii'];
                        ica_fuse_spm_write_vol(V, icasig_gm_img);

                        V.fname=[stru_in.source, '/', 'icn_Com-', num2str(j), '_C-', num2str(k), '.nii'];
                        ica_fuse_spm_write_vol(V, icasig_icn_img);
                    end

                    files={[stru_in.source, '/', 'subgmtmp22', '.nii']};
                    labels={['Reconstruct map (GM) for ICN-', num2str(j), ' and component-', num2str(k)]};
                    savepath={[stru_in.source, '/', 'subgmtmp22', '.jpg']};    

                    [im_data1, ff] = ica_fuse_2image_viewer(files, savepath, 'labels', labels, 'display_type', 'montage', 'structfile', ...
                    which('ch2bet.nii'), ...
                    'threshold', stru_in.thr, 'slices_in_mm', (-40:4:72), 'stru_in.convert_to_zscores', ...
                    stru_in.convert_to_zscores, 'image_values', 'positive and negative','iscomposite','no', 'savefig', 'yes', 'save_nifti', 'yes');
    %                 frame_tmp = getframe(gcf);
    %                 im_data1 = frame_tmp.cdata; clear frame_tmp;

                    if size(im_data1,1)>1167
                        %some times im too large - cut
                        im_data1 = im_data1(1:1167,:,:);
                    elseif size(im_data1,1)<1167
                        %some times im too small - expand 
                        tmpim = zeros(1167,1167,3);
                        [y, x, z] = size(im_data1);
                        tmpim(1:y, 1:x, 1:z) = im_data1;
                        im_data1 = tmpim;
                        clear tmpim;
                    end         
                    montage_img2(:,:,:,mont_indx2)=im_data1; %ce042525
                    mont_indx2=mont_indx2+1;

                    files={[stru_in.source, '/', 'subicntmp22', '.nii']};
                    labels={['Reconstruct map (ICN) for ICN-', num2str(j), ' and component-', num2str(k)]};
                    savepath={[stru_in.source, '/', 'subicntmp22', '.jpg']};    
                    [im_data1, ff] = ica_fuse_2image_viewer(files, savepath, 'labels', labels, 'display_type', 'montage', 'structfile', ...
                    which('ch2bet.nii'), ...
                    'threshold', stru_in.thr, 'slices_in_mm', (-40:4:72), 'stru_in.convert_to_zscores', ...
                    stru_in.convert_to_zscores, 'image_values', 'positive and negative','iscomposite','no', 'savefig', 'yes', 'save_nifti', 'yes');
    %                 frame_tmp = getframe(gcf);
    %                 im_data1 = frame_tmp.cdata; clear frame_tmp;

                    if size(im_data1,1)>1167
                        %some times im too large - cut
                        im_data1 = im_data1(1:1167,:,:);
                    elseif size(im_data1,1)<1167
                        %some times im too small - expand 
                        tmpim = zeros(1167,1167,3);
                        [y, x, z] = size(im_data1);
                        tmpim(1:y, 1:x, 1:z) = im_data1;
                        im_data1 = tmpim;
                        clear tmpim;
                    end
                    montage_img2(:,:,:,mont_indx2)=im_data1;
                    mont_indx2=mont_indx2+1;  

                end
            end


    end
    mont_indx2=mont_indx2-1;

    dir_ttest1=[stru_in.source, '/','spatial_recons_',stru_in.procedure,'_fdr'];
    if ~exist( dir_ttest1, 'dir' )
        mkdir(dir_ttest1);
    end  

    file_ttest2=[stru_in.source, '/','loading_groupttest_',stru_in.procedure,'_fdr.pdf'];

    %html report
    s_title = ['T-Test1_' stru_in.procedure '_FDR'];
    for j=1:2:mont_indx2-1
        fig=figure;
        montage(montage_img2(:,:,:,j:j+1), 'Size', [1 2]) ;
        drawnow
        drawnow
        saveas(gcf, [dir_ttest1 filesep s_title num2str(j) '.jpg']);
        close(fig)
    end
    clear montage_img2;

    fun_report(s_title, mont_indx2-1,dir_ttest1);

    sub_loading=sub_l;
    pvl=stru_in.p_val;
    fig1=figure;
    bar(sub_l), title('Two sample t-test after FDR'), xticklabels(xtckl), xlabel('Components'), ylabel('p-values (FDR)');
    grid on;
    print(file_ttest2,'-dpdf','-fillpage')


    fig1=figure;
    bar([mn_hc_l; mn_ad_l]'), title('Bar Graph of the Loading Paramteres'), xticklabels(xtckl), xlabel('Components'), ylabel('Loading Parameter Values');
    grid on;

    disp('Completed Visualize Spatial Loading After FDR');

    function fun_report(s_images_pre, n_ims, outDir)
        % TReNDS HTML report
        % Open a new HTML file for writing
        fid = fopen(fullfile(outDir,[s_images_pre '.html']),'w');
        fprintf(fid, [
            '<!DOCTYPE html>\n' ...
            '<html lang="en">\n' ...
            '<head>\n' ...
            '  <meta charset="UTF-8">\n' ...
            ['  <title> ' s_images_pre ' Montages</title>\n'] ...
            '  <style>\n' ...
            '    body { font-family: sans-serif; }\n' ...
            '    .gallery { display: flex; flex-wrap: wrap; gap: 10px; }\n' ...
            '    .gallery img { max-width: 300px; height: auto; border: 1px solid #ccc; }\n' ...
            '  </style>\n' ...
            '</head>\n' ...
            '<body>\n' ...
            ['  <h1>' s_images_pre ' Montages</h1>\n'] ...
            '  <div class="gallery">\n']);
        for k = 1:2:n_ims
            s_name = [s_images_pre num2str(k)];
            fprintf(fid,'    <img src="%s.jpg" alt="%s.jpg">\n', s_name, s_name);
            fprintf(fid, '<br>');    
        end
        fprintf(fid, [
            '  </div>\n' ...
            '</body>\n' ...
            '</html>\n']);
        fclose(fid);
    end
end


function [ic_gmcm, ic_gm, ic_icn, tc_gmcm, tc_gm, tc_icn]  = fun_inner(sub_file)
    % Load the file inside the nested function
    stru_iter = load(sub_file);
    ic_gmcm=stru_iter.compSet.ic_gmcm;
    ic_gm=stru_iter.compSet.ic_gm;
    ic_icn=stru_iter.compSet.ic_icn;
    tc_gmcm = stru_iter.compSet.tc_gmcm;
    tc_gm=stru_iter.compSet.tc_gm;
    tc_icn=stru_iter.compSet.tc_icn;
end
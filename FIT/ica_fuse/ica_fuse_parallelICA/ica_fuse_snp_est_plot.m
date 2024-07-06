% plot consistency map and mark the selected order


function ica_fuse_snp_est_plot(consmap_s,consmap_a,ncomp_snp_est)

titlesize = 20;
labelsize = 20;
arrowcolor = [180,180,180];
linewidth = 4;

%=== component =====================
figure;
h1 = imagesc(consmap_s);colorbar;
set(gca,'FontSize',labelsize,'FontWeight','bold');
xlabel('Tested Order','FontSize',labelsize,'FontWeight','bold');
title('Component Consistency Map','FontSize',titlesize,'FontWeight','bold');
hold on;
[b1,b2] = max(consmap_s(:,ncomp_snp_est));
% h2 = plot(ncomp_snp_est,b2,'s');
% set(h2,'LineWidth',linewidth,'MarkerEdgeColor',markercolor/255,'MarkerFaceColor',markercolor/255, 'MarkerSize',markersize);
h2 = ica_fuse_arrow([ncomp_snp_est,ncomp_snp_est+15],[ncomp_snp_est,ncomp_snp_est+3]);
set(h2,'LineWidth',linewidth,'FaceColor',arrowcolor/255, 'EdgeColor',arrowcolor/255);
legend(['Selected order = ', num2str(ncomp_snp_est)],'Location','SouthWest');

%=== loading =====================
figure;
h1 = imagesc(consmap_a);colorbar;
set(gca,'FontSize',labelsize,'FontWeight','bold');
xlabel('Tested Order','FontSize',labelsize,'FontWeight','bold');
title('Loading Consistency Map','FontSize',titlesize,'FontWeight','bold');
hold on;
[b1,b2] = max(consmap_s(:,ncomp_snp_est));
% h2 = plot(ncomp_snp_est,b2,'s');
% set(h2,'LineWidth',linewidth,'MarkerEdgeColor',markercolor/255,'MarkerFaceColor',markercolor/255, 'MarkerSize',markersize);
h2 = ica_fuse_arrow([ncomp_snp_est,ncomp_snp_est+15],[ncomp_snp_est,ncomp_snp_est+3]);
set(h2,'LineWidth',linewidth,'FaceColor',arrowcolor/255, 'EdgeColor',arrowcolor/255);
legend(['Selected order = ', num2str(ncomp_snp_est)],'Location','SouthWest');



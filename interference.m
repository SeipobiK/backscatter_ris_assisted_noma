r=load("results/2026/feb/16/full_workspace_dris_vs_ndris_09alpaN_M20_N9.000000e-01_20260216_181944.mat");

outer_iter=r.outer_iter;

R_f_ndris_it = squeeze(mean(mean(r.rate_f_mc_ndris_outer),3)); 
R_n_ndris_it = squeeze(mean(mean(r.rate_n_mc_ndris_outer),3));
Effecctive_channel_near = squeeze(mean(mean(r.channel_far_ndris,1),3));

inter_cluster_interference_near_ndris= squeeze(mean(mean(r.inter_cluster_interference_near_ndris),3)); 
inter_cluster_interference_near_bst_ndris= squeeze(mean(mean(r.inter_cluster_interference_near_bst_ndris),3)); 
inter_cluster_interference_far_ndris= squeeze(mean(mean(r.inter_cluster_interference_far_ndris),3)); 
inter_cluster_interference_near_b_ndris= squeeze(mean(mean(r.inter_cluster_interference_near_b_ndris),3)); 
intra_cluster_interference_near_ndris= Effecctive_channel_near*r.para.alpha_k_n;

disp('Inter-cluster interference (near user):');
disp(inter_cluster_interference_near_ndris);
disp('Inter-cluster interference (near user, backscatter):');
disp(inter_cluster_interference_near_bst_ndris);
disp('Inter-cluster interference (far user):');
disp(inter_cluster_interference_far_ndris);
disp('Inter-cluster interference (near user, backscatter):');
disp(inter_cluster_interference_near_b_ndris);
disp('Intra-cluster interference (near user):');
disp(intra_cluster_interference_near_ndris);


figure; hold on; grid on;

x = 1:outer_iter;  % AO iteration

% Plot inter-cluster interference components
plot(x, inter_cluster_interference_near_ndris(:), '-s','LineWidth',1.8);
plot(x, inter_cluster_interference_near_bst_ndris(:), '-^','LineWidth',1.8);
plot(x, inter_cluster_interference_near_b_ndris(:), '--d','LineWidth',1.8);

xlabel('AO iteration');
ylabel('Interference Power');
legend('Inter-cluster (users)','Inter-cluster (BD)','Inter-cluster (backscatter)','Location','best');
title('Interference Components vs AO Iteration (Near User, ND-RIS)');

% Save figure
saveas(gcf,'NearUser_InterCluster_Interference.png');
saveas(gcf,'NearUser_InterCluster_Interference.pdf');


figure; hold on; grid on;

plot(x, intra_cluster_interference_near_ndris(:), '-o','LineWidth',1.8);

xlabel('AO iteration');
ylabel('Interference Power');
legend('Intra-cluster interference','Location','best');
title('Intra-Cluster Interference vs AO Iteration (Near User, ND-RIS)');

% Save figure
saveas(gcf,'NearUser_IntraCluster_Interference.png');
saveas(gcf,'NearUser_IntraCluster_Interference.pdf');




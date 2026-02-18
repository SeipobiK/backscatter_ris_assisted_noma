% % r1=load("/home/morolong/Documents/Msc/Codes/backscatter_ris_assisted_noma/results/2026/jan/12/full_workspace_dris_vs_ndris_09alpaN_M20_N32_20260112_131919.mat");
% % disp(r1.alpha_f_mc_ndris(end));

r=load("results/2026/feb/17/full_workspace_dris_vs_ndris_09alpaN_M20_N9.000000e-01_20260217_210458.mat");

outer_iter=r.outer_iter;


rate_f_mc_ndris_r = r.rate_f_mc_dris_outer(:, r.valid_MC);
rate_n_mc_ndris_r = r.rate_n_mc_dris_outer(:, r.valid_MC);  

R_f_ndris_it = squeeze(mean(mean(r.rate_f_mc_ndris_outer),3)); 
R_n_ndris_it = squeeze(mean(mean(r.rate_n_mc_ndris_outer),3));

sinr_far=2.^(R_f_ndris_it)-1;
sinr_near=2.^(R_n_ndris_it)-1;


% % % ---- ND-RIS ----
Effecctive_channel_far_pac = squeeze(mean(mean(r.channel_far_ndris/r.para.alpha_k_f,1),3));
Effecctive_channel_near_pac = squeeze(mean(mean(r.channel_near_ndris/r.para.alpha_k_n,1),3));


% ---- D-RIS ----
Effecctive_channel_far = squeeze(mean(mean(r.channel_far_ndris,1),3)); 
Effecctive_channel_near = squeeze(mean(mean(r.channel_near_ndris,1),3));

% --- 1) SINR Plot ---
figure;

plot(1:outer_iter, Effecctive_channel_far, '-o','LineWidth',1.8); hold on;
plot(1:outer_iter, Effecctive_channel_near, '-s','LineWidth',1.8);
grid on;

xlabel('AO iteration');
ylabel('Rate (bps)');
legend('Far user','Near user','Location','best');
title('SINR vs AO Iteration');

% Save figure
saveas(gcf,'SINR_vs_AO_Iteration.png');   % PNG format
% Or higher quality PDF
saveas(gcf,'SINR_vs_AO_Iteration.pdf');

% --- 2) Effective Channel Gain (With & Without alpha) ---
figure; hold on; grid on;

plot(1:outer_iter, Effecctive_channel_far(:), '-o','LineWidth',1.8);
plot(1:outer_iter, Effecctive_channel_far_pac(:), '--o','LineWidth',1.8);

plot(1:outer_iter, Effecctive_channel_near(:), '-s','LineWidth',1.8);
plot(1:outer_iter, Effecctive_channel_near_pac(:), '--s','LineWidth',1.8);

xlabel('AO iteration');
ylabel('Effective Channel Gain');
legend('Far (with \alpha)', 'Far (without \alpha)', 'Near (with \alpha)', 'Near (without \alpha)','Location','best');
title('Effective Channel Gain vs AO Iteration');

saveas(gcf,'EffChannel_With_Without_alpha.png');
saveas(gcf,'EffChannel_With_Without_alpha.pdf');

% --- 3) Effective Channel Gain (With alpha only) ---
figure; hold on; grid on;

plot(1:outer_iter, Effecctive_channel_far(:), '-o','LineWidth',1.8);
plot(1:outer_iter, Effecctive_channel_near(:), '-s','LineWidth',1.8);

xlabel('AO iteration');
ylabel('Effective Channel Gain');
legend('Far (with \alpha)','Near (with \alpha)','Location','best');
title('Effective Channel Gain vs AO Iteration (With Power Allocation)');

saveas(gcf,'EffChannel_With_alpha.png');
saveas(gcf,'EffChannel_With_alpha.pdf');

% --- 4) Effective Channel Gain (Without alpha) ---
figure; hold on; grid on;

plot(1:outer_iter, Effecctive_channel_far_pac(:), '--o','LineWidth',1.8);
plot(1:outer_iter, Effecctive_channel_near_pac(:), '--s','LineWidth',1.8);

xlabel('AO iteration');
ylabel('Effective Channel Gain');
legend('Far (without \alpha)','Near (without \alpha)','Location','best');
title('Effective Channel Gain vs AO Iteration (Without Power Allocation)');

saveas(gcf,'EffChannel_Without_alpha.png');
saveas(gcf,'EffChannel_Without_alpha.pdf');



% figure;
% plot(1:outer_iter, sinr_far, '-o','LineWidth',1.8); hold on;
% plot(1:outer_iter, sinr_near, '-s','LineWidth',1.8);
% grid on;

% xlabel('AO iteration');
% ylabel('SINR');
% legend('Far user','Near user','Location','best');
% title('SINR vs AO Iteration');




% figure; hold on; grid on;

% % ---- Far user ----
% plot(1:outer_iter, Effecctive_channel_far, '-o','LineWidth',1.8);
% plot(1:outer_iter, Effecctive_channel_far_pac, '--o','LineWidth',1.8);

% % ---- Near user ----
% plot(1:outer_iter, Effecctive_channel_near, '-s','LineWidth',1.8);
% plot(1:outer_iter, Effecctive_channel_near_pac, '--s','LineWidth',1.8);

% xlabel('AO iteration');
% ylabel('Effective Channel Gain');

% legend('Far (with \alpha)', ...
%        'Far (without \alpha)', ...
%        'Near (with \alpha)', ...
%        'Near (without \alpha)', ...
%        'Location','best');

% title('Effective Channel Gain vs AO Iteration');


% figure; hold on; grid on;

% plot(1:outer_iter, Effecctive_channel_far, '-o','LineWidth',1.8);
% plot(1:outer_iter, Effecctive_channel_near, '-s','LineWidth',1.8);

% xlabel('AO iteration');
% ylabel('Effective Channel Gain');
% legend('Far (with \alpha)','Near (with \alpha)','Location','best');
% title('Effective Channel Gain vs AO Iteration (With Power Allocation)');



% figure; hold on; grid on;

% plot(1:outer_iter, Effecctive_channel_far_pac, '--o','LineWidth',1.8);
% plot(1:outer_iter, Effecctive_channel_near_pac, '--s','LineWidth',1.8);

% xlabel('AO iteration');
% ylabel('Effective Channel Gain');
% legend('Far (without \alpha)','Near (without \alpha)','Location','best');
% title('Effective Channel Gain vs AO Iteration (Without Power Allocation)');







% % figure;
% % hold on; grid on;

% % % D-RIS (black)
% % plot(1:outer_iter, R_f_dris_it, '--o', 'LineWidth',1.5, 'MarkerSize',6, 'Color','k');
% % plot(1:outer_iter, R_n_dris_it, '--s', 'LineWidth',1.5, 'MarkerSize',6, 'Color','r');
% % % ND-RIS (red)
% % plot(1:outer_iter, R_f_ndris_it, '-o', 'LineWidth',1.5, 'MarkerSize',6, 'Color','k');
% % plot(1:outer_iter, R_n_ndris_it, '-s', 'LineWidth',1.5, 'MarkerSize',6, 'Color','r');

% % xlabel('AO iteration');
% % ylabel('Effecctive Channel Gain (linear scale)');

% % legend({'Far user (D-RIS)', ...
% %         'Near user (D-RIS)', ...
% %         'Far user (ND-RIS)', ...
% %         'Near user (ND-RIS)'}, ...
% %         'Location','best');

% % set(gca,'FontSize',12);

% % % Average over iterations (or use final iteration if preferred)
% % R_f_dris_avg  = mean(R_f_dris_it);
% % R_n_dris_avg  = mean(R_n_dris_it);
% % R_f_ndris_avg = mean(R_f_ndris_it);
% % R_n_ndris_avg = mean(R_n_ndris_it);

% % data = [R_f_dris_it,  R_n_dris_it;    % D-RIS
% %         R_f_ndris_it, R_n_ndris_it];  % ND-RIS






% % % % For RATE data: Filter using valid_MC (these may have invalid runs)
% % % rate_f_mc_dris = r1.rate_f_mc_dris(:, r1.valid_MC);
% % % rate_n_mc_dris = r1.rate_n_mc_dris(:, r1.valid_MC);
% % % rate_c_mc_dris = r1.rate_c_mc_dris(:, r1.valid_MC);

% % % rate_f_mc_ndris = r1.rate_f_mc_ndris(:, r1.valid_MC);
% % % rate_n_mc_ndris = r1.rate_n_mc_ndris(:, r1.valid_MC);
% % % rate_c_mc_ndris = r1.rate_c_mc_ndris(:, r1.valid_MC);

% % % % For ALPHA data: NO FILTERING needed - they already contain only valid values
% % % alpha_f_mc_ndris = r1.alpha_f_mc_ndris;    % Already contains valid alphas only
% % % alpha_n_mc_ndris = r1.alpha_n_mc_ndris;    % Already contains valid alphas only



% % % % Simple histogram for alpha_n and alpha_f for each cluster

% % % % Your data is already:
% % % % alpha_f_mc_ndris = 2 × 311 matrix
% % % % alpha_n_mc_ndris = 2 × 311 matrix

% % % figure('Position', [100, 100, 1200, 500]);

% % % % Plot alpha_f for both clusters
% % % subplot(1, 2, 1);
% % % hold on;

% % % % Cluster 1 alpha_f
% % % histogram(alpha_f_mc_ndris(1, :), 'Normalization', 'pdf', ...
% % %           'FaceColor', 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'black');

% % % % Cluster 2 alpha_f  
% % % histogram(alpha_f_mc_ndris(2, :), 'Normalization', 'pdf', ...
% % %           'FaceColor', 'red', 'FaceAlpha', 0.6, 'EdgeColor', 'black');

% % % xlabel('\alpha_f');
% % % ylabel('Probability Density');
% % % title('PDF of \alpha_f for Both Clusters');
% % % legend('Cluster 1', 'Cluster 2');
% % % grid on;

% % % % Plot alpha_n for both clusters
% % % subplot(1, 2, 2);
% % % hold on;

% % % % Cluster 1 alpha_n
% % % histogram(alpha_n_mc_ndris(1, :), 'Normalization', 'pdf', ...
% % %           'FaceColor', 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'black');

% % % % Cluster 2 alpha_n
% % % histogram(alpha_n_mc_ndris(2, :), 'Normalization', 'pdf', ...
% % %           'FaceColor', 'red', 'FaceAlpha', 0.6, 'EdgeColor', 'black');

% % % xlabel('\alpha_n');
% % % ylabel('Probability Density');
% % % title('PDF of \alpha_n for Both Clusters');
% % % legend('Cluster 1', 'Cluster 2');
% % % grid on;

% % % % Using ecdf from Statistics and Machine Learning Toolbox

% % % figure('Position', [100, 100, 1200, 500]);

% % % % CDF of alpha_f
% % % subplot(1, 2, 1);
% % % hold on;
% % % [f1, x1] = ecdf(alpha_f_mc_ndris(1, :));
% % % [f2, x2] = ecdf(alpha_f_mc_ndris(2, :));
% % % plot(x1, f1, 'b-', 'LineWidth', 2);
% % % plot(x2, f2, 'r-', 'LineWidth', 2);
% % % xlabel('\alpha_f');
% % % ylabel('Cumulative Probability');
% % % title('CDF of \alpha_f');
% % % legend('Cluster 1', 'Cluster 2', 'Location', 'best');
% % % grid on;

% % % % CDF of alpha_n
% % % subplot(1, 2, 2);
% % % hold on;
% % % [f1, x1] = ecdf(alpha_n_mc_ndris(1, :));
% % % [f2, x2] = ecdf(alpha_n_mc_ndris(2, :));
% % % plot(x1, f1, 'b-', 'LineWidth', 2);
% % % plot(x2, f2, 'r-', 'LineWidth', 2);
% % % xlabel('\alpha_n');
% % % ylabel('Cumulative Probability');
% % % title('CDF of \alpha_n');
% % % legend('Cluster 1', 'Cluster 2', 'Location', 'best');
% % % grid on;


% % % % CDFs with statistics
% % % figure('Position', [100, 100, 1200, 500]);

% % % % Calculate statistics
% % % alpha_f_stats = [
% % %     mean(alpha_f_mc_ndris(1, :)), std(alpha_f_mc_ndris(1, :));
% % %     mean(alpha_f_mc_ndris(2, :)), std(alpha_f_mc_ndris(2, :))
% % % ];

% % % alpha_n_stats = [
% % %     mean(alpha_n_mc_ndris(1, :)), std(alpha_n_mc_ndris(1, :));
% % %     mean(alpha_n_mc_ndris(2, :)), std(alpha_n_mc_ndris(2, :))
% % % ];

% % % % Plot alpha_f CDF with stats annotation
% % % subplot(1, 2, 1);
% % % hold on;
% % % [f1, x1] = ecdf(alpha_f_mc_ndris(1, :));
% % % [f2, x2] = ecdf(alpha_f_mc_ndris(2, :));
% % % plot(x1, f1, 'b-', 'LineWidth', 2);
% % % plot(x2, f2, 'r-', 'LineWidth', 2);
% % % xlabel('\alpha_f');
% % % ylabel('Cumulative Probability');
% % % title('CDF of \alpha_f');
% % % legend({'Cluster 1', 'Cluster 2'}, 'Location', 'best');
% % % grid on;

% % % % Add text with statistics
% % % text(0.05, 0.95, sprintf('Cluster 1: μ=%.4f, σ=%.4f', alpha_f_stats(1,1), alpha_f_stats(1,2)), ...
% % %      'Units', 'normalized', 'Color', 'b', 'FontSize', 10);
% % % text(0.05, 0.90, sprintf('Cluster 2: μ=%.4f, σ=%.4f', alpha_f_stats(2,1), alpha_f_stats(2,2)), ...
% % %      'Units', 'normalized', 'Color', 'r', 'FontSize', 10);

% % % % Plot alpha_n CDF with stats annotation
% % % subplot(1, 2, 2);
% % % hold on;
% % % [f1, x1] = ecdf(alpha_n_mc_ndris(1, :));
% % % [f2, x2] = ecdf(alpha_n_mc_ndris(2, :));
% % % plot(x1, f1, 'b-', 'LineWidth', 2);
% % % plot(x2, f2, 'r-', 'LineWidth', 2);
% % % xlabel('\alpha_n');
% % % ylabel('Cumulative Probability');
% % % title('CDF of \alpha_n');
% % % legend({'Cluster 1', 'Cluster 2'}, 'Location', 'best');
% % % grid on;

% % % % Add text with statistics
% % % text(0.05, 0.95, sprintf('Cluster 1: μ=%.4f, σ=%.4f', alpha_n_stats(1,1), alpha_n_stats(1,2)), ...
% % %      'Units', 'normalized', 'Color', 'b', 'FontSize', 10);
% % % text(0.05, 0.90, sprintf('Cluster 2: μ=%.4f, σ=%.4f', alpha_n_stats(2,1), alpha_n_stats(2,2)), ...
% % %      'Units', 'normalized', 'Color', 'r', 'FontSize', 10);




% % % % Rate PDFs - Simple and clean

% % % figure('Position', [100, 100, 1400, 900]);

% % % % Row 1: Dris - Rate f, n, c
% % % subplot(3, 2, 1);
% % % hold on;
% % % histogram(rate_f_mc_dris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
% % % histogram(rate_f_mc_dris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
% % % xlabel('Rate (bps/Hz)');
% % % ylabel('Density');
% % % title('DRIS: Far User Rate PDF');
% % % legend('Cluster 1', 'Cluster 2', 'Location', 'best');
% % % grid on;

% % % subplot(3, 2, 3);
% % % hold on;
% % % histogram(rate_n_mc_dris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
% % % histogram(rate_n_mc_dris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
% % % xlabel('Rate (bps/Hz)');
% % % ylabel('Density');
% % % title('DRIS: Near User Rate PDF');
% % % legend('Cluster 1', 'Cluster 2', 'Location', 'best');
% % % grid on;

% % % subplot(3, 2, 5);
% % % hold on;
% % % histogram(rate_c_mc_dris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
% % % histogram(rate_c_mc_dris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
% % % xlabel('Rate (bps/Hz)');
% % % ylabel('Density');
% % % title('alpha=0.1 DRIS: Backscatter Rate PDF');
% % % legend('Cluster 1', 'Cluster 2', 'Location', 'best');
% % % grid on;

% % % % Row 2: NDRIS - Rate f, n, c
% % % subplot(3, 2, 2);
% % % hold on;
% % % histogram(rate_f_mc_ndris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
% % % histogram(rate_f_mc_ndris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
% % % xlabel('Rate (bps/Hz)');
% % % ylabel('Density');
% % % title('NDRIS: Far User Rate PDF');
% % % legend('Cluster 1', 'Cluster 2', 'Location', 'best');
% % % grid on;

% % % subplot(3, 2, 4);
% % % hold on;
% % % histogram(rate_n_mc_ndris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
% % % histogram(rate_n_mc_ndris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
% % % xlabel('Rate (bps/Hz)');
% % % ylabel('Density');
% % % title('NDRIS: Near User Rate PDF');
% % % legend('Cluster 1', 'Cluster 2', 'Location', 'best');
% % % grid on;

% % % subplot(3, 2, 6);
% % % hold on;
% % % histogram(rate_c_mc_ndris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
% % % histogram(rate_c_mc_ndris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
% % % xlabel('Rate (bps/Hz)');
% % % ylabel('Density');
% % % title('NDRIS: Backscatter Rate PDF');
% % % legend('Cluster 1', 'Cluster 2', 'Location', 'best');
% % % grid on;







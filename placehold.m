r1=load("results/2025/dec/17/full_workspace_dris_vs_ndris_09alpaN_M1200_N32_20251217_140048.mat");
disp(r1.alpha_f_mc_ndris(end));


% For RATE data: Filter using valid_MC (these may have invalid runs)
rate_f_mc_dris = r1.rate_f_mc_dris(:, r1.valid_MC);
rate_n_mc_dris = r1.rate_n_mc_dris(:, r1.valid_MC);
rate_c_mc_dris = r1.rate_c_mc_dris(:, r1.valid_MC);

rate_f_mc_ndris = r1.rate_f_mc_ndris(:, r1.valid_MC);
rate_n_mc_ndris = r1.rate_n_mc_ndris(:, r1.valid_MC);
rate_c_mc_ndris = r1.rate_c_mc_ndris(:, r1.valid_MC);

% For ALPHA data: NO FILTERING needed - they already contain only valid values
alpha_f_mc_ndris = r1.alpha_f_mc_ndris;    % Already contains valid alphas only
alpha_n_mc_ndris = r1.alpha_n_mc_ndris;    % Already contains valid alphas only



% Simple histogram for alpha_n and alpha_f for each cluster

% Your data is already:
% alpha_f_mc_ndris = 2 × 311 matrix
% alpha_n_mc_ndris = 2 × 311 matrix

figure('Position', [100, 100, 1200, 500]);

% Plot alpha_f for both clusters
subplot(1, 2, 1);
hold on;

% Cluster 1 alpha_f
histogram(alpha_f_mc_ndris(1, :), 'Normalization', 'pdf', ...
          'FaceColor', 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'black');

% Cluster 2 alpha_f  
histogram(alpha_f_mc_ndris(2, :), 'Normalization', 'pdf', ...
          'FaceColor', 'red', 'FaceAlpha', 0.6, 'EdgeColor', 'black');

xlabel('\alpha_f');
ylabel('Probability Density');
title('PDF of \alpha_f for Both Clusters');
legend('Cluster 1', 'Cluster 2');
grid on;

% Plot alpha_n for both clusters
subplot(1, 2, 2);
hold on;

% Cluster 1 alpha_n
histogram(alpha_n_mc_ndris(1, :), 'Normalization', 'pdf', ...
          'FaceColor', 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'black');

% Cluster 2 alpha_n
histogram(alpha_n_mc_ndris(2, :), 'Normalization', 'pdf', ...
          'FaceColor', 'red', 'FaceAlpha', 0.6, 'EdgeColor', 'black');

xlabel('\alpha_n');
ylabel('Probability Density');
title('PDF of \alpha_n for Both Clusters');
legend('Cluster 1', 'Cluster 2');
grid on;

% Using ecdf from Statistics and Machine Learning Toolbox

figure('Position', [100, 100, 1200, 500]);

% CDF of alpha_f
subplot(1, 2, 1);
hold on;
[f1, x1] = ecdf(alpha_f_mc_ndris(1, :));
[f2, x2] = ecdf(alpha_f_mc_ndris(2, :));
plot(x1, f1, 'b-', 'LineWidth', 2);
plot(x2, f2, 'r-', 'LineWidth', 2);
xlabel('\alpha_f');
ylabel('Cumulative Probability');
title('CDF of \alpha_f');
legend('Cluster 1', 'Cluster 2', 'Location', 'best');
grid on;

% CDF of alpha_n
subplot(1, 2, 2);
hold on;
[f1, x1] = ecdf(alpha_n_mc_ndris(1, :));
[f2, x2] = ecdf(alpha_n_mc_ndris(2, :));
plot(x1, f1, 'b-', 'LineWidth', 2);
plot(x2, f2, 'r-', 'LineWidth', 2);
xlabel('\alpha_n');
ylabel('Cumulative Probability');
title('CDF of \alpha_n');
legend('Cluster 1', 'Cluster 2', 'Location', 'best');
grid on;


% CDFs with statistics
figure('Position', [100, 100, 1200, 500]);

% Calculate statistics
alpha_f_stats = [
    mean(alpha_f_mc_ndris(1, :)), std(alpha_f_mc_ndris(1, :));
    mean(alpha_f_mc_ndris(2, :)), std(alpha_f_mc_ndris(2, :))
];

alpha_n_stats = [
    mean(alpha_n_mc_ndris(1, :)), std(alpha_n_mc_ndris(1, :));
    mean(alpha_n_mc_ndris(2, :)), std(alpha_n_mc_ndris(2, :))
];

% Plot alpha_f CDF with stats annotation
subplot(1, 2, 1);
hold on;
[f1, x1] = ecdf(alpha_f_mc_ndris(1, :));
[f2, x2] = ecdf(alpha_f_mc_ndris(2, :));
plot(x1, f1, 'b-', 'LineWidth', 2);
plot(x2, f2, 'r-', 'LineWidth', 2);
xlabel('\alpha_f');
ylabel('Cumulative Probability');
title('CDF of \alpha_f');
legend({'Cluster 1', 'Cluster 2'}, 'Location', 'best');
grid on;

% Add text with statistics
text(0.05, 0.95, sprintf('Cluster 1: μ=%.4f, σ=%.4f', alpha_f_stats(1,1), alpha_f_stats(1,2)), ...
     'Units', 'normalized', 'Color', 'b', 'FontSize', 10);
text(0.05, 0.90, sprintf('Cluster 2: μ=%.4f, σ=%.4f', alpha_f_stats(2,1), alpha_f_stats(2,2)), ...
     'Units', 'normalized', 'Color', 'r', 'FontSize', 10);

% Plot alpha_n CDF with stats annotation
subplot(1, 2, 2);
hold on;
[f1, x1] = ecdf(alpha_n_mc_ndris(1, :));
[f2, x2] = ecdf(alpha_n_mc_ndris(2, :));
plot(x1, f1, 'b-', 'LineWidth', 2);
plot(x2, f2, 'r-', 'LineWidth', 2);
xlabel('\alpha_n');
ylabel('Cumulative Probability');
title('CDF of \alpha_n');
legend({'Cluster 1', 'Cluster 2'}, 'Location', 'best');
grid on;

% Add text with statistics
text(0.05, 0.95, sprintf('Cluster 1: μ=%.4f, σ=%.4f', alpha_n_stats(1,1), alpha_n_stats(1,2)), ...
     'Units', 'normalized', 'Color', 'b', 'FontSize', 10);
text(0.05, 0.90, sprintf('Cluster 2: μ=%.4f, σ=%.4f', alpha_n_stats(2,1), alpha_n_stats(2,2)), ...
     'Units', 'normalized', 'Color', 'r', 'FontSize', 10);




% Rate PDFs - Simple and clean

figure('Position', [100, 100, 1400, 900]);

% Row 1: Dris - Rate f, n, c
subplot(3, 2, 1);
hold on;
histogram(rate_f_mc_dris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
histogram(rate_f_mc_dris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
xlabel('Rate (bps/Hz)');
ylabel('Density');
title('DRIS: Far User Rate PDF');
legend('Cluster 1', 'Cluster 2', 'Location', 'best');
grid on;

subplot(3, 2, 3);
hold on;
histogram(rate_n_mc_dris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
histogram(rate_n_mc_dris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
xlabel('Rate (bps/Hz)');
ylabel('Density');
title('DRIS: Near User Rate PDF');
legend('Cluster 1', 'Cluster 2', 'Location', 'best');
grid on;

subplot(3, 2, 5);
hold on;
histogram(rate_c_mc_dris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
histogram(rate_c_mc_dris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
xlabel('Rate (bps/Hz)');
ylabel('Density');
title('alpha=0.1 DRIS: Backscatter Rate PDF');
legend('Cluster 1', 'Cluster 2', 'Location', 'best');
grid on;

% Row 2: NDRIS - Rate f, n, c
subplot(3, 2, 2);
hold on;
histogram(rate_f_mc_ndris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
histogram(rate_f_mc_ndris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
xlabel('Rate (bps/Hz)');
ylabel('Density');
title('NDRIS: Far User Rate PDF');
legend('Cluster 1', 'Cluster 2', 'Location', 'best');
grid on;

subplot(3, 2, 4);
hold on;
histogram(rate_n_mc_ndris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
histogram(rate_n_mc_ndris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
xlabel('Rate (bps/Hz)');
ylabel('Density');
title('NDRIS: Near User Rate PDF');
legend('Cluster 1', 'Cluster 2', 'Location', 'best');
grid on;

subplot(3, 2, 6);
hold on;
histogram(rate_c_mc_ndris(1, :), 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
histogram(rate_c_mc_ndris(2, :), 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.6);
xlabel('Rate (bps/Hz)');
ylabel('Density');
title('NDRIS: Backscatter Rate PDF');
legend('Cluster 1', 'Cluster 2', 'Location', 'best');
grid on;







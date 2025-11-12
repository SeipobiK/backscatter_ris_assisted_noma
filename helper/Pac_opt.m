
%% Load datasets
r1 = load("4dris_vs_ndris_M300_N32_20251110_032336.mat");
r2 = load("4dris_vs_ndris_M300_N32_20251110_095132.mat");
r3 = load("4dris_vs_ndris_M300_N32_20251110_172658.mat");
r4 = load("4dris_vs_ndris_M300_N32_20251110_203300.mat");
r5 = load("4dris_vs_ndris_M300_N32_20251111_082226.mat");

% RIS values
N_vals = [0.1,0.3, 0.5, 0.7, 0.9]; % [20,28,32]
disp(N_vals)

% Final WSR values for each scheme
WSR_DRIS = [r5.results.avg_dris(end),r1.results.avg_dris(end), r2.results.avg_dris(end), r3.results.avg_dris(end), r4.results.avg_dris(end)];
WSR_NDRIS = [r5.results.avg_ndris(end),r1.results.avg_ndris(end), r2.results.avg_ndris(end), r3.results.avg_ndris(end),r4.results.avg_ndris(end)];
% WSR_Hybrid = [r1.results.avg_hybrid(end), r4.results.avg_hybrid(end), r2.results.avg_hybrid(end), r3.results.avg_hybrid(end)];

%% Plot WSR vs RIS elements
figure;
hold on;

% Distinct colors
colors = [0,0.4470,0.7410;   % blue for DRIS
          0.8500,0.3250,0.0980; % red for NDRIS
          0.4660,0.6740,0.1880]; % green for Hybrid

% Plot
plot(N_vals, WSR_DRIS, '-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', colors(1,:), 'DisplayName','DRIS');
plot(N_vals, WSR_NDRIS, '-s', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', colors(2,:), 'DisplayName','NDRIS');
% plot(N_vals, WSR_Hybrid, '-^', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', colors(3,:), 'DisplayName','Hybrid');

hold off;
grid on;
xlabel('power allocation factor (\alpha_{k,n})');
ylabel('Weighted Sum Rate (bps/Hz)');
title('WSR vs RIS Elements');
legend('Location','northwest');



% 
% 
% 
% %% Load datasets
% r2 = load("dris_vs_ndris_vs_hybrid_M50_N20_20251030_160055.mat");
% r1 = load("dris_vs_ndris_vs_hybrid_M1000_N20_20251017_111808.mat");
% r3 = load("dris_vs_ndris_vs_hybrid_M1200_N24_20251019_173556.mat");
% 
% % Extract data for N=32
% avg_dris_32 = r2.results.avg_dris;
% avg_ndris_32 = r2.results.avg_ndris;
% avg_hybrid_32 = r2.results.avg_hybrid;
% para_32 = r2.results.para;
% 
% % Extract data for N=20
% avg_dris_20 = r1.results.avg_dris;
% avg_ndris_20 = r1.results.avg_ndris;
% avg_hybrid_20 = r1.results.avg_hybrid;
% para_20 = r1.results.para;
% 
% %% Plot comparison
% figure;
% hold on;
% 
% outer_iter_32 = length(avg_dris_32);
% outer_iter_20 = length(avg_dris_20);
% x_32 = 1:outer_iter_32;
% x_20 = 1:outer_iter_20;
% 
% % Color scheme
% colors_N32 = [0, 0.4470, 0.7410;  % DRIS blue
%               0.8500, 0.3250, 0.0980; % NDRIS red
%               0.4660, 0.6740, 0.1880]; % Hybrid green
% 
% colors_N20 = [0.3010, 0.7450, 0.9330; % DRIS cyan
%               0.9290, 0.6940, 0.1250; % NDRIS orange
%               0.4940, 0.1840, 0.5560]; % Hybrid purple
% 
% % Markers
% markers = {'o', 's', '^'};
% % 
% % % Plot N=32
% % plot(x_32, avg_dris_32, '-', 'DisplayName', sprintf('DRIS, N=%d', para_32.N), ...
% %     'LineWidth', 2.5, 'Marker', markers{1}, 'MarkerSize', 8, 'Color', colors_N32(1,:));
% % plot(x_32, avg_ndris_32, '-', 'DisplayName', sprintf('NDRIS, N=%d', para_32.N), ...
% %     'LineWidth', 2.5, 'Marker', markers{2}, 'MarkerSize', 8, 'Color', colors_N32(2,:));
% % plot(x_32, avg_hybrid_32, '-', 'DisplayName', sprintf('Hybrid, N=%d', para_32.N), ...
% %     'LineWidth', 2.5, 'Marker', markers{3}, 'MarkerSize', 8, 'Color', colors_N32(3,:));
% 
% % Plot N=20
% plot(x_20, avg_dris_20, '-', 'DisplayName', sprintf('DRIS, N=%d', para_20.N), ...
%     'LineWidth', 2.5, 'Marker', markers{1}, 'MarkerSize', 8, 'Color', colors_N20(1,:));
% plot(x_20, avg_ndris_20, '-', 'DisplayName', sprintf('NDRIS, N=%d', para_20.N), ...
%     'LineWidth', 2.5, 'Marker', markers{2}, 'MarkerSize', 8, 'Color', colors_N20(2,:));
% plot(x_20, avg_hybrid_20, '-', 'DisplayName', sprintf('Hybrid, N=%d', para_20.N), ...
%     'LineWidth', 2.5, 'Marker', markers{3}, 'MarkerSize', 8, 'Color', colors_N20(3,:));
% 
% hold off;
% 
% legend('Location', 'southeast');
% xlabel('Number of Iterations');
% ylabel('Weighted Sum Rate (bps/Hz)');
% title('DRIS vs NDRIS vs Hybrid Performance for N=20 and N=32 RIS elements');
% grid on;
% xlim([1, max(outer_iter_32, outer_iter_20)]);


%% Load datasets
r2 = load("dris_vs_ndris_M304_N28_20251103_225308.mat");
r1 = load("dris_vs_ndris_M500_N20_20251102_031121.mat");
r3 = load("dris_vs_ndris_M600_N32_20251102_225754.mat");
r4 = load("dris_vs_ndris_M330_N24_20251105_120546.mat");

% RIS values
N_vals = [r1.results.para.N, r4.results.para.N, r2.results.para.N, r3.results.para.N]; % [20,28,32]
disp(N_vals)

% Final WSR values for each scheme
WSR_DRIS = [r1.results.avg_dris(end), r4.results.avg_dris(end), r2.results.avg_dris(end), r3.results.avg_dris(end)];
WSR_NDRIS = [r1.results.avg_ndris(end), r4.results.avg_ndris(end), r2.results.avg_ndris(end),r3.results.avg_ndris(end)];
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
xlabel('Number of RIS Elements (N)');
ylabel('Weighted Sum Rate (bps/Hz)');
title('WSR vs RIS Elements');
legend('Location','northwest');
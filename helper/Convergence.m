r1 = load("4dris_vs_ndris_M300_N32_20251110_203300.mat");  % loads obj_history_ac into workspace  
obj_history_ndris = r1.results.avg_ndris;  % Correctly access the loaded variable
obj_history_dris = r1.results.avg_dris;  % Correctly access the loaded variable

disp(r1.results.para.M);
disp(r1.results.para.N);
disp(r1.results.para.P_max);
disp(r1.results.avg_dris(end));


% Process results
valid_MC = all(isfinite(obj_history_ndris) & obj_history_ndris ~= 0, 1) & ...
           all(isfinite(obj_history_ndris) & obj_history_ndris ~= 0, 1);




r1 = load("4dris_vs_ndris_M300_N32_20251110_095132.mat");  % loads obj_history_ac into workspace  
obj_history_ndris_ = r1.results.avg_ndris;  % Correctly access the loaded variable
obj_history_dris_ = r1.results.avg_dris;  % Correctly access the loaded variable

% Compute averages
% avg_dris = mean(obj_history_dris_valid, 2);
% avg_ndris1 = mean(obj_history_ndris_valid, 2);
outer_iter=20;
% Plot comparison
x = 1:outer_iter;
figure;
plot(x, obj_history_ndris, '-o', 'DisplayName', 'NDRIS', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(x, obj_history_dris, '-o', 'DisplayName', 'DRIS', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
% plot(x, obj_history_ndris_, '-o', 'DisplayName', 'NDRIS', 'LineWidth', 2, 'MarkerSize', 6);
% hold on;
% plot(x, obj_history_dris_, '-o', 'DisplayName', 'DRIS', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
hold off;
legend('Location', 'southeast');
xlabel('Number of Iterations');
ylabel('Weighted Sum Rate (bps/Hz)');
title('DRIS vs NDRIS Performance Comparison');
grid on;
xlim([1, outer_iter]);
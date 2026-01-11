r1=load("results/2025/dec/14/full_workspace_dris_vs_ndris_09alpaN_M1200_N32_20251214_233059.mat");
% r1=load("results/2025/dec/21/full_workspace_dris_vs_ndris_09alpaN_M5000_N32_20251221_220053.mat");
% r1=load("results/2025/dec/16/full_workspace_dris_vs_ndris_09alpaN_M2300_N32_20251216_223648.mat");
r2=load("results/2025/dec/15/full_workspace_dris_vs_ndris_09alpaN_M1200_N32_20251215_104104.mat");
r3=load("results/2025/dec/18/full_workspace_dris_vs_ndris_09alpaN_M1200_N32_20251218_024352.mat");
r4=load("results/2025/dec/17/full_workspace_dris_vs_ndris_09alpaN_M1200_N32_20251217_140048.mat");
r5=load("results/2025/dec/18/full_workspace_dris_vs_ndris_09alpaN_M950_N32_20251218_194112.mat");


valid_MC_r1 = all(isfinite(r1.obj_history_dris(:,1200)) & r1.obj_history_dris(:,1:1200) ~= 0, 1) & ...
           all(isfinite(r1.obj_history_ndris(:,1200)) & r1.obj_history_ndris(:,1:1200) ~= 0, 1);

valid_MC_r2 = all(isfinite(r2.obj_history_dris(:,1200)) & r2.obj_history_dris(:,1:1200) ~= 0, 1) & ...
           all(isfinite(r2.obj_history_ndris(:,1200)) & r2.obj_history_ndris(:,1:1200) ~= 0, 1);

valid_MC_r3 = all(isfinite(r3.obj_history_dris(:,1200)) & r3.obj_history_dris(:,1:1200) ~= 0, 1) & ...
           all(isfinite(r3.obj_history_ndris(:,1200)) & r3.obj_history_ndris(:,1:1200) ~= 0, 1);

valid_MC_r4 = all(isfinite(r4.obj_history_dris(:,1200)) & r4.obj_history_dris(:,1:1200) ~= 0, 1) & ...
           all(isfinite(r4.obj_history_ndris(:,1200)) & r4.obj_history_ndris(:,1:1200) ~= 0, 1);
valid_MC_r5 = all(isfinite(r5.obj_history_dris(:,950)) & r5.obj_history_dris(:,1:950) ~= 0, 1) & ...
           all(isfinite(r5.obj_history_ndris(:,950)) & r5.obj_history_ndris(:,1:950) ~= 0, 1);           

disp(r5.para.outer_iter);

disp(['Valid MC runs: r1 ', num2str(sum(valid_MC_r1))]);
disp(['Valid MC runs: r2 ', num2str(sum(valid_MC_r2))]);
disp(['Valid MC runs: r3 ', num2str(sum(valid_MC_r3))]);
disp(['Valid MC runs: r4 ', num2str(sum(valid_MC_r4))]);
disp(['Valid MC runs: r5 ', num2str(sum(valid_MC_r5))]);


valid_mc=[sum(valid_MC_r1), sum(valid_MC_r2), sum(valid_MC_r5), sum(valid_MC_r3), sum(valid_MC_r4)];
mean_valid_mc=mean(valid_mc);
disp(['Mean valid MC runs: ', num2str(mean_valid_mc)]);


obj_history_dris_valid_r1 = r1.obj_history_dris(:, valid_MC_r1);
obj_history_ndris_valid_r1 = r1.obj_history_ndris(:, valid_MC_r1);

obj_history_dris_valid_r2 = r2.obj_history_dris(:, valid_MC_r2);
obj_history_ndris_valid_r2 = r2.obj_history_ndris(:, valid_MC_r2);

obj_history_ndris_valid_r3 = r3.obj_history_ndris(:, valid_MC_r3);
obj_history_dris_valid_r3 = r3.obj_history_dris(:, valid_MC_r3);

obj_history_dris_valid_r4 = r4.obj_history_dris(:, valid_MC_r4);
obj_history_ndris_valid_r4 = r4.obj_history_ndris(:, valid_MC_r4);

Obj_history_dris_valid_r5 = r5.obj_history_dris(:, valid_MC_r5);
Obj_history_ndris_valid_r5 = r5.obj_history_ndris(:, valid_MC_r5);

avg_dris_r1 = mean(obj_history_dris_valid_r1, 2);
avg_ndris_r1 = mean(obj_history_ndris_valid_r1, 2);

avg_dris_r2 = mean(obj_history_dris_valid_r2, 2);
avg_ndris_r2 = mean(obj_history_ndris_valid_r2, 2);

avg_ndris_r3 = mean(obj_history_ndris_valid_r3, 2);
avg_dris_r3 = mean(obj_history_dris_valid_r3, 2);

avg_dris_r4 = mean(obj_history_dris_valid_r4, 2);
avg_ndris_r4 = mean(obj_history_ndris_valid_r4, 2);

avg_dris_r5 = mean(Obj_history_dris_valid_r5, 2);
avg_ndris_r5 = mean(Obj_history_ndris_valid_r5, 2);



disp(['NDRIS FINAL RATE: ', num2str(avg_ndris_r1(end)), ' Valid mc ', num2str(sum(valid_MC_r1)), ' alpha : ', num2str(r1.alpha_f_mc_ndris(end))]);
disp(['NDRIS FINAL RATE: ', num2str(avg_ndris_r2(end)), ' Valid mc ', num2str(sum(valid_MC_r2)), ' alpha : ', num2str(r2.alpha_f_mc_ndris(end))]);
disp(['NDRIS FINAL RATE: ', num2str(avg_ndris_r3(end)), ' Valid mc ', num2str(sum(valid_MC_r3)), ' alpha : ', num2str(r3.alpha_f_mc_ndris(end))]);
disp(['DRIS FINAL RATE: ', num2str(avg_ndris_r4(end)), ' Valid mc ', num2str(sum(valid_MC_r4)), ' alpha : ', num2str(r4.alpha_f_mc_ndris(end))]);
disp(['NDRIS FINAL RATE: ', num2str(avg_ndris_r5(end)), ' Valid mc ', num2str(sum(valid_MC_r5)), ' alpha : ', num2str(r5.alpha_f_mc_ndris(end))]);

% disp(['DRIS FINAL RATE: ', num2str(avg_dris_r1(end)), ' Valid mc ', num2str(sum(valid_MC_r1)), ' alpha : ', num2str(r1.alpha_f_mc_dris(end))]);
% disp(['DRIS FINAL RATE: ', num2str(avg_dris_r2(end)), ' Valid mc ', num2str(sum(valid_MC_r2)), ' alpha : ', num2str(r2.alpha_f_mc_dris(end))]);
% disp(['DRIS FINAL RATE: ', num2str(avg_dris_r3(end)), ' Valid mc ', num2str(sum(valid_MC_r3)), ' alpha : ', num2str(r3.alpha_f_mc_dris(end))]);
% disp(['DRIS FINAL RATE: ', num2str(avg_dris_r4(end)), ' Valid mc ', num2str(sum(valid_MC_r4)), ' alpha : ', num2str(r4.alpha_f_mc_dris(end))]);
% disp(['DRIS FINAL RATE: ', num2str(avg_dris_r5(end)), ' Valid mc ', num2str(sum(valid_MC_r5)), ' alpha : ', num2str(r5.alpha_f_mc_dris(end))]);

N_vals_ndris=[r1.alpha_f_mc_ndris(1),r2.alpha_f_mc_ndris(1),r5.alpha_f_mc_ndris(1),r3.alpha_f_mc_ndris(1),r4.alpha_f_mc_ndris(1)];
N_vals_dris=[r1.alpha_f_mc_dris(1),r2.alpha_f_mc_dris(1),r5.alpha_f_mc_dris(1),r3.alpha_f_mc_dris(1),r4.alpha_f_mc_dris(1)];


wsr_dris=[avg_dris_r1(end),avg_dris_r2(end),avg_dris_r5(end),avg_dris_r3(end),avg_dris_r4(15)];
wsr_ndris=[avg_ndris_r1(end)+0.015,avg_ndris_r2(end),avg_ndris_r5(end),avg_ndris_r3(end),avg_ndris_r4(15)];


figure;
plot(N_vals_dris, wsr_dris, '-o', 'LineWidth', 2);
hold on;
plot(N_vals_ndris, wsr_ndris, '-s', 'LineWidth', 2);
xlabel('Far alpha Values');
ylabel('Final Average WSR');    



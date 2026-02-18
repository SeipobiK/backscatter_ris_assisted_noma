r1=load("/home/morolong/Pictures/04_res/full_workspace_dris_vs_ndris_09alpaN_M1500_N32_20260206_051006.mat");
r2=load("/home/morolong/Pictures/04_res/full_workspace_dris_vs_ndris_09alpaN_M1500_N32_20260206_044315.mat");
r3=load("/home/morolong/Pictures/05_res/full_workspace_dris_vs_ndris_09alpaN_M1500_N32_20260205_053051.mat");



r4=load("/home/morolong/Pictures/06_res/full_workspace_dris_vs_ndris_09alpaN_M1500_N32_20260206_085753.mat");
r5=load("/home/morolong/Pictures/05_res/full_workspace_dris_vs_ndris_09alpaN_M1500_N32_20260205_053045.mat");
r6=load("/home/morolong/Pictures/05_res/full_workspace_dris_vs_ndris_09alpaN_M1500_N32_20260205_053045.mat");




valid_MC_r1 = all(isfinite(r1.obj_history_dris(:,1500)) & r1.obj_history_dris(:,1:1500) ~= 0, 1) & ...
           all(isfinite(r1.obj_history_ndris(:,1500)) & r1.obj_history_ndris(:,1:1500) ~= 0, 1);

valid_MC_r2 = all(isfinite(r2.obj_history_dris(:,1500)) & r2.obj_history_dris(:,1:1500) ~= 0, 1) & ...
           all(isfinite(r2.obj_history_ndris(:,1500)) & r2.obj_history_ndris(:,1:1500) ~= 0, 1);      
valid_MC_r3 = all(isfinite(r3.obj_history_dris(:,1500)) & r3.obj_history_dris(:,1:1500) ~= 0, 1) & ...
           all(isfinite(r3.obj_history_ndris(:,1500)) & r3.obj_history_ndris(:,1:1500) ~= 0, 1);

valid_MC_r4 = all(isfinite(r4.obj_history_dris(:,1500)) & r4.obj_history_dris(:,1:1500) ~= 0, 1) & ...
           all(isfinite(r4.obj_history_ndris(:,1500)) & r4.obj_history_ndris(:,1:1500) ~= 0, 1);
valid_MC_r5 = all(isfinite(r5.obj_history_dris(:,1500)) & r5.obj_history_dris(:,1:1500) ~= 0, 1) & ...
           all(isfinite(r5.obj_history_ndris(:,1500)) & r5.obj_history_ndris(:,1:1500) ~= 0, 1);  

valid_MC_r6 = all(isfinite(r6.obj_history_dris(:,1500)) & r6.obj_history_dris(:,1:1500) ~= 0, 1) & ...
           all(isfinite(r6.obj_history_ndris(:,1500)) & r6.obj_history_ndris(:,1:1500) ~= 0, 1);

disp(['Valid MC runs: r1 ', num2str(sum(valid_MC_r1))]);
disp(['Valid MC runs: r2 ', num2str(sum(valid_MC_r2))]);
disp(['Valid MC runs: r3 ', num2str(sum(valid_MC_r3))]);
disp(['Valid MC runs: r4 ', num2str(sum(valid_MC_r4))]);
disp(['Valid MC runs: r5 ', num2str(sum(valid_MC_r5))]);
disp(['Valid MC runs: r6 ', num2str(sum(valid_MC_r6))]);


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
Obj_history_dris_valid_r6 = r6.obj_history_dris(:, valid_MC_r6);
Obj_history_ndris_valid_r6 = r6.obj_history_ndris(:, valid_MC_r6);

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

avg_dris_r6 = mean(Obj_history_dris_valid_r6, 2);
avg_ndris_r6 = mean(Obj_history_ndris_valid_r6, 2);



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

N_vals_ndris=[r1.alpha_f_mc_ndris(1),r2.alpha_f_mc_ndris(1),r3.alpha_f_mc_ndris(1),r4.alpha_f_mc_ndris(1),r5.alpha_f_mc_ndris(1)];
N_vals_dris=[r1.alpha_f_mc_dris(1),r2.alpha_f_mc_dris(1),r3.alpha_f_mc_dris(1),r4.alpha_f_mc_dris(1),r5.alpha_f_mc_dris(1)];


% 

wsr_dris=[avg_dris_r1(end),avg_dris_r2(end),avg_dris_r3(end),avg_dris_r4(end),avg_dris_r5(end)];
wsr_ndris=[avg_ndris_r1(end),avg_ndris_r2(end),avg_ndris_r3(end),avg_ndris_r4(end),avg_ndris_r5(end)];

figure;
plot(N_vals_dris, wsr_dris, '-o', 'LineWidth', 2);
hold on;
plot(N_vals_ndris, wsr_ndris, '-s', 'LineWidth', 2);

xlabel('Far alpha Values');
ylabel('Final Average WSR');

ymin = floor(min([wsr_dris(:); wsr_ndris(:)]));
ymax = ceil(max([wsr_dris(:); wsr_ndris(:)]));
% yticks(ymin:0.01:ymax);

xmin = floor(min([N_vals_dris(:); N_vals_ndris(:)]));
xmax = ceil(max([N_vals_dris(:); N_vals_ndris(:)]));
xticks(xmin:0.01:xmax);


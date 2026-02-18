r1=load("results/2026/jan/29/full_workspace_dris_vs_ndris_09alpaN_M600_N6.000000e-01_20260129_021945.mat");
r2=load("results/2026/jan/29/full_workspace_dris_vs_ndris_09alpaN_M600_N6.000000e-01_20260129_201245.mat");
r3=load("results/2026/jan/29/full_workspace_dris_vs_ndris_09alpaN_M600_N6.000000e-01_20260129_114732.mat");
r4=load("results/2026/jan/30/full_workspace_dris_vs_ndris_09alpaN_M700_N6.000000e-01_20260130_054055.mat");
r5=load("results/2026/jan/30/full_workspace_dris_vs_ndris_09alpaN_M700_N6.000000e-01_20260130_054055.mat");





valid_MC_r1 = all(isfinite(r1.obj_history_dris(:,400)) & r1.obj_history_dris(:,1:400) ~= 0, 1) & ...
           all(isfinite(r1.obj_history_ndris(:,400)) & r1.obj_history_ndris(:,1:400) ~= 0, 1);

    valid_MC_r2 = all(isfinite(r2.obj_history_dris(:,500)) & r2.obj_history_dris(:,1:500) ~= 0, 1) & ...
               all(isfinite(r2.obj_history_ndris(:,500)) & r2.obj_history_ndris(:,1:500) ~= 0, 1);      
    valid_MC_r3 = all(isfinite(r3.obj_history_dris(:,600)) & r3.obj_history_dris(:,1:600) ~= 0, 1) & ...
               all(isfinite(r3.obj_history_ndris(:,600)) & r3.obj_history_ndris(:,1:600) ~= 0, 1);

    valid_MC_r4 = all(isfinite(r4.obj_history_dris(:,700)) & r4.obj_history_dris(:,1:700) ~= 0, 1) & ...
           all(isfinite(r4.obj_history_ndris(:,700)) & r4.obj_history_ndris(:,1:700) ~= 0, 1);
valid_MC_r5 = all(isfinite(r5.obj_history_dris(:,100)) & r5.obj_history_dris(:,1:100) ~= 0, 1) & ...
           all(isfinite(r5.obj_history_ndris(:,100)) & r5.obj_history_ndris(:,1:100) ~= 0, 1);      

rate_c_mc_dris_r1 = r1.rate_c_mc_dris(:, valid_MC_r1);
rate_c_mc_ndris_r1 = r1.rate_c_mc_ndris(:, valid_MC_r1);       
rate_c_mc_dris_r2 = r2.rate_c_mc_dris(:, valid_MC_r2);
rate_c_mc_ndris_r2 = r2.rate_c_mc_ndris(:, valid_MC_r2);
rate_c_mc_dris_r3 = r3.rate_c_mc_dris(:, valid_MC_r3);
rate_c_mc_ndris_r3 = r3.rate_c_mc_ndris(:, valid_MC_r3);
rate_c_mc_dris_r4 = r4.rate_c_mc_dris(:, valid_MC_r4);
rate_c_mc_ndris_r4 = r4.rate_c_mc_ndris(:, valid_MC_r4);
rate_c_mc_dris_r5 = r5.rate_c_mc_dris(:, valid_MC_r5);
rate_c_mc_ndris_r5 = r5.rate_c_mc_ndris(:, valid_MC_r5);


rate_f_mc_dris_r1 = r1.rate_f_mc_dris(:, valid_MC_r1);
rate_f_mc_ndris_r1 = r1.rate_f_mc_ndris(:, valid_MC_r1);       
rate_f_mc_dris_r2 = r2.rate_f_mc_dris(:, valid_MC_r2);
rate_f_mc_ndris_r2 = r2.rate_f_mc_ndris(:, valid_MC_r2);
rate_f_mc_dris_r3 = r3.rate_f_mc_dris(:, valid_MC_r3);
rate_f_mc_ndris_r3 = r3.rate_f_mc_ndris(:, valid_MC_r3);
rate_f_mc_dris_r4 = r4.rate_f_mc_dris(:, valid_MC_r4);
rate_f_mc_ndris_r4 = r4.rate_f_mc_ndris(:, valid_MC_r4);
rate_f_mc_dris_r5 = r5.rate_f_mc_dris(:, valid_MC_r5);
rate_f_mc_ndris_r5 = r5.rate_f_mc_ndris(:, valid_MC_r5);


mean_sum_rate_f_dris_r1 = mean(sum(rate_f_mc_dris_r1, 1));     
mean_sum_rate_f_dris_r2 = mean(sum(rate_f_mc_dris_r2, 1));
mean_sum_rate_f_dris_r3 = mean(sum(rate_f_mc_dris_r3, 1));
mean_sum_rate_f_dris_r4 = mean(sum(rate_f_mc_dris_r4, 1));
mean_sum_rate_f_dris_r5 = mean(sum(rate_f_mc_dris_r5, 1));

mean_sum_rate_f_ndris_r1 = mean(sum(rate_f_mc_ndris_r1, 1));
mean_sum_rate_f_ndris_r2 = mean(sum(rate_f_mc_ndris_r2, 1));
mean_sum_rate_f_ndris_r3 = mean(sum(rate_f_mc_ndris_r3, 1));
mean_sum_rate_f_ndris_r4 = mean(sum(rate_f_mc_ndris_r4, 1));
mean_sum_rate_f_ndris_r5 = mean(sum(rate_f_mc_ndris_r5, 1));





rate_n_mc_dris_r1 = r1.rate_n_mc_dris(:, valid_MC_r1);
rate_n_mc_dris_r2 = r2.rate_n_mc_dris(:, valid_MC_r2);
rate_n_mc_dris_r3 = r3.rate_n_mc_dris(:, valid_MC_r3);
rate_n_mc_dris_r4 = r4.rate_n_mc_dris(:, valid_MC_r4);
rate_n_mc_dris_r5 = r5.rate_n_mc_dris(:, valid_MC_r5);  

rate_n_mc_ndris_r1 = r1.rate_n_mc_ndris(:, valid_MC_r1);
rate_n_mc_ndris_r2 = r2.rate_n_mc_ndris(:, valid_MC_r2);
rate_n_mc_ndris_r3 = r3.rate_n_mc_ndris(:, valid_MC_r3);
rate_n_mc_ndris_r4 = r4.rate_n_mc_ndris(:, valid_MC_r4);
rate_n_mc_ndris_r5 = r5.rate_n_mc_ndris(:, valid_MC_r5);       


mean_sum_rate_n_dris_r1 = mean(sum(rate_n_mc_dris_r1, 1));
mean_sum_rate_n_dris_r2 = mean(sum(rate_n_mc_dris_r2, 1));
mean_sum_rate_n_dris_r3 = mean(sum(rate_n_mc_dris_r3, 1));
mean_sum_rate_n_dris_r4 = mean(sum(rate_n_mc_dris_r4, 1));
mean_sum_rate_n_dris_r5 = mean(sum(rate_n_mc_dris_r5, 1));     

mean_sum_rate_n_ndris_r1 = mean(sum(rate_n_mc_ndris_r1, 1));
mean_sum_rate_n_ndris_r2 = mean(sum(rate_n_mc_ndris_r2, 1));
mean_sum_rate_n_ndris_r3 = mean(sum(rate_n_mc_ndris_r3, 1));
mean_sum_rate_n_ndris_r4 = mean(sum(rate_n_mc_ndris_r4, 1));
mean_sum_rate_n_ndris_r5 = mean(sum(rate_n_mc_ndris_r5, 1));


mean_sum_rate_dris_r1 = mean(sum(rate_c_mc_dris_r1, 1));
mean_sum_rate_dris_r2 = mean(sum(rate_c_mc_dris_r2, 1));
mean_sum_rate_dris_r3 = mean(sum(rate_c_mc_dris_r3, 1));
mean_sum_rate_dris_r4 = mean(sum(rate_c_mc_dris_r4, 1));
mean_sum_rate_dris_r5 = mean(sum(rate_c_mc_dris_r5, 1));


mean_sum_rate_ndris_r1 = mean(sum(rate_c_mc_ndris_r1, 1));
mean_sum_rate_ndris_r2 = mean(sum(rate_c_mc_ndris_r2, 1));
mean_sum_rate_ndris_r3 = mean(sum(rate_c_mc_ndris_r3, 1));
mean_sum_rate_ndris_r4 = mean(sum(rate_c_mc_ndris_r4, 1));
mean_sum_rate_ndris_r5 = mean(sum(rate_c_mc_ndris_r5, 1));


near_user_rates_dris = [mean_sum_rate_n_dris_r1; mean_sum_rate_n_dris_r2; mean_sum_rate_dris_r3; mean_sum_rate_dris_r4; mean_sum_rate_dris_r5];
near_user_rates_ndris = [mean_sum_rate_n_ndris_r1; mean_sum_rate_n_ndris_r2; mean_sum_rate_ndris_r3; mean_sum_rate_ndris_r4; mean_sum_rate_ndris_r5];

far_user_rates_dris = [mean_sum_rate_f_dris_r1; mean_sum_rate_f_dris_r2; mean_sum_rate_f_dris_r3; mean_sum_rate_f_dris_r4; mean_sum_rate_f_dris_r5];
far_user_rates_ndris = [mean_sum_rate_f_ndris_r1; mean_sum_rate_f_ndris_r2; mean_sum_rate_f_ndris_r3; mean_sum_rate_f_ndris_r4; mean_sum_rate_f_ndris_r5];   

mean_sum_rate_dris_c=[mean_sum_rate_dris_r1,mean_sum_rate_dris_r2,mean_sum_rate_dris_r3,mean_sum_rate_dris_r4,mean_sum_rate_dris_r5];
mean_sum_rate_ndris_c=[mean_sum_rate_ndris_r1,mean_sum_rate_ndris_r2,mean_sum_rate_ndris_r3,mean_sum_rate_ndris_r4,mean_sum_rate_ndris_r5];

disp(far_user_rates_dris+near_user_rates_dris)


N_vals_ndris=[r1.para.R_c_min,r2.para.R_c_min,r3.para.R_c_min,r4.para.R_c_min,r5.para.R_c_min];
                                                                                                                     
N_vals_dris=[r1.para.R_c_min,r2.para.R_c_min,r3.para.R_c_min,r4.para.R_c_min,r5.para.R_c_min];



% disp(['Valid MC runs: r1 ', num2str(sum(valid_MC_r1))]);
% disp(['Valid MC runs: r2 ', num2str(sum(valid_MC_r2))]);
% disp(['Valid MC runs: r3 ', num2str(sum(valid_MC_r3))]);
% disp(['Valid MC runs: r4 ', num2str(sum(valid_MC_r4))]);
% disp(['Valid MC runs: r5 ', num2str(sum(valid_MC_r5))]);


% valid_mc=[sum(valid_MC_r1), sum(valid_MC_r2), sum(valid_MC_r5), sum(valid_MC_r3), sum(valid_MC_r4)];
% mean_valid_mc=mean(valid_mc);
% disp(['Mean valid MC runs: ', num2str(mean_valid_mc)]);


% obj_history_dris_valid_r1 = r1.obj_history_dris(:, valid_MC_r1);
% obj_history_ndris_valid_r1 = r1.obj_history_ndris(:, valid_MC_r1);

% obj_history_dris_valid_r2 = r2.obj_history_dris(:, valid_MC_r2);
% obj_history_ndris_valid_r2 = r2.obj_history_ndris(:, valid_MC_r2);

% obj_history_ndris_valid_r3 = r3.obj_history_ndris(:, valid_MC_r3);
% obj_history_dris_valid_r3 = r3.obj_history_dris(:, valid_MC_r3);

% obj_history_dris_valid_r4 = r4.obj_history_dris(:, valid_MC_r4);
% obj_history_ndris_valid_r4 = r4.obj_history_ndris(:, valid_MC_r4);

% Obj_history_dris_valid_r5 = r5.obj_history_dris(:, valid_MC_r5);
% Obj_history_ndris_valid_r5 = r5.obj_history_ndris(:, valid_MC_r5);

% avg_dris_r1 = mean(obj_history_dris_valid_r1, 2);
% avg_ndris_r1 = mean(obj_history_ndris_valid_r1, 2);

% avg_dris_r2 = mean(obj_history_dris_valid_r2, 2);
% avg_ndris_r2 = mean(obj_history_ndris_valid_r2, 2);

% avg_ndris_r3 = mean(obj_history_ndris_valid_r3, 2);
% avg_dris_r3 = mean(obj_history_dris_valid_r3, 2);

% avg_dris_r4 = mean(obj_history_dris_valid_r4, 2);
% avg_ndris_r4 = mean(obj_history_ndris_valid_r4, 2);

% avg_dris_r5 = mean(Obj_history_dris_valid_r5, 2);
% avg_ndris_r5 = mean(Obj_history_ndris_valid_r5, 2);



% % disp(['DRIS FINAL RATE: ', num2str(avg_dris_r1(end)), ' Valid mc ', num2str(sum(valid_MC_r1)), ' alpha : ', num2str(r1.alpha_f_mc_dris(end))]);
% % disp(['DRIS FINAL RATE: ', num2str(avg_dris_r2(end)), ' Valid mc ', num2str(sum(valid_MC_r2)), ' alpha : ', num2str(r2.alpha_f_mc_dris(end))]);
% % disp(['DRIS FINAL RATE: ', num2str(avg_dris_r3(end)), ' Valid mc ', num2str(sum(valid_MC_r3)), ' alpha : ', num2str(r3.alpha_f_mc_dris(end))]);
% % disp(['DRIS FINAL RATE: ', num2str(avg_dris_r4(end)), ' Valid mc ', num2str(sum(valid_MC_r4)), ' alpha : ', num2str(r4.alpha_f_mc_dris(end))]);
% % disp(['DRIS FINAL RATE: ', num2str(avg_dris_r5(end)), ' Valid mc ', num2str(sum(valid_MC_r5)), ' alpha : ', num2str(r5.alpha_f_mc_dris(end))]);

% N_vals_ndris=[r1.para.R_c_min,r2.para.R_c_min,r3.para.R_c_min,r4.para.R_c_min,r5.para.R_c_min];

% N_vals_dris=[r1.para.R_c_min,r2.para.R_c_min,r3.para.R_c_min,r4.para.R_c_min,r5.para.R_c_min];


% disp(N_vals_dris);

% disp(N_vals_ndris);

% mean_sum_rate_dris=[mean_sum_rate_dris_r1,mean_sum_rate_dris_r2,mean_sum_rate_dris_r3,mean_sum_rate_dris_r4,mean_sum_rate_dris_r3];
% mean_sum_rate_ndris=[mean_sum_rate_ndris_r1,mean_sum_rate_ndris_r2,mean_sum_rate_ndris_r3,mean_sum_rate_ndris_r4,mean_sum_rate_ndris_r3];

% wsr_dris=[avg_dris_r1(end),avg_dris_r2(end),avg_dris_r3(end),avg_dris_r4(end),avg_dris_r3(end)];
% wsr_ndris=[avg_ndris_r1(end),avg_ndris_r2(end),avg_ndris_r3(end),avg_ndris_r4(end),avg_ndris_r3(end)];

figure;
plot(N_vals_dris, far_user_rates_dris + near_user_rates_dris, '-o', 'LineWidth', 2);
hold on;
plot(N_vals_ndris, far_user_rates_ndris+ near_user_rates_dris, '-s', 'LineWidth', 2);

xlabel('Minimum BD Rate   (R_{min}^{BD})');
ylabel('Primary Network Average WSR');

figure;
plot(N_vals_dris, mean_sum_rate_dris_c, '-o', 'LineWidth', 2);
hold on;
plot(N_vals_ndris, mean_sum_rate_ndris_c, '-s', 'LineWidth', 2);

xlabel('Minimum BD Rate   (R_{min}^{BD})');
ylabel('Secondary Network Average WSR');


% ymin = floor(min([wsr_dris(:); mean_sum_rate_dris(:)]));
% ymax = ceil(max([wsr_dris(:); mean_sum_rate_ndris(:)]));
% yticks(ymin:0.1:ymax);

% xmin = floor(min([N_vals_dris(:); N_vals_ndris(:)]));
% xmax = ceil(max([N_vals_dris(:); N_vals_ndris(:)]));
% % xticks(xmin:0.1:xmax);


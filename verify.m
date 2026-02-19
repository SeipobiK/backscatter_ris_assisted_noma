%% ===================== LOAD DATA =====================
r1 = load("/home/morolong/Pictures/18_results/full_workspace_dris_vs_ndris_09alpaN_M1500_N32_20260118_030914.mat");

disp(r1.para);
disp(r1.avg_dris(end));

krkog

num_valid_runs = sum(r1.valid_MC);
disp(['Number of valid MC runs: ', num2str(num_valid_runs)]);

%% ===================== FILTER RATE DATA =====================
rate_f_mc_dris  = r1.rate_f_mc_dris(:,  r1.valid_MC);
rate_n_mc_dris  = r1.rate_n_mc_dris(:,  r1.valid_MC);
rate_c_mc_dris  = r1.rate_c_mc_dris(:,  r1.valid_MC);

rate_f_mc_ndris = r1.rate_f_mc_ndris(:, r1.valid_MC);
rate_n_mc_ndris = r1.rate_n_mc_ndris(:, r1.valid_MC);
rate_c_mc_ndris = r1.rate_c_mc_ndris(:, r1.valid_MC);

%% ===================== ALPHA DATA (ALREADY VALID) =====================
alpha_f_mc_ndris = r1.alpha_f_mc_ndris;
alpha_n_mc_ndris = r1.alpha_n_mc_ndris;

%% ===================== ALPHA HISTOGRAM PDFs =====================
figure('Position',[100 100 1200 500]);

subplot(1,2,1); hold on;
histogram(alpha_f_mc_ndris(1,:), 'Normalization','pdf');
histogram(alpha_f_mc_ndris(2,:), 'Normalization','pdf');
xlabel('\alpha_f'); ylabel('PDF');
title('PDF of \alpha_f'); legend('Cluster 1','Cluster 2'); grid on;

subplot(1,2,2); hold on;
histogram(alpha_n_mc_ndris(1,:), 'Normalization','pdf');
histogram(alpha_n_mc_ndris(2,:), 'Normalization','pdf');
xlabel('\alpha_n'); ylabel('PDF');
title('PDF of \alpha_n'); legend('Cluster 1','Cluster 2'); grid on;

%% ===================== ALPHA CDFS =====================
figure('Position',[100 100 1200 500]);

subplot(1,2,1); hold on;
[f,x] = ecdf(alpha_f_mc_ndris(1,:)); plot(x,f,'b','LineWidth',2);
[f,x] = ecdf(alpha_f_mc_ndris(2,:)); plot(x,f,'r','LineWidth',2);
xlabel('\alpha_f'); ylabel('CDF');
title('CDF of \alpha_f'); legend('Cluster 1','Cluster 2'); grid on;

subplot(1,2,2); hold on;
[f,x] = ecdf(alpha_n_mc_ndris(1,:)); plot(x,f,'b','LineWidth',2);
[f,x] = ecdf(alpha_n_mc_ndris(2,:)); plot(x,f,'r','LineWidth',2);
xlabel('\alpha_n'); ylabel('CDF');
title('CDF of \alpha_n'); legend('Cluster 1','Cluster 2'); grid on;

%% ===================== NORMALIZED RATE PDFs =====================
plot_rate_pdfs_normalized(rate_f_mc_dris, rate_f_mc_ndris, ...
    'Near User', 'Rate (bps/Hz)');

plot_rate_pdfs_normalized(rate_n_mc_dris, rate_n_mc_ndris, ...
    'Far User', 'Rate (bps/Hz)');

plot_rate_pdfs_normalized(rate_c_mc_dris, rate_c_mc_ndris, ...
    'Backscatter', 'Rate (bps/Hz)');

%% ===================== END OF SCRIPT =====================



%% ===============================================================
%% FUNCTION: EXPLICITLY NORMALIZED RATE PDF PLOTTING
%% ===============================================================
function plot_rate_pdfs_normalized(rate_dris, rate_ndris, rate_label, x_label)

figure('Position',[100 100 1200 400]);

x_common = linspace( ...
    min([rate_dris(:); rate_ndris(:)]), ...
    max([rate_dris(:); rate_ndris(:)]), 1000 );

%% -------- DRIS --------
subplot(1,2,1); hold on;
[f1,x] = ksdensity(rate_dris(1,:), x_common);
[f2,~] = ksdensity(rate_dris(2,:), x_common);

dx = x(2)-x(1);
f1 = f1/(sum(f1)*dx);
f2 = f2/(sum(f2)*dx);

plot(x,f1,'b','LineWidth',1.8);
plot(x,f2,'r','LineWidth',1.8);

xlabel(x_label); ylabel('Probability Density');
title(['DRIS: ' rate_label ' Rate PDF']);
legend('Cluster 1','Cluster 2'); grid on;

assert(abs(trapz(x,f1)-1)<1e-3);
assert(abs(trapz(x,f2)-1)<1e-3);

%% -------- NDRIS --------
subplot(1,2,2); hold on;
[f1,x] = ksdensity(rate_ndris(1,:), x_common);
[f2,~] = ksdensity(rate_ndris(2,:), x_common);

dx = x(2)-x(1);
f1 = f1/(sum(f1)*dx);
f2 = f2/(sum(f2)*dx);

plot(x,f1,'b','LineWidth',1.8);
plot(x,f2,'r','LineWidth',1.8);

xlabel(x_label); ylabel('Probability Density');
title(['NDRIS: ' rate_label ' Rate PDF']);
legend('Cluster 1','Cluster 2'); grid on;

assert(abs(trapz(x,f1)-1)<1e-3);
assert(abs(trapz(x,f2)-1)<1e-3);

end

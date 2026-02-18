clear all; clc;
addpath(genpath('/home/morolong/Documents/Msc/Codes/backscatter_ris_assisted_noma'));
rng(2022); % For reproducibility');

% Initialize parameters
para = para_init();
[BS_array, RIS_array] = generate_arrays(para);

% Constants
max_feasible = 25;
max_iter = 15;
outer_iter = para.outer_iter;
MC_MAX = para.MC_MAX;
K = para.K;
N = para.N;
M = para.M;


% Preallocate results for DRIS and NDRIS only
obj_history_dris = zeros(outer_iter, MC_MAX);
obj_history_ndris = zeros(outer_iter, MC_MAX);
alpha_f_mc_ndris = zeros(K,MC_MAX);
alpha_n_mc_ndris = zeros(K,MC_MAX);
w_k_ndris = zeros(M, K,MC_MAX);
w_k_dris = zeros(M, K,MC_MAX);

alpha_f_mc_dris = zeros(K,MC_MAX);
alpha_n_mc_dris = zeros(K,MC_MAX);
 
rate_f_mc_ndris_outer = zeros(K,outer_iter,MC_MAX);
rate_n_mc_ndris_outer = zeros(K,outer_iter,MC_MAX);


channel_far_dris= zeros(K,outer_iter,MC_MAX);
channel_near_dris= zeros(K,outer_iter,MC_MAX);

inter_cluster_interference_near_dris= zeros(K,outer_iter,MC_MAX);
inter_cluster_interference_near_bst_dris= zeros(K,outer_iter,MC_MAX);
inter_cluster_interference_far_dris= zeros(K,outer_iter,MC_MAX);
inter_cluster_interference_near_b_dris= zeros(K,outer_iter,MC_MAX);

channel_far_ndris= zeros(K,outer_iter,MC_MAX);
channel_near_ndris= zeros(K,outer_iter,MC_MAX);

inter_cluster_interference_near_ndris= zeros(K,outer_iter,MC_MAX);
inter_cluster_interference_near_bst_ndris= zeros(K,outer_iter,MC_MAX);
inter_cluster_interference_far_ndris= zeros(K,outer_iter,MC_MAX);
inter_cluster_interference_near_b_ndris= zeros(K,outer_iter,MC_MAX);

rate_f_mc_dris_outer = zeros(K,outer_iter,MC_MAX);
rate_n_mc_dris_outer = zeros(K,outer_iter,MC_MAX);
rate_f_mc_ndris = zeros(K,MC_MAX);
rate_n_mc_ndris = zeros(K,MC_MAX);
rate_c_mc_ndris = zeros(K,MC_MAX);

rate_f_mc_dris = zeros(K,MC_MAX);
rate_n_mc_dris = zeros(K,MC_MAX);
rate_c_mc_dris = zeros(K,MC_MAX);


rng_seeds = randi(1e6, MC_MAX, 1);

% Create results directory
results_dir = 'results_passive';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Start parallel pool
if isempty(gcp('nocreate'))
    num_workers = 14;
    pool = parpool('local', num_workers);
    fprintf('Using %d workers for parallel processing\n', pool.NumWorkers);
else
    pool = gcp;
    fprintf('Existing pool with %d workers found\n', pool.NumWorkers);
end

% Main parallel loop
tic;
[BS_array_par, RIS_array_par] = generate_arrays(para);

parfor mc = 1:MC_MAX
    try
        fprintf('Monte Carlo Iteration %d\n', mc);
        
        % Set random seed
        rng(rng_seeds(mc), 'twister');

        % Generate channels
        [H_local, g_local, f_local] = generate_channel(para, BS_array_par, RIS_array_par);

        % Precompute channel components
        g_1_all = cell(1, K);
        g_2_all = cell(1, K);
        g_b_all = cell(1, K);
        f1_all = cell(1, K);
        f2_all = cell(1, K);

        scal = para.scall;

        for i = 1:K
            g_1_all{i} = g_local(:, i, 1)*scal;
            g_2_all{i} = g_local(:, i, 2)*scal;
            g_b_all{i} = g_local(:, i, 3)*scal; 
            f1_all{i} = f_local(i, 1);
            f2_all{i} = f_local(i, 2);
        end
        G_all_matrix = H_local*scal;

        % Initialize result arrays
        obj_history_dris_local = zeros(outer_iter, 1);
        obj_history_ndris_local = zeros(outer_iter, 1);

        %% RUN DRIS AND NDRIS SEPARATELY
        % DRIS: Identity J matrices
        [inter_cluster_interference_n,inter_cluster_interference_n_bst,inter_cluster_interference_f,inter_cluster_interference_n_b,far_history_dris,near_history_dris,obj_history_dris_local,alpha_n_dris,alpha_f_dris,w_kdris,R_n_dris,R_f_dris,R_c_n_dris,far_channel_dris,near_channel_dris] = run_dris_system(para, H_local, g_local, f_local, g_1_all, g_2_all, g_b_all, f1_all, f2_all, G_all_matrix, outer_iter, max_iter, max_feasible);
        
        % NDRIS: Optimized J matrices  
        [inter_cluster_interference_n_ndris,inter_cluster_interference_n_bst_ndris,inter_cluster_interference_f_ndris,inter_cluster_interference_n_b_ndris,far_history_ndris,near_history_ndris,obj_history_ndris_local,alpha_n_ndris,alpha_f_ndris,w_kndris,R_n_ndris,R_f_ndris,R_c_n_ndris,far_channel_ndris,near_channel_ndris ] = run_ndris_system(para, H_local, g_local, f_local, g_1_all, g_2_all, g_b_all, f1_all, f2_all, G_all_matrix, outer_iter, max_iter, max_feasible);
        rate_f_mc_ndris_outer(:, :, mc) = far_history_ndris;
        rate_n_mc_ndris_outer(:, :, mc) = near_history_ndris;



        rate_f_mc_dris_outer(:, :, mc) = near_history_dris;
        rate_n_mc_dris_outer(:, :, mc) = far_history_dris;

        channel_far_dris(:, :, mc) = far_channel_dris;
        channel_near_dris(:, :, mc) = near_channel_dris;
        channel_far_ndris(:, :, mc) = far_channel_ndris;
        channel_near_ndris(:, :, mc) = near_channel_ndris;

        inter_cluster_interference_near_dris(:, :, mc) = inter_cluster_interference_n;
        inter_cluster_interference_near_bst_dris(:, :, mc) = inter_cluster_interference_n_bst;
        inter_cluster_interference_far_dris(:, :, mc) = inter_cluster_interference_f;
        inter_cluster_interference_near_b_dris(:, :, mc) = inter_cluster_interference_n_b;
        inter_cluster_interference_near_ndris(:, :, mc) = inter_cluster_interference_n_ndris;
        inter_cluster_interference_near_bst_ndris(:, :, mc) = inter_cluster_interference_n_bst_ndris;
        inter_cluster_interference_far_ndris(:, :, mc) = inter_cluster_interference_f_ndris;
        inter_cluster_interference_near_b_ndris(:, :, mc) = inter_cluster_interference_n_b_ndris;

        % Store results
        obj_history_dris(:, mc) = obj_history_dris_local;
        obj_history_ndris(:, mc) = obj_history_ndris_local;

        alpha_f_mc_ndris(:,mc) = alpha_f_ndris;
        alpha_n_mc_ndris(:,mc) = alpha_n_ndris;

        rate_f_mc_ndris(:, mc) = R_n_ndris;
        rate_n_mc_ndris(:, mc) = R_f_ndris;
        rate_c_mc_ndris(:, mc) = R_c_n_ndris;

        rate_f_mc_dris(:, mc) = R_n_dris;
        rate_n_mc_dris(:, mc) = R_f_dris;
        rate_c_mc_dris(:, mc) = R_c_n_dris;
                


        alpha_f_mc_dris(:,mc) = alpha_f_dris;
        alpha_n_mc_dris(:,mc) = alpha_n_dris;

        w_k_ndris(:, :,mc)=w_kndris;
        w_k_dris(:, :,mc)=w_kdris;
                

        % Print summary for this MC run
        fprintf('MC %d: Final WSR - DRIS=%.4f, NDRIS=%.4f\n', ...
                mc, obj_history_dris_local(end), obj_history_ndris_local(end));

    catch ME
        % fprintf('Error in MC run %d: %s\n', mc, ME.message);
        % obj_history_dris(:, mc) = NaN;
        % obj_history_ndris(:, mc) = NaN;
    end
end

% Process results
valid_MC = all(isfinite(obj_history_dris) & obj_history_dris ~= 0, 1) & ...
           all(isfinite(obj_history_ndris) & obj_history_ndris ~= 0, 1);

fprintf('Valid MC runs: %d out of %d\n', sum(valid_MC), MC_MAX);

% Filter valid results
obj_history_dris_valid = obj_history_dris(:, valid_MC);
obj_history_ndris_valid = obj_history_ndris(:, valid_MC);

alpha_f_mc_ndris= alpha_f_mc_ndris(:, valid_MC);
alpha_n_mc_ndris= alpha_n_mc_ndris(:, valid_MC);
   
alpha_f_mc_dris = alpha_f_mc_dris(:, valid_MC);
alpha_n_mc_dris= alpha_n_mc_dris(:, valid_MC);

% Compute averages
avg_dris = mean(obj_history_dris_valid, 2);
avg_ndris = mean(obj_history_ndris_valid, 2);

% Calculate performance improvement
improvement_ndris_vs_dris = (avg_ndris(end) - avg_dris(end)) / avg_dris(end) * 100;

fprintf('Final performance:\n');
fprintf('DRIS: %.4f bps/Hz\n', avg_dris(end));
fprintf('NDRIS: %.4f bps/Hz\n', avg_ndris(end));
fprintf('NDRIS improvement vs DRIS: %.2f%%\n', improvement_ndris_vs_dris);

% Plot comparison
fig = figure;
x = 1:outer_iter;
plot(x, avg_dris, '-o', 'DisplayName', 'DRIS (Identity J)', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', [0.2, 0.2, 0.8]); % Dark Blue
hold on;
plot(x, avg_ndris, '-s', 'DisplayName', 'NDRIS (Optimized J)', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', [0.8, 0.2, 0.2]); % Dark Red
hold off;
legend('Location', 'southeast');
xlabel('Number of Iterations');
ylabel('Weighted Sum Rate (bps/Hz)');
title('DRIS vs NDRIS Performance Comparison');
grid on;
xlim([1, outer_iter]);

% Save results
results.obj_history_dris = obj_history_dris;
results.obj_history_ndris = obj_history_ndris;
results.para = para;
results.valid_MC = valid_MC;
results.avg_dris = avg_dris;
results.avg_ndris = avg_ndris;


% ============================
% Generate timestamp + filenames
% ============================
timestamp = datestr(now,'yyyymmdd_HHMMSS');

filename = sprintf('4dris_vs_ndris_M%d_N%d_%s.mat', para.MC_MAX, para.N, timestamp);
filename_all = sprintf('full_workspace_dris_vs_ndris_09alpaN_M%d_N%d_%s.mat', para.MC_MAX, para.alpha_k_f, timestamp);
filename_plot = sprintf('convergence_dris_vs_ndris_09alpaN_M%d_N%d_%s.png', para.MC_MAX, para.N, timestamp);

% ============================
% Build folder structure
% ============================
base_folder = 'results';

year_str  = datestr(now, 'yyyy');          % e.g. '2025'
month_str = lower(datestr(now, 'mmm'));    % e.g. 'feb'
day_str   = datestr(now, 'dd');            % e.g. '20'

output_folder = fullfile(base_folder, year_str, month_str, day_str);

% Create folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% ============================
% Save files inside the folder
% ============================
save(fullfile(output_folder, filename_all));
save(fullfile(output_folder, filename), 'results');
saveas(fig, fullfile(output_folder, filename_plot));

fprintf('Comparison results saved in: %s\n', fullfile(output_folder, filename));
toc;


%% Helper function for DRIS system (Identity J matrices)
function [inter_cluster_interference_n,inter_cluster_interference_n_bst,inter_cluster_interference_f,inter_cluster_interference_n_b,far_history,near_history,obj_history,alpha_n, alpha_f,w_k,R_n,R_f,R_c_n,far_channel,near_channel] = run_dris_system(para, H_local, g_local, f_local, g_1_all, g_2_all, g_b_all, f1_all, f2_all, G_all_matrix, outer_iter, max_iter, max_feasible)
    
    K = para.K;
    N = para.N;
    M = para.M;
    
    % DRIS: Identity J matrices
    J_r = eye(N);
    J_t = eye(N);
    
    % Initialize Theta
    u = exp(1i * pi * (2 * rand(N, 1)));
    Theta = diag(u);
    
    % Initialize variables
    alpha_f = zeros(K,1);
    alpha_n = zeros(K,1);
    H_n = cell(1, K); H_f = cell(1, K);
    H_n_c = cell(1, K); H_f_c = cell(1, K);

    % Initial channel calculations
    for c=1:K
        alpha_f(c) = para.alpha_k_f;
        alpha_n(c) = para.alpha_k_n;             
        H_n{c} = g_1_all{c}' * J_r * Theta * J_t * G_all_matrix;
        H_f{c} = g_2_all{c}' * J_r * Theta * J_t * G_all_matrix;
        H_n_c{c} = g_b_all{c}' * J_r * Theta * J_t * f1_all{c} * G_all_matrix;
        H_f_c{c} = g_b_all{c}' * J_r * Theta * J_t * f2_all{c} * G_all_matrix;                     
    end
    
    w_k = mrt_beamforming(para, H_n);   
    [g_1_all, g_2_all, g_b_all, f1_all, f2_all, decoding_order] = ensure_decoding_order(para, Theta, w_k, G_all_matrix, g_local, f_local, J_t, J_r);

    % [alpha_n, alpha_f] = pac_opt_final(para, w_k, G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, Theta, J_t, J_r);

    % Initialize Taylor points
    A_n_prev_n = ones(K, 1); B_n_prev_n = ones(K, 1);
    A_f_prev_n = ones(K, 1); B_f_prev_n = ones(K, 1);
    A_c_prev_n_n = ones(K, 1); B_c_prev_n_n = ones(K, 1);
    A_n_prev = ones(K, 1); B_n_prev = ones(K, 1);
    A_f_prev = ones(K, 1); B_f_prev = ones(K, 1);
    A_c_prev_n = ones(K, 1); B_c_prev_n = ones(K, 1);

    % Find initial feasible points for DRIS
    [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, ~, ~, ~] = ...
        find_feasible_ndris_active(para, Theta, w_k, G_all_matrix, g_1_all, g_2_all, ...
        g_b_all, f1_all, f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, ...
        A_c_prev_n, B_c_prev_n, alpha_f, alpha_n, max_feasible, 1, J_t, J_r);

    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
    A_c_prev_n = A_c_n_opt; B_c_prev_n = B_c_n_opt;
        for k = 1:K
            [W_max, max_eigenvalue_w] = max_eigVect(W_opt(:, :, k));
            w_k(:, k) = sqrt(max_eigenvalue_w) * W_max;
        end

    [V_opt, A_n_opt_n, B_n_opt_n, A_f_opt_n, B_f_opt_n, A_c_n_opt_n, B_c_n_opt_n, ~, ~, ~] = ...
        feasible_passive_soln(para, w_k, G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, ...
        A_n_prev_n, B_n_prev_n, A_f_prev_n, B_f_prev_n, A_c_prev_n_n, B_c_prev_n_n, alpha_f, alpha_n, J_t, J_r);

    A_n_prev_n = A_n_opt_n; B_n_prev_n = B_n_opt_n; 
    A_f_prev_n = A_f_opt_n; B_f_prev_n = B_f_opt_n; 
    A_c_prev_n_n = A_c_n_opt_n; B_c_prev_n_n = B_c_n_opt_n;

      % Extract optimal phase shifts
        [V_max, max_eigenvalue_v] = max_eigVect(V_opt);
        v_k = sqrt(max_eigenvalue_v) * V_max;

        % Test candidate solutions
        cand = {exp(1j * angle(v_k)), exp(-1j * angle(v_k)), ...
                conj(exp(1j * angle(v_k))), conj(exp(-1j * angle(v_k)))};

        bestWSR = 0;
        for t = 1:4
            theta_test = diag(cand{t});
            WSR = Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                                g_b_all, f1_all, f2_all, alpha_n, ...
                                alpha_f, theta_test, J_t, J_r);
            if WSR > bestWSR
                bestWSR = WSR;
                Theta = theta_test;
            end
        end
    inter_cluster_interference_n = zeros(K, outer_iter);
    inter_cluster_interference_n_bst = zeros(K, outer_iter);
    inter_cluster_interference_f = zeros(K, outer_iter);
    inter_cluster_interference_n_b = zeros(K, outer_iter);


    obj_history = zeros(outer_iter, 1);
    far_history = zeros(K, outer_iter);
    near_history = zeros(K, outer_iter);

    far_channel = zeros(K, outer_iter);
    near_channel = zeros(K, outer_iter);   
     
        [WSR,R_n,R_f,R_c_n,A_f,A_n] = Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                    g_b_all, f1_all, f2_all, alpha_n, alpha_f, Theta, J_t, J_r);

            [inter_cluster_interference_near,inter_cluster_interference_near_bst,inter_cluster_interference_far,inter_cluster_interference_near_b,A_f,A_n] = sinr_terms(para,w_k,G_all_matrix, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, alpha_n, alpha_f, Theta,J_t,J_r);

        inter_cluster_interference_n(:,1) = inter_cluster_interference_near;
        inter_cluster_interference_n_bst(:,1) = inter_cluster_interference_near_bst;
        inter_cluster_interference_f(:,1) = inter_cluster_interference_far;
        inter_cluster_interference_n_b(:,1) = inter_cluster_interference_near_b;

        obj_history(1) = WSR;
        far_history(:,1) = R_f;
        near_history(:,1) = R_n;

        far_channel(:,1) = A_f;
        near_channel(:,1) = A_n;
    % Main optimization loop for DRIS
    outer_iter=outer_iter-1;
    for tau_2 = 1:outer_iter
        % Update channels
        % [g_1_all, g_2_all, g_b_all, f1_all, f2_all, decoding_order] = ensure_decoding_order(para, Theta, w_k, G_all_matrix, g_local, f_local, J_t, J_r);

        for c=1:K  
            % alpha_f(c) = para.alpha_k_f;
            % alpha_n(c) = para.alpha_k_n;        
            H_n{c} = g_1_all{c}' * J_r * Theta * J_t * G_all_matrix;
            H_f{c} = g_2_all{c}' * J_r * Theta * J_t * G_all_matrix;
            H_n_c{c} = g_b_all{c}' * J_r * Theta * J_t * f1_all{c} * G_all_matrix;
            H_f_c{c} = g_b_all{c}' * J_r * Theta * J_t * f2_all{c} * G_all_matrix;                     
        end

        disp(['======== DRIS Outer Iteration ', num2str(tau_2), ' ========']);

        if tau_2 > 1
        %  [alpha_n, alpha_f] = pac_opt_final(para, w_k, G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, Theta, J_t, J_r);
        %    disp(far_history);
        %    disp(near_history);
        end
        % [alpha_n, alpha_f] = pac_opt_final(para, w_k, G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, Theta, J_t, J_r);

        % Active BF optimization
        [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, ~, ~, ~, cvx_status] = ...
            active_BF_ndrs(para, H_n, H_f, H_n_c, H_f_c, A_n_prev, B_n_prev, ...
            A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, alpha_f, alpha_n, max_iter);

        A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
        A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
        A_c_prev_n = A_c_n_opt; B_c_prev_n = B_c_n_opt;


        % Extract beamforming vectors
        for k = 1:K
            [W_max, max_eigenvalue_w] = max_eigVect(W_opt(:, :, k));
            w_k(:, k) = sqrt(max_eigenvalue_w) * W_max;
        end

        % Passive BF optimization
        [V_opt, A_n_opt_n, B_n_opt_n, A_f_opt_n, B_f_opt_n, A_c_n_opt_n, B_c_n_opt_n, ~, ~, ~] = ...
            passive_bf_ndris_opt(para, w_k, G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, ...
            A_n_prev_n, B_n_prev_n, A_f_prev_n, B_f_prev_n, A_c_prev_n_n, B_c_prev_n_n, alpha_f, alpha_n, max_iter, 1, J_t, J_r);

        A_n_prev_n = A_n_opt_n; B_n_prev_n = B_n_opt_n; 
        A_f_prev_n = A_f_opt_n; B_f_prev_n = B_f_opt_n; 
        A_c_prev_n_n = A_c_n_opt_n; B_c_prev_n_n = B_c_n_opt_n;

        % Extract optimal phase shifts
        [V_max, max_eigenvalue_v] = max_eigVect(V_opt);
        v_k = sqrt(max_eigenvalue_v) * V_max;

        % Test candidate solutions
        cand = {exp(1j * angle(v_k)), exp(-1j * angle(v_k)), ...
                conj(exp(1j * angle(v_k))), conj(exp(-1j * angle(v_k)))};

        bestWSR = 0;
        for t = 1:4
            theta_test = diag(cand{t});
            WSR = Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                                g_b_all, f1_all, f2_all, alpha_n, ...
                                alpha_f, theta_test, J_t, J_r);
            if WSR > bestWSR
                bestWSR = WSR;
                Theta = theta_test;
            end
        end

        [WSR,R_n,R_f,R_c_n,A_f,A_n]= Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                    g_b_all, f1_all, f2_all, alpha_n, alpha_f, Theta, J_t, J_r);
            [inter_cluster_interference_near,inter_cluster_interference_near_bst,inter_cluster_interference_far,inter_cluster_interference_near_b,A_f,A_n] = sinr_terms(para,w_k,G_all_matrix, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, alpha_n, alpha_f, Theta,J_t,J_r)

        inter_cluster_interference_n(:,tau_2+1) = inter_cluster_interference_near;
        inter_cluster_interference_n_bst(:,tau_2+1) = inter_cluster_interference_near_bst;
        inter_cluster_interference_f(:,tau_2+1) = inter_cluster_interference_far;
        inter_cluster_interference_n_b(:,tau_2+1) = inter_cluster_interference_near_b;

        obj_history(tau_2+1) = WSR;
        % disp(size(R_n));
        far_history(:,tau_2+1) = R_f;
        near_history(:,tau_2+1) = R_n;

        far_channel(:,tau_2+1) = A_f;
        near_channel(:,tau_2+1) = A_n;
        [g_1_all, g_2_all, g_b_all, f1_all, f2_all, decoding_order] = ensure_decoding_order(para, Theta, w_k, G_all_matrix, g_local, f_local, J_t, J_r);
            

        
    end
end

%% Helper function for NDRIS system (Optimized J matrices)
function [inter_cluster_interference_n,inter_cluster_interference_n_bst,inter_cluster_interference_f,inter_cluster_interference_n_b,far_history,near_history,obj_history,alpha_n, alpha_f,w_k,R_n,R_f,R_c_n,far_channel,near_channel] = run_ndris_system(para, H_local , g_local, f_local, g_1_all, g_2_all, g_b_all, f1_all, f2_all, G_all_matrix, outer_iter, max_iter, max_feasible)
    
    K = para.K;
    N = para.N;
    M = para.M;
    
    % NDRIS: Optimized J matrices
    g_LOS_reshaped = reshape(g_local, para.N, []);
    J_r = design_J_r(g_LOS_reshaped);
    J_t = design_J_t(H_local);
    
    % Initialize Theta
    u = exp(1i * pi * (2 * rand(N, 1)));
    Theta = diag(u);
    
    % Initialize variables
    alpha_f = zeros(K,1);
    alpha_n = zeros(K,1);
    H_n = cell(1, K); H_f = cell(1, K);
    H_n_c = cell(1, K); H_f_c = cell(1, K);

    % Initial channel calculations
    for c=1:K
        alpha_f(c) = para.alpha_k_f;
        alpha_n(c) = para.alpha_k_n;             
        H_n{c} = g_1_all{c}' * J_r * Theta * J_t * G_all_matrix;
        H_f{c} = g_2_all{c}' * J_r * Theta * J_t * G_all_matrix;
        H_n_c{c} = g_b_all{c}' * J_r * Theta * J_t * f1_all{c} * G_all_matrix;
        H_f_c{c} = g_b_all{c}' * J_r * Theta * J_t * f2_all{c} * G_all_matrix;                     
    end
    
    w_k = mrt_beamforming(para, H_n);   
    
    [g_1_all, g_2_all, g_b_all, f1_all, f2_all, decoding_order] = ensure_decoding_order(para, Theta, w_k, G_all_matrix, g_local, f_local, J_t, J_r);

    % [alpha_n, alpha_f] = pac_opt_final(para, w_k, G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, Theta, J_t, J_r);

    % Initialize Taylor points
    A_n_prev_n = ones(K, 1); B_n_prev_n = ones(K, 1);
    A_f_prev_n = ones(K, 1); B_f_prev_n = ones(K, 1);
    A_c_prev_n_n = ones(K, 1); B_c_prev_n_n = ones(K, 1);
    A_n_prev = ones(K, 1); B_n_prev = ones(K, 1);
    A_f_prev = ones(K, 1); B_f_prev = ones(K, 1);
    A_c_prev_n = ones(K, 1); B_c_prev_n = ones(K, 1);

    % Find initial feasible points for NDRIS
    [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, ~, ~, ~] = ...
        find_feasible_ndris_active(para, Theta, w_k, G_all_matrix, g_1_all, g_2_all, ...
        g_b_all, f1_all, f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, ...
        A_c_prev_n, B_c_prev_n, alpha_f, alpha_n, max_feasible, 1, J_t, J_r);

    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
    A_c_prev_n = A_c_n_opt; B_c_prev_n = B_c_n_opt;
        for k = 1:K
            [W_max, max_eigenvalue_w] = max_eigVect(W_opt(:, :, k));
            w_k(:, k) = sqrt(max_eigenvalue_w) * W_max;
        end

    [V_opt, A_n_opt_n, B_n_opt_n, A_f_opt_n, B_f_opt_n, A_c_n_opt_n, B_c_n_opt_n, ~, ~, ~] = ...
        feasible_passive_soln(para, w_k, G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, ...
        A_n_prev_n, B_n_prev_n, A_f_prev_n, B_f_prev_n, A_c_prev_n_n, B_c_prev_n_n, alpha_f, alpha_n, J_t, J_r);

    A_n_prev_n = A_n_opt_n; B_n_prev_n = B_n_opt_n; 
    A_f_prev_n = A_f_opt_n; B_f_prev_n = B_f_opt_n; 
    A_c_prev_n_n = A_c_n_opt_n; B_c_prev_n_n = B_c_n_opt_n;

    obj_history = zeros(outer_iter, 1);

    far_history = zeros(K, outer_iter);
    near_history = zeros(K, outer_iter);

    far_channel = zeros(K, outer_iter);
    near_channel = zeros(K, outer_iter);  

    inter_cluster_interference_n = zeros(K, outer_iter);
    inter_cluster_interference_n_bst = zeros(K, outer_iter);
    inter_cluster_interference_f = zeros(K, outer_iter);
    inter_cluster_interference_n_b = zeros(K, outer_iter);

    [WSR,R_n,R_f,R_c_n,A_f,A_n] = Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                    g_b_all, f1_all, f2_all, alpha_n, alpha_f, Theta, J_t, J_r);
    [inter_cluster_interference_near,inter_cluster_interference_near_bst,inter_cluster_interference_far,inter_cluster_interference_near_b,A_f,A_n] = sinr_terms(para,w_k,G_all_matrix, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, alpha_n, alpha_f, Theta,J_t,J_r)

        obj_history(1) = WSR;
        far_history(:,1) = R_f;
        near_history(:,1) = R_n;

        far_channel(:,1) = A_f;
        near_channel(:,1) = A_n;
        inter_cluster_interference_n(:,1) = inter_cluster_interference_near;
        inter_cluster_interference_n_bst(:,1) = inter_cluster_interference_near_bst;
        inter_cluster_interference_f(:,1) = inter_cluster_interference_far;
        inter_cluster_interference_n_b(:,1) = inter_cluster_interference_near_b;
        
    % Main optimization loop for NDRIS
    outer_iter=outer_iter-1;
    for tau_2 = 1:outer_iter
        % Update channels
        % [g_1_all, g_2_all, g_b_all, f1_all, f2_all, decoding_order] = ensure_decoding_order(para, Theta, w_k, G_all_matrix, g_local, f_local, J_t, J_r);

        for c=1:K           
            H_n{c} = g_1_all{c}' * J_r * Theta * J_t * G_all_matrix;
            H_f{c} = g_2_all{c}' * J_r * Theta * J_t * G_all_matrix;
            H_n_c{c} = g_b_all{c}' * J_r * Theta * J_t * f1_all{c} * G_all_matrix;
            H_f_c{c} = g_b_all{c}' * J_r * Theta * J_t * f2_all{c} * G_all_matrix;                     
        end
        disp(['======== NDRIS Outer Iteration ', num2str(tau_2), ' ========']);

        if tau_2 > 1
            % [alpha_n, alpha_f] = pac_opt_final(para, w_k, G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, Theta, J_t, J_r);
        end

        % Active BF optimization
        [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, ~, ~, ~, ~] = ...
            active_BF_ndrs(para, H_n, H_f, H_n_c, H_f_c, A_n_prev, B_n_prev, ...
            A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, alpha_f, alpha_n, max_iter);

        A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
        A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
        A_c_prev_n = A_c_n_opt; B_c_prev_n = B_c_n_opt;

        % Extract beamforming vectors
        for k = 1:K
            [W_max, max_eigenvalue_w] = max_eigVect(W_opt(:, :, k));
            w_k(:, k) = sqrt(max_eigenvalue_w) * W_max;
        end

        % Passive BF optimization
        [V_opt, A_n_opt_n, B_n_opt_n, A_f_opt_n, B_f_opt_n, A_c_n_opt_n, B_c_n_opt_n, ~, ~, ~] = ...
            passive_bf_ndris_opt(para, w_k, G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, ...
            A_n_prev_n, B_n_prev_n, A_f_prev_n, B_f_prev_n, A_c_prev_n_n, B_c_prev_n_n, alpha_f, alpha_n, max_iter, 1, J_t, J_r);

        A_n_prev_n = A_n_opt_n; B_n_prev_n = B_n_opt_n; 
        A_f_prev_n = A_f_opt_n; B_f_prev_n = B_f_opt_n; 
        A_c_prev_n_n = A_c_n_opt_n; B_c_prev_n_n = B_c_n_opt_n;

        % Extract optimal phase shifts
        [V_max, max_eigenvalue_v] = max_eigVect(V_opt);
        v_k = sqrt(max_eigenvalue_v) * V_max;

        % Test candidate solutions
        cand = {exp(1j * angle(v_k)), exp(-1j * angle(v_k)), ...
                conj(exp(1j * angle(v_k))), conj(exp(-1j * angle(v_k)))};

        bestWSR = 0;
        for t = 1:4
            theta_test = diag(cand{t});
            WSR = Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                                g_b_all, f1_all, f2_all, alpha_n, ...
                                alpha_f, theta_test, J_t, J_r);
            if WSR > bestWSR
                bestWSR = WSR;
                Theta = theta_test;
            end
        end

        [WSR,R_n,R_f,R_c_n,A_f,A_n] = Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                    g_b_all, f1_all, f2_all, alpha_n, alpha_f, Theta, J_t, J_r);

    [inter_cluster_interference_near,inter_cluster_interference_near_bst,inter_cluster_interference_far,inter_cluster_interference_near_b,A_f,A_n] = sinr_terms(para,w_k,G_all_matrix, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, alpha_n, alpha_f, Theta,J_t,J_r)
            inter_cluster_interference_n(:,tau_2+1) = inter_cluster_interference_near;
        inter_cluster_interference_n_bst(:,tau_2+1) = inter_cluster_interference_near_bst;
        inter_cluster_interference_f(:,tau_2+1) = inter_cluster_interference_far;
        inter_cluster_interference_n_b(:,tau_2+1) = inter_cluster_interference_near_b;


        obj_history(tau_2+1) = WSR;
        far_history(:,tau_2+1) = R_f;
        near_history(:,tau_2+1) = R_n;

        far_channel(:,tau_2+1) = A_f;
        near_channel(:,tau_2+1) = A_n;

        [g_1_all, g_2_all, g_b_all, f1_all, f2_all, decoding_order] = ensure_decoding_order(para, Theta, w_k, G_all_matrix, g_local, f_local, J_t, J_r);

    end
end
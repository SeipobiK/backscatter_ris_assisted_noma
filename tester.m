clear all; clc;
addpath(genpath('/home/morolong/Documents/MATLAB/backscatter_ris_assisted_noma'));
rng(2022);

% Initialize parameters
para = para_init();
[BS_array, RIS_array] = generate_arrays(para);

% Constants
max_feasible = 25;
max_iter = 25;
outer_iter = para.outer_iter;
MC_MAX = para.MC_MAX;
K = para.K;
N = para.N;
M = para.M;

% Preallocate results
obj_history_all = zeros(outer_iter, MC_MAX);
obj_history_ac = zeros(outer_iter, MC_MAX);
rng_seeds = randi(1e6, MC_MAX, 1); % one seed per MC run

obj_history_dris = zeros(max_iter, MC_MAX);
obj_history_ndris = zeros(max_iter, MC_MAX);

% Create results directory if it doesn't exist
results_dir = 'results_passive';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end


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

% Precompute arrays once since they don't change across MC runs
[BS_array_par, RIS_array_par] = generate_arrays(para);
  
       
for mc = 1:MC_MAX
    % try
        fprintf('Monte Carlo Iteration %d\n', mc);
        
        % Set random seed for reproducibility
        % stream = RandStream('mt19937ar', 'Seed', rng_seeds(mc));
        % RandStream.setGlobalStream(stream);
        P_max=para.P_max;
        rng(rng_seeds(mc), 'twister');   % fix seed for this MC run



        % Initialize Theta(tau0)
        u = exp(1i * pi * (2 * rand(N, 1)));
        Theta = diag(u);
        Theta_1 = diag(u);

        
        % Generate channels for this MC run
        [H_local, g_local, f_local] = generate_channel(para, BS_array_par, RIS_array_par);
        
        
        % Precompute channel components as local variables
        g_1_all = cell(1, K);
        g_2_all = cell(1, K);
        g_b_all = cell(1, K);
        f1_all = cell(1, K);
        f2_all = cell(1, K);

        N_u=3*K;

        G_alll=cell(N, N_u);     
        g_LOS_reshaped = reshape(g_local, para.N, []);
      
        



        % Find Permutation matrices
        scal=para.scall;
        % 
        J_r=eye(N);
        J_t=eye(N);
        for i = 1:K
            g_1_all{i} = g_local(:, i, 1)*scal;
            g_2_all{i} = g_local(:, i, 2)*scal;
            g_b_all{i} = g_local(:, i, 3)*scal; 
            f1_all{i} = f_local(i, 1);
            f2_all{i} = f_local(i, 2);            
        end
        
        G_all_matrix = H_local*scal;

        g_LOS_reshaped=[g_1_all{1},g_2_all{1},g_b_all{1},g_1_all{2},g_2_all{2},g_b_all{2};];

        % J_r = design_J_r(g_LOS_reshaped);
        % J_t = design_J_t(G_all_matrix);
       
        G_RIS_users = zeros(K, N);
 
        % Initialize w_k(tau0) and pac using MRT
        alpha_f=zeros(K,1);
        alpha_n=zeros(K,1);
        H_n = cell(1, K); H_f = cell(1, K);
        H_n_c = cell(1, K); H_f_c = cell(1, K);


        for c=1:K
            % disp('========================Active BFFFF======================================');
              alpha_f(c)=para.alpha_k_f;
              alpha_n(c)=para.alpha_k_n;            
              H_n{c}  = g_1_all{c}'*J_r*Theta*J_t*G_all_matrix;
              H_f{c}  = g_2_all{c}'*J_r*Theta*J_t*G_all_matrix;

              H_n_c{c}  = g_b_all{c}'*J_r*Theta*J_t*f1_all{c}*G_all_matrix;
              H_f_c{c}  = g_b_all{c}'*J_r*Theta*J_t*f2_all{c}*G_all_matrix;                     
        end
        w_k = mrt_beamforming(para, H_n);   

        % [alpha_n, alpha_f] = pac_opt_final(para,w_k,G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, Theta,J_t,J_r);

 
                       [WSR,~,~,~,~,~] = Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                       g_b_all, f1_all, f2_all, alpha_n, ...
                       alpha_f, Theta,J_t,J_r);


                       disp(['Initial WSR at MC ', num2str(mc), ': ', num2str(WSR)]);
                       

                        A_n_prev_n = ones(K, 1); 
                        B_n_prev_n = ones(K, 1);
                        A_f_prev_n = ones(K, 1) ; 
                        B_f_prev_n = ones(K, 1);
                        A_c_prev_n_n = ones(K, 1); 
                        B_c_prev_n_n = ones(K, 1);

                        A_n_prev = ones(K, 1); 
                        B_n_prev = ones(K, 1);
                        A_f_prev = ones(K, 1); 
                        B_f_prev = ones(K, 1);
                        A_c_prev_n = ones(K, 1); 
                        B_c_prev_n = ones(K, 1);


                       [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, ~, ~, ~] = ...
                            find_feasible_ndris_active(para, Theta,w_k,G_all_matrix, g_1_all, g_2_all, ...
                            g_b_all, f1_all, f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, ...
                            A_c_prev_n, B_c_prev_n,alpha_f,alpha_n, max_feasible, mc,J_t,J_r);


                        


                        % Update Taylor points
                        A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                        A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                        A_c_prev_n = A_c_n_opt; B_c_prev_n = B_c_n_opt;


                        djkfjf

                     % 
                     % % extract the solution(Update Active BF vector for next iteration)
                     for k = 1:K
                         [W_max, max_eigenvalue_w] = max_eigVect(W_opt(:, :, k));
                         w_k(:, k) = sqrt(max_eigenvalue_w) * W_max;
                     end
                        
    
                        [V_opt, A_n_opt_n, B_n_opt_n, A_f_opt_n, B_f_opt_n, A_c_n_opt_n, B_c_n_opt_n, obj_prev, cvx_status, ~] = ...
                            feasible_passive_soln(para,w_k,G_all_matrix, g_1_all,...
                             g_2_all,g_b_all,f1_all,f2_all,A_n_prev_n, B_n_prev_n, A_f_prev_n, B_f_prev_n,  A_c_prev_n_n, B_c_prev_n_n,alpha_f,alpha_n,J_t,J_r);
                    
                            A_n_prev_n = A_n_opt_n; B_n_prev_n = B_n_opt_n; 
                            A_f_prev_n = A_f_opt_n; B_f_prev_n = B_f_opt_n; 
                            A_c_prev_n_n = A_c_n_opt_n; B_c_prev_n_n = B_c_n_opt_n;


                   % % % Extract optimal phase shifts(DRIS)
                    [V_max, max_eigenvalue_v] = max_eigVect(V_opt);
                     v_k= sqrt(max_eigenvalue_v) * V_max;
                     v_opt=exp(1j * angle(v_k));

                    % Test candidate solutions
                    cand = {exp(1j * angle(v_k)), exp(-1j * angle(v_k)), ...
                            conj(exp(1j * angle(v_k))), conj(exp(-1j * angle(v_k)))};

                        Theta = diag(cand{2});


                            

        obj_history_ac_local = zeros(outer_iter, 1);
                        

        for tau_2 = 1:para.outer_iter

               for c=1:K
                     % disp('========================Active BFFFF======================================');           
                     H_n{c}  = g_1_all{c}'*J_r*Theta*J_t*G_all_matrix;
                     H_f{c}  = g_2_all{c}'*J_r*Theta*J_t*G_all_matrix;

                     H_n_c{c}  = g_b_all{c}'*J_r*Theta*J_t*f1_all{c}*G_all_matrix;
                     H_f_c{c}  = g_b_all{c}'*J_r*Theta*J_t*f2_all{c}*G_all_matrix;                     
               end

            %    [alpha_n, alpha_f] = pac_opt_final(para,w_k,G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, Theta,J_t,J_r);


               % %  % % % ======================================Active BF=================
                    [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, ...
                     obj_history_active, obj_history_, ~, ~] = active_BF_ndrs(para, H_n, ...
                     H_f,H_n_c,H_f_c, A_n_prev, B_n_prev, ...
                     A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n,alpha_f,alpha_n, max_iter);

                     % Update Taylor points and extract beamforming vectors
                     A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                     A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                     A_c_prev_n = A_c_n_opt; B_c_prev_n = B_c_n_opt;




                     



                    % extract the solution(Update Active BF vector for next iteration)
                     for k = 1:K
                         [W_max, max_eigenvalue_w] = max_eigVect(W_opt(:, :, k));
                         w_k(:, k) = sqrt(max_eigenvalue_w) * W_max;
                     end
            % % % %=====================Active BF optimization================================


     
            %=====================Passive BF optimization====================================            
               
                % % % Passive BF optimization(NDRIS)
                [V_opt, A_n_opt_n, B_n_opt_n, A_f_opt_n, B_f_opt_n, A_c_n_opt_n, B_c_n_opt_n,obj_history_n, obj_history_mc_n, ~] =passive_bf_ndris_opt(para,w_k,G_all_matrix, g_1_all,...
                     g_2_all,g_b_all,f1_all,f2_all, A_n_prev_n, B_n_prev_n, A_f_prev_n, B_f_prev_n, ...
                     A_c_prev_n_n, B_c_prev_n_n,alpha_f,alpha_n,max_iter,mc,J_t,J_r);



                 % Update Taylor points
                    A_n_prev_n = A_n_opt_n; B_n_prev_n = B_n_opt_n; 
                    A_f_prev_n = A_f_opt_n; B_f_prev_n = B_f_opt_n; 
                    A_c_prev_n_n = A_c_n_opt_n; B_c_prev_n_n = B_c_n_opt_n;

                   % % 
                   % % % Extract optimal phase shifts(DRIS)
                     [V_max, max_eigenvalue_v] = max_eigVect(V_opt);
                     v_k= sqrt(max_eigenvalue_v) * V_max;
                     v_opt=exp(1j * angle(v_k));

                    % Test candidate solutions
                    cand = {exp(1j * angle(v_k)), exp(-1j * angle(v_k)), ...
                            conj(exp(1j * angle(v_k))), conj(exp(-1j * angle(v_k)))};

                    bestWSR = 0;
                    mm=zeros(outer_iter);
                    for t = 1:4
                        theta_test = diag(cand{t});
                        WSR = Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                                            g_b_all, f1_all, f2_all, alpha_n, ...
                                            alpha_f, theta_test,J_t,J_r);
                        if WSR > bestWSR
                            bestWSR = WSR;
                            Theta = theta_test;
                            mm(tau_2)=t;
                        end
                    end
                    
                 [WSR,~,~,~,~,~] = Compute_WSR_NDRIS(para, w_k, G_all_matrix, g_1_all, g_2_all, ...
                            g_b_all, f1_all, f2_all, alpha_n, ...
                            alpha_f, Theta,J_t,J_r);

                
                   
                    disp(['WSR at Tau   ','At tau_1 :',num2str(tau_2),' ',num2str(norm(Theta, 'fro')^2 )]);
                    disp(['WSR at Tau   ',num2str(WSR)]);  

                    obj_history_ac_local(tau_2) = WSR;

                    % [g_1_all, g_2_all, g_b_all, f1_all, f2_all, decoding_order] = ensure_decoding_order(para, Theta, w_k, G_all_matrix, g_local, f_local, J_t, J_r);




        end
        obj_history_ac(:, mc) = obj_history_ac_local;

        


       % obj_history_ac(:, mc) = obj_history_n;
       % obj_history_dris(:, mc) = obj_history;
      % disp(obj_history_n);
    % % 
    % catch ME
    %    obj_history_ndris(:, mc) = NaN;
    %    obj_history_dris(:, mc) = NaN;    
    % end
end 

disp(obj_history_ac);

% %% Step 1: Identify valid MC runs
valid_dris  = all(isfinite(obj_history_ac) & obj_history_ac ~= 0, 1);   % 1 × MC_MAX
valid_ndris = all(isfinite(obj_history_ac) & obj_history_ac ~= 0, 1); % 1 × MC_MAX

% % % Only keep MC runs valid in both systems
valid_MC = valid_dris & valid_ndris;

% % Filter history matrices
obj_history_dris  = obj_history_ac(:, valid_MC);
obj_history_ndris = obj_history_ac(:, valid_MC);

MC_valid = size(obj_history_dris);  % number of valid MC runs

% %% Step 2: Compute mean per iteration
avg_dris_iter  = mean(obj_history_dris, 2);   % mean across valid MC runs
avg_ndris_iter = mean(obj_history_ndris, 2);
% disp(valid_dris);
% disp(valid_ndris);
% disp(valid_MC);




    x = 1:outer_iter;
    figure;
    plot(x, avg_dris_iter, '-o', 'DisplayName','P=1W/30dBm, M=10, N=4' , 'LineWidth', 2, 'MarkerSize', 6);
    hold on;
    plot(x, avg_ndris_iter, '-s', 'DisplayName','P=1W/30dBm, M=10, N=4 ', 'LineWidth', 2, 'MarkerSize', 6);
    hold off;
    legend('Location', 'southeast');

    xlabel('number of iterations');
    ylabel('Weighted Sum Rate (bps/Hz)');
    title('');
    grid on;
    xlim([1, outer_iter]);
%     png_file = fullfile(results_dir, ['PassiveBF_sims', timestamp, '.png']);
disp(['Valid Mcs:  ',num2str(sum(valid_MC))]);

%     saveas(gcf, png_file);

% After Monte Carlo loop finishes
results.obj_history_ac = obj_history_ac;          % from your AO/SCA loop
results.para = para;                        % system parameters

% Add timestamp for unique save
timestamp = datestr(now,'yyyymmdd_HHMMSS');
filename = sprintf('c_ndrisresults_M%d_N%d_%s.mat', para.MC_MAX, para.N, timestamp);

save(filename, 'results');
fprintf('Simulation results saved in %s\n', filename);
toc
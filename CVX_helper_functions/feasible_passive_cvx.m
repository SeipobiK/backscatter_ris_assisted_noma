% function [V_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, status] = feasible_passive_cvx(para, w_k, G_all, g_1_all, ...
%     g_2_all, g_b_all, f1_all, f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, alpha_f, alpha_n, J_t, J_r)
   
%    numClusters = para.K;
%    N = para.N;
%    R_n_min = para.R_min_n;
%    R_f_min = para.R_min_f;
%    R_c_min = para.R_c_min;
%    eta_k = para.eta;
   
%    % Precompute all channel matrices outside CVX
%    H_n = cell(numClusters, 1);
%    H_f = cell(numClusters, 1);
%    H_n_c = cell(numClusters, 1);
%    H_f_c = cell(numClusters, 1);
   
%    % Precompute interference terms
%    interference_near = cell(numClusters, 1);
%    interference_far = cell(numClusters, 1);
%    interference_near_b = cell(numClusters, 1);
   
%    for c = 1:numClusters
%        % Compute channel matrices
%        H_n_temp = diag(g_1_all{c}' * J_r) * J_t * G_all * w_k(:, c);
%        H_f_temp = diag(g_2_all{c}' * J_r) * J_t * G_all * w_k(:, c);
%        H_n_c_temp = diag(g_b_all{c}' * J_r) * J_t * f1_all{c} * G_all * w_k(:, c);
%        H_f_c_temp = diag(g_b_all{c}' * J_r) * J_t * f2_all{c} * G_all * w_k(:, c);
       
%        % Store for later use
%        H_n{c} = H_n_temp;
%        H_f{c} = H_f_temp;
%        H_n_c{c} = H_n_c_temp;
%        H_f_c{c} = H_f_c_temp;
       
%        % Precompute interference terms
%        interference_near{c} = 0;
%        interference_far{c} = 0;
%        interference_near_b{c} = 0;
       
%        for j = 1:numClusters
%            if j ~= c
%                H_n_j = diag(g_1_all{c}' * J_r) * J_t * G_all * w_k(:, j);
%                H_f_j = diag(g_2_all{c}' * J_r) * J_t * G_all * w_k(:, j);
%                H_n_b_j = diag(g_b_all{c}' * J_r) * J_t * f1_all{c} * G_all * w_k(:, j);
               
%                interference_near{c} = interference_near{c} + (H_n_j * H_n_j');
%                interference_far{c} = interference_far{c} + (H_f_j * H_f_j');
%                interference_near_b{c} = interference_near_b{c} + (H_n_b_j * H_n_b_j');
%            end
%        end
%    end
   
%    % Precompute Taylor approximation constants
%    log2_e = log2(exp(1));
%    taylor_const_far = zeros(numClusters, 2);
%    taylor_const_near = zeros(numClusters, 2);
%    taylor_const_back = zeros(numClusters, 2);
   
%    for c = 1:numClusters
%        denom_far = 1 + A_f_prev(c) * B_f_prev(c);
%        denom_near = 1 + A_n_prev(c) * B_n_prev(c);
%        denom_back = 1 + A_c_prev_n(c) * B_c_prev_n(c);
       
%        taylor_const_far(c, 1) = log2_e / (A_f_prev(c) * denom_far);
%        taylor_const_far(c, 2) = log2_e / (B_f_prev(c) * denom_far);
       
%        taylor_const_near(c, 1) = log2_e / (A_n_prev(c) * denom_near);
%        taylor_const_near(c, 2) = log2_e / (B_n_prev(c) * denom_near);
       
%        taylor_const_back(c, 1) = log2_e / (A_c_prev_n(c) * denom_back);
%        taylor_const_back(c, 2) = log2_e / (B_c_prev_n(c) * denom_back);
%    end

%    cvx_begin quiet
%        cvx_solver mosek
%         % 
%         % % Optimized solver settings
%         % cvx_solver_settings('MSK_IPAR_NUM_THREADS', 4);
%         % cvx_solver_settings('MSK_IPAR_INTPNT_BASIS', 0);
%         % cvx_solver_settings('MSK_DPAR_INTPNT_TOL_PFEAS', 1e-8);
%         % cvx_solver_settings('MSK_DPAR_INTPNT_TOL_DFEAS', 1e-8);
%         % cvx_solver_settings('MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-6);
       
%        % Variable declarations
%        variable V(N, N) Hermitian semidefinite
%        variable A_n(numClusters) nonnegative
%        variable B_n(numClusters) nonnegative
%        variable A_f(numClusters) nonnegative
%        variable B_f(numClusters) nonnegative
%        variable A_c_n(numClusters) nonnegative
%        variable B_c_n(numClusters) nonnegative
%        variable R_n(numClusters) nonnegative
%        variable R_f(numClusters) nonnegative
%        variable R_c_n(numClusters) nonnegative
%        variable delta_p nonnegative

%        minimize(delta_p)

%        subject to
%            for c = 1:numClusters
%                % Precomputed Taylor approximations
%                base_far = log2(1 + 1/(A_f_prev(c) * B_f_prev(c)));
%                base_near = log2(1 + 1/(A_n_prev(c) * B_n_prev(c)));
%                base_back = log2(1 + 1/(A_c_prev_n(c) * B_c_prev_n(c)));
               
%                % Rate constraints
%                R_f(c) <= base_far - taylor_const_far(c, 1) * (A_f(c) - A_f_prev(c)) - ...
%                          taylor_const_far(c, 2) * (B_f(c) - B_f_prev(c)) + delta_p;
               
%                R_n(c) <= base_near - taylor_const_near(c, 1) * (A_n(c) - A_n_prev(c)) - ...
%                         taylor_const_near(c, 2) * (B_n(c) - B_n_prev(c)) + delta_p;
               
%                R_c_n(c) <= base_back - taylor_const_back(c, 1) * (A_c_n(c) - A_c_prev_n(c)) - ...
%                           taylor_const_back(c, 2) * (B_c_n(c) - B_c_prev_n(c)) + delta_p;
               
%                % Minimum rate constraints
%                R_f(c) >= R_f_min - delta_p;
%                R_n(c) >= R_n_min - delta_p;
%                R_c_n(c) >= R_c_min - delta_p;
               
%                % Main constraints using precomputed values
%                inv_pos(A_n(c)) <= real(trace(V * (H_n{c} * H_n{c}'))) * alpha_n(c) + delta_p;
               
%                B_n(c) >= real(trace(V * interference_near{c})) + ...
%                        real(trace(V * (H_n_c{c} * H_n_c{c}'))) * eta_k + para.noise - delta_p;
               
%                inv_pos(A_f(c)) <= real(trace(V * (H_f{c} * H_f{c}'))) * alpha_f(c) + delta_p;
               
%                B_f(c) >= real(trace(V * interference_far{c})) + ...
%                        real(trace(V * (H_f{c} * H_f{c}'))) * alpha_n(c) + ...
%                        real(trace(V * (H_f_c{c} * H_f_c{c}'))) * eta_k + para.noise - delta_p;
               
%                inv_pos(A_c_n(c)) <= real(trace(V * (H_n_c{c} * H_n_c{c}'))) * eta_k + delta_p;
               
%                B_c_n(c) >= real(trace(V * interference_near_b{c})) + para.noise - delta_p;
%            end
           
%            % Unit modulus constraints with relaxation
%            for m = 1:N
%                V(m, m) == 1 + delta_p;
%                % V(m, m) >= 1 - delta_p;
%            end
           
%            % PSD constraint
%           for k = 1:numClusters
%             V + delta_p * eye(N) == hermitian_semidefinite(N);
%           end           
%    cvx_end

%    % Return results
%    obj_prev = cvx_optval;
%    A_n_opt = A_n;
%    B_n_opt = B_n;
%    A_f_opt = A_f;
%    B_f_opt = B_f;
%    A_c_n_opt = A_c_n;
%    B_c_n_opt = B_c_n;
%    V_opt = V;
%    status = cvx_status;
% end










function [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = feasible_passive_cvx(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,alpha_f,alpha_n,J_t,J_r)

   numClusters = para.K; % Number of clusters
   N = para.N; % Number of BS antennas
   M = para.M; % Number of BS antennas
   R_n_min = para.R_min_n; % Minimum rate for near user
   R_f_min = para.R_min_f; % Minimum rate for far user
   R_c_min = para.R_c_min; % Minimum rate for backscatter user
   eta_k = para.eta; % Backscatter coefficient
   para.P_max = para.P_max;

   H_n = cell(numClusters,1); H_f = cell(numClusters,1);  
   H_n_c = cell(numClusters,1); H_f_c = cell(numClusters,1); 
   for c=1:numClusters
         H_n{c}  = diag(g_1_all{c}'*J_r)*J_t*G_all*w_k(:, c);
         H_f{c}  = diag(g_2_all{c}'*J_r)*J_t*G_all*w_k(:, c);

         H_n_c{c}  = diag(g_b_all{c}'*J_r)*J_t*f1_all{c}*G_all*w_k(:, c);
         H_f_c{c}  = diag(g_b_all{c}'*J_r)*J_t*f2_all{c}*G_all*w_k(:, c);
   end

   cvx_begin quiet sdp
       cvx_solver mosek
       % cvx_precision medium
       %  % cvx_precision high
       %  cvx_solver_settings( ...
       %      'MSK_DPAR_INTPNT_TOL_PFEAS', 1e-13, ...
       %      'MSK_DPAR_INTPNT_TOL_DFEAS', 1e-13, ...
       %      'MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-13 ...
       %  );

       %  (59h)
       variable V(N,N) Hermitian semidefinite  
       variable A_n(numClusters) nonnegative % Slack variable for near users
       variable B_n(numClusters) nonnegative % Slack variable for near users
       variable A_f(numClusters) nonnegative % Slack variable for far users
       variable B_f(numClusters) nonnegative % Slack variable for far users
       variable A_c_f(numClusters) nonnegative % Slack variable for backscatter devices at far user
       variable B_c_f(numClusters) nonnegative % Slack variable for backscatter devices at far user
       variable A_c_n(numClusters) nonnegative % Slack variable for backscatter devices at near user
       variable B_c_n(numClusters) nonnegative % Slack variable for backscatter devices at near user
       variable R_n(numClusters)  nonnegative% Slack variable for backscatter devices at near user
       variable R_f(numClusters)  nonnegative% Slack variable for backscatter devices at near user
       variable R_c_n(numClusters)  nonnegative% Slack variable for backscatter devices at near user
    %    variable delta_p% Slack variable for backscatter devices at near user
       variable delta_p nonnegative  % This should be a variable, not an expression

       expressions taylor_approx_far(numClusters, 1) taylor_approx_n(numClusters, 1) taylor_approx_backscatter_n(numClusters, 1)





       % Objective function: Maximize weighted sum rate (59a)
       minimise delta_p

        subject to
           for c = 1:numClusters

                    taylor_approx_far(c) = log2(1 + 1 ./ (A_f_prev(c) * B_f_prev(c))) -  ...
                    (log2(exp(1)) * 1 ./ (A_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (A_f(c) - A_f_prev(c)) - ...
                    (log2(exp(1)) * 1 ./ (B_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (B_f(c) - B_f_prev(c));


                    taylor_approx_n(c) = log2(1 + 1 ./ (A_n_prev(c) * B_n_prev(c))) -  ...
                    (log2(exp(1)) * 1 ./ (A_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (A_n(c) - A_n_prev(c)) - ...
                    (log2(exp(1)) * 1 ./ (B_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (B_n(c) - B_n_prev(c));



                    taylor_approx_backscatter_n(c)= log2(1 + 1 ./ (A_c_prev_n(c) * B_c_prev_n(c))) - ...
                    (log2(exp(1)) *1 ./  (A_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (A_c_n(c) - A_c_prev_n(c)) - ...
                    (log2(exp(1)) * 1 ./  (B_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (B_c_n(c) - B_c_prev_n(c)); 



                        R_f(c)- taylor_approx_far(c)<=delta_p;
                        R_n(c) - taylor_approx_n(c) <=delta_p;
                        R_c_n(c) - taylor_approx_backscatter_n(c)<=delta_p;

                       % (54b) and (54c)
                        delta_p>= R_f_min-R_f(c);  
                        delta_p>= R_n_min-R_n(c);
                        delta_p>= R_c_min- R_c_n(c);

                   inter_cluster_interference_near = 0;
                   inter_cluster_interference_far = 0;
                   inter_cluster_interference_near_b=0;
                %    inter_cluster_interference_far_b=0;

                   for j = 1:numClusters
                       if j ~= c

                           inter_cluster_interference_near = inter_cluster_interference_near + ...
                               real(trace(V * (diag(g_1_all{c}'*J_r)*J_t*G_all*w_k(:, j)) * (diag(g_1_all{c}'*J_r)*J_t*G_all*w_k(:, j))'));

                           inter_cluster_interference_far= inter_cluster_interference_far + ...
                               real(trace(V * (diag(g_2_all{c}'*J_r)*J_t*G_all*w_k(:, j)) * (diag(g_2_all{c}'*J_r)*J_t*G_all*w_k(:, j))')); 

                           inter_cluster_interference_near_b= inter_cluster_interference_near_b + ...
                               real(trace(V* (diag(g_b_all{c}'*J_r)*J_t*G_all*f1_all{c}*w_k(:, j)) * (diag(g_b_all{c}'*J_r)*J_t*G_all*f1_all{c}*w_k(:, j))'))*eta_k; 

                        %    inter_cluster_interference_far_b= inter_cluster_interference_far_b + ...
                        %     real(trace(V* (diag(g_b_all{c}'*J_r)*J_t*G_all*f2_all{c}*w_k(:, j)) * (diag(g_b_all{c}'*J_r)*J_t*G_all*f2_all{c}*w_k(:, j))'))*eta_k; 
                       end

                   end

                       % Define slack variables based on cascaded channel

                   delta_p >= inv_pos(A_n(c))-real(trace(V * H_n{c} * H_n{c}')) * alpha_n(c); % (59b)

                    % inter cluster interference  + backscatter interference + noise power
                   delta_p>=inter_cluster_interference_near + inter_cluster_interference_near_b+...
                            real(trace(V * H_n_c{c} * H_n_c{c}')) * eta_k + para.noise- B_n(c) ;  %% (59c)


                   delta_p>= inv_pos(A_f(c)) - real(trace(V * H_f{c} * H_f{c}')) * alpha_f(c); %% (59d)

                   % %% (59e)
                   delta_p>= inter_cluster_interference_far + ...
                           real(trace(V * H_f{c} * H_f{c}'))  * alpha_n(c)  + para.noise-B_f(c) ;

                   %% (59f)  
                   delta_p>=inv_pos(A_c_n(c)) - real(trace(V * H_n_c{c} * H_n_c{c}')) * eta_k;

                    %% (50g)  
                   delta_p >= inter_cluster_interference_near + inter_cluster_interference_near_b + para.noise- B_c_n(c) ;
                   %% (59j)

           end
        %    diag(V) == 1 + delta_p; 

          for m=1:N
            V(m,m) == 1 + delta_p;
          end

          for k = 1:numClusters
            V + delta_p * eye(N) == hermitian_semidefinite(N);
          end
   cvx_end

       obj_prev = cvx_optval;
       A_n_opt = A_n;
       B_n_opt = B_n;
       A_f_opt = A_f;
       B_f_opt = B_f;
       A_c_n_opt = A_c_n;
       B_c_n_opt = B_c_n;
       V_opt = V;
       status = cvx_status;
end
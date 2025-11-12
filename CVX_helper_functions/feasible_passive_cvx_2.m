function [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = feasible_passive_cvx_2(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,alpha_f,alpha_n,J_t,J_r)
   
   numClusters = para.K; % Number of clusters
   N = para.N; % Number of BS antennas
   M = para.M; % Number of BS antennas
   R_n_min = para.R_min_n; % Minimum rate for near user
   R_f_min = para.R_min_f; % Minimum rate for far user
   R_c_min = para.R_c_min; % Minimum rate for backscatter user
   eta_k = para.eta; % Backscatter coefficient
   para.P_max = para.P_max;
   K=para.K;

   H_n = cell(numClusters,1); H_f = cell(numClusters,1);  
   for c=1:numClusters
         H_n{c}  = diag(g_1_all{c}'*J_r)*J_t*G_all*w_k(:, c);
         H_f{c}  = diag(g_2_all{c}'*J_r)*J_t*G_all*w_k(:, c);
         H_n_{c}  = diag(g_1_all{c}'*J_r)*J_t*G_all*w_k(:, c);
         H_f_{c}  = diag(g_2_all{c}'*J_r)*J_t*G_all*w_k(:, c);

         disp(['With J Near QWE at : ',num2str(c) ,num2str(trace(H_n_{c}*H_n_{c}'))]);
         disp(['With J Far : ',num2str(c),num2str(trace(H_f_{c}*H_f_{c}'))]);
   end

   cvx_begin quiet sdp
       cvx_solver mosek


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
       

        minimise delta_p

        subject to
            for c = 1:numClusters
               
                       R_f(c) <= log2(1 + 1 ./ (A_f_prev(c) * B_f_prev(c))) -  ...
                       (log2(exp(1)) * 1 ./ (A_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (A_f(c) - A_f_prev(c)) - ...
                       (log2(exp(1)) * 1 ./ (B_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (B_f(c) - B_f_prev(c))+delta_p;

                       % (53b) for i=n
                       R_n(c) <= log2(1 + 1 ./ (A_n_prev(c) * B_n_prev(c))) -  ...
                       (log2(exp(1)) * 1 ./ (A_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (A_n(c) - A_n_prev(c)) - ...
                       (log2(exp(1)) * 1 ./ (B_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (B_n(c) - B_n_prev(c))+delta_p;

                       % (54b) and (54c)
                        delta_p>= R_f_min-R_f(c);  
                        delta_p>= R_n_min-R_n(c);

                        
                   inter_cluster_interference_near = 0;
                   inter_cluster_interference_far = 0;
                   for j = 1:numClusters
                       if j ~= c
                           inter_cluster_interference_near = inter_cluster_interference_near + ...
                               real(trace(V * (diag(g_1_all{c}'*J_r)*J_t*G_all*w_k(:, j)) * (diag(g_1_all{c}'*J_r)*J_t*G_all*w_k(:, j))'));

                           inter_cluster_interference_far= inter_cluster_interference_far + ...
                               real(trace(V * (diag(g_2_all{c}'*J_r)*J_t*G_all*w_k(:, j)) * (diag(g_2_all{c}'*J_r)*J_t*G_all*w_k(:, j))')); 
                       end

                   end
    

                   % Define slack variables based on cascaded channel
                   delta_p >= inv_pos(A_n(c))-real(trace(V * H_n{c} * H_n{c}')) * alpha_n(c); % (59b)

                   % inter cluster interference  + backscatter interference + noise power
                   delta_p>=inter_cluster_interference_near + para.noise- B_n(c) ;  %% (59c)

                       
                   delta_p>= inv_pos(A_f(c)) - real(trace(V * H_f{c} * H_f{c}')) * alpha_f(c); %% (59d)

                   %% (59e)
                   delta_p>= inter_cluster_interference_far + ...
                           real(trace(V * H_f{c} * H_f{c}'))  * alpha_n(c) + para.noise-B_f(c) ;

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
       A_c_n_opt = 0;
       B_c_n_opt = 0;
       V_opt = V;
       status = cvx_status;
end
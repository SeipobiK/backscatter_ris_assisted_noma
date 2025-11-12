function [alpha_n, alpha_f] = pac_opt_final(para,w_k,G_all_matrix, g_1_all, g_2_all, g_b_all, f1_all, f2_all, Theta,J_t,J_r)
    numClusters = para.K; % Number of cluster
    alpha_f = zeros(numClusters,1);
    alpha_n = zeros(numClusters,1);

          numClusters = para.K; % Number of 
          K = para.K;
          eta_k=para.eta;
          N=para.N; % Number of reflecting elements at RIS
          A_n = zeros(numClusters, 1); % Initialize A_n vector
          B_n = zeros(numClusters, 1); % Initialize B_n vector
          A_f = zeros(numClusters, 1); % Initialize A_f vector
          B_f = zeros(numClusters, 1); % Initialize B_f vector
          % A_c_n = zeros(numClusters, 1); % Initialize A_c_n vector
          % B_c_n = zeros(numClusters, 1); % Initialize B_c_n vector
          % R_n=zeros(numClusters, 1); % Initialize R_n vector
          % R_f=zeros(numClusters, 1); % Initialize R_f vector    
          % R_c_n=zeros(numClusters, 1); % Initialize R_c_n vector
          noise = para.noise;  % Noise power
        H_n = cell(1, K); H_f = cell(1, K);
        H_n_c = cell(1, K); H_f_c = cell(1, K);

        for c=1:K
            % disp('========================Active BFFFF======================================');
              alpha_f(c)=para.alpha_k_f;
              alpha_n(c)=para.alpha_k_n;            
              H_n{c}  = g_1_all{c}'*J_r*Theta*J_t*G_all_matrix*w_k(:, c);
              H_f{c}  = g_2_all{c}'*J_r*Theta*J_t*G_all_matrix*w_k(:, c);

              H_n_c{c}  = g_b_all{c}'*J_r*Theta*J_t*f1_all{c}*G_all_matrix*w_k(:, c);
              H_f_c{c}  = g_b_all{c}'*J_r*Theta*J_t*f2_all{c}*G_all_matrix*w_k(:, c);                     
        end




           for c = 1:numClusters
                    inter_cluster_interference_near = 0;
                    inter_cluster_interference_far = 0;
                    inter_cluster_interference_near_b=0;
                    
                    for j = 1:numClusters
                        if j ~= c
                            % Near user inter cluster interference  
                            inter_cluster_interference_near = inter_cluster_interference_near + ...
                                abs(g_1_all{c}'*J_r*Theta*J_t*G_all_matrix*w_k(:, j)).^2;

                            inter_cluster_interference_far= inter_cluster_interference_far + ...
                               abs(g_2_all{c}'*J_r*Theta*J_t*G_all_matrix*w_k(:, j)).^2;

                           inter_cluster_interference_near_b= inter_cluster_interference_near_b + ...
                                abs(g_b_all{c}'*J_r*Theta*J_t*f1_all{c}*G_all_matrix*w_k(:, j)).^2;
                            
                        end

                    end
                    eta_k=0;

                    A_n(c)= abs(H_n{c}).^2;

                    B_n(c) =inter_cluster_interference_near + ...
                            abs(H_n_c{c}).^2 * eta_k + noise;
        
                    A_f(c) = abs( H_f{c}).^2;

                    B_f(c) = inter_cluster_interference_far + ...
                            abs(H_f_c{c}).^2  * eta_k + noise;

                    

                    % Compute SINR-based fraction
                    gamma_f= A_f(c) / B_f(c);
                    gamma_n= A_n(c) / B_n(c);

                            % Minimum rate requirements
                    r_min = 2^para.R_min_f - 1;
                    
                    % NOMA Power Allocation (assuming fixed: near=strong, far=weak)
                    % Far user (weak) gets more power, Near user (strong) gets less power
                    alpha_n(c) = (r_min / (1 + r_min))*(1 + 1/gamma_f);
                    alpha_f(c) = 1 - alpha_n(c);
                    
                    % Verify decoding order is maintained
                    if gamma_n < gamma_f
                        warning(['Cluster ', num2str(c), ': Decoding order violated! Gamma_n=', ...
                                 num2str(gamma_n), ', Gamma_f=', num2str(gamma_f)]);
                    end
                    
                    % disp(['Cluster ', num2str(c), ': Gamma_n=', num2str(gamma_n), ...
                    %       ', Gamma_f=', num2str(gamma_f), ', Alpha_n=', num2str(alpha_n(c)), ...
                    %       ', Alpha_f=', num2str(alpha_f(c))]);
                    % 
                    % 
                    disp(gamma_n);

                    disp(['near user gamma  :',num2str(gamma_f),'   R min  : ',num2str(r_min),' PAC far user  :',num2str(alpha_f(c))])
                    alpha_f(c)=0.1;
                    alpha_n(c)=0.9;

                    
                              
           end

end

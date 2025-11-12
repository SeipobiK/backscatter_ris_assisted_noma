function [WSR,R_n,R_f,R_c_n,A_n,B_n] = Compute_WSR_NDRIS(para,w_k,G_all_matrix, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, alpha_n, alpha_f, Theta,J_t,J_r)
         % power= abs(w_k(:,1))^2 + abs(w_k(:,2))^2;
         % disp(['Total power in Compute WSR: ', num2str(power)]);         
          numClusters = para.K; % Number of 
          K = para.K;
          eta_k=para.eta;
          N=para.N; % Number of reflecting elements at RIS
          A_n = zeros(numClusters, 1); % Initialize A_n vector
          B_n = zeros(numClusters, 1); % Initialize B_n vector
          A_f = zeros(numClusters, 1); % Initialize A_f vector
          B_f = zeros(numClusters, 1); % Initialize B_f vector
          A_c_n = zeros(numClusters, 1); % Initialize A_c_n vector
          B_c_n = zeros(numClusters, 1); % Initialize B_c_n vector
          R_n=zeros(numClusters, 1); % Initialize R_n vector
          R_f=zeros(numClusters, 1); % Initialize R_f vector    
          R_c_n=zeros(numClusters, 1); % Initialize R_c_n vector
          noise = para.noise;  % Noise power
          WSR=0;

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

                    A_n(c)= abs(H_n{c}).^2 * alpha_n(c);

                    B_n(c) =inter_cluster_interference_near + ...
                            abs(H_n_c{c}).^2 * eta_k + noise;
        
                    A_f(c) = abs( H_f{c}).^2 * alpha_f(c);

                    B_f(c) = inter_cluster_interference_far + ...
                            abs(H_f{c}).^2  * alpha_n(c) + ...
                            abs(H_f_c{c}).^2  * eta_k + noise;

                    A_c_n(c) = abs(H_n_c{c}).^2 * eta_k;

                    B_c_n(c) = inter_cluster_interference_near_b + noise;


                  


                    % compute rates
                    R_f(c) = log2(1 + A_f(c) / B_f(c)); % Far user rate
                    R_n(c) = log2(1 + A_n(c) / B_n(c)); % Near user rate
                    R_c_n(c) = log2(1 + A_c_n(c) / B_c_n(c));

                    WSR =WSR+ para.weights_n*R_n(c) + para.weights_f*R_f(c)+R_c_n(c); % Weighted sum rate

           end
            % Compute WSR    
end

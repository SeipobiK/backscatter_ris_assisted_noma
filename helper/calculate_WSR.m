function [WSR,R_n,R_f,R_c_n,A_f,A_n] = calculate_WSR(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, alpha_n, alpha_f, Theta)
         % power= abs(w_k(:,1))^2 + abs(w_k(:,2))^2;
         % disp(['Total power in Compute WSR: ', num2str(power)]);         
          numClusters = para.K; % Number of 
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

           for c = 1:numClusters
                    inter_cluster_interference_near = 0;
                    inter_cluster_interference_far = 0;
                    for j = 1:numClusters
                        if j ~= c
                            % Near user inter cluster interference  
                            inter_cluster_interference_near = inter_cluster_interference_near + ...
                                abs(g_1_all{c}'*Theta*G_all*w_k(:, j)).^2;

                            inter_cluster_interference_far= inter_cluster_interference_far + ...
                               abs(g_2_all{c}'*Theta*G_all*w_k(:, j)).^2;
                            % 
                            % inter_cluster_interference_near_b= inter_cluster_interference_near_b + ...
                            %     abs(g_b_all{c}'*Theta*G_all*f1_all{c}*w_k(:, j)).^2;
                            
                        end

                    end

                    A_n(c) = abs(g_1_all{c}'*Theta*G_all*w_k(:, c)).^2 * alpha_n(c); % Near user


                    B_n(c) =inter_cluster_interference_near ++ noise;

                        
                    A_f(c) = abs(g_2_all{c}'*Theta*G_all*w_k(:, c)).^2  * alpha_f(c);


                    B_f(c) = inter_cluster_interference_far + ...
                            abs(g_2_all{c}'*Theta*G_all*w_k(:, c)).^2  * alpha_n(c) + noise;


                    % A_c_n(c) = abs(g_b_all{c}'*Theta*G_all*f1_all{c}*w_k(:, c)).^2* para.eta;
                    % 
                    % B_c_n(c) = inter_cluster_interference_near_b + noise;

                    % compute rates
                    R_f(c) = log2(1 + A_f(c) / B_f(c)); % Far user rate
                    R_n(c) = log2(1 + A_n(c) / B_n(c)); % Near user rate
                    WSR =WSR+ para.weights_n.*R_n(c) + para.weights_f.*R_f(c); % Weighted sum rate

           end
            % Compute WSR    
end

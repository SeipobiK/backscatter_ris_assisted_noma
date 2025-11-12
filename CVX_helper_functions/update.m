function [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, status] = update(para, H_n, H_f, H_n_c, H_f_c, A_n_prev, ...
    B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, alpha_f, alpha_n)
    % Extract configuration parameters
    numClusters = para.K;
    M = para.M;
    R_n_min = para.R_min_n;
    R_f_min = para.R_min_f;
    R_c_min = para.R_c_min;
    eta_k = para.eta;
    P_max = para.P_max;
    noise = para.noise;
    
    % Precompute constants and channel products
    log2_e = log2(exp(1));
    
    % Precompute channel products to avoid repeated calculations
    H_n_H_n = cell(numClusters, 1);
    H_f_H_f = cell(numClusters, 1);
    H_n_c_H_n_c = cell(numClusters, 1);
    H_f_c_H_f_c = cell(numClusters, 1);
    
    for c = 1:numClusters
        H_n_H_n{c} = H_n{c}' * H_n{c};
        H_f_H_f{c} = H_f{c}' * H_f{c};
        H_n_c_H_n_c{c} = H_n_c{c}' * H_n_c{c};
        H_f_c_H_f_c{c} = H_f_c{c}' * H_f_c{c};
    end
    
    % Precompute Taylor approximation constants
    taylor_const_far = zeros(numClusters, 3);
    taylor_const_near = zeros(numClusters, 3);
    taylor_const_back = zeros(numClusters, 3);
    
    for c = 1:numClusters
        % Far user constants
        denom_far = 1 + A_f_prev(c) * B_f_prev(c);
        taylor_const_far(c, 1) = log2(1 + 1/(A_f_prev(c) * B_f_prev(c)));
        taylor_const_far(c, 2) = log2_e / (A_f_prev(c) * denom_far);
        taylor_const_far(c, 3) = log2_e / (B_f_prev(c) * denom_far);
        
        % Near user constants
        denom_near = 1 + A_n_prev(c) * B_n_prev(c);
        taylor_const_near(c, 1) = log2(1 + 1/(A_n_prev(c) * B_n_prev(c)));
        taylor_const_near(c, 2) = log2_e / (A_n_prev(c) * denom_near);
        taylor_const_near(c, 3) = log2_e / (B_n_prev(c) * denom_near);
        
        % Backscatter constants
        denom_back = 1 + A_c_prev_n(c) * B_c_prev_n(c);
        taylor_const_back(c, 1) = log2(1 + 1/(A_c_prev_n(c) * B_c_prev_n(c)));
        taylor_const_back(c, 2) = log2_e / (A_c_prev_n(c) * denom_back);
        taylor_const_back(c, 3) = log2_e / (B_c_prev_n(c) * denom_back);
    end

    cvx_begin quiet
        cvx_solver mosek
        
        % % Optimized solver settings
        % cvx_solver_settings('MSK_IPAR_NUM_THREADS', 4);
        % cvx_solver_settings('MSK_IPAR_INTPNT_BASIS', 0);
        % cvx_solver_settings('MSK_DPAR_INTPNT_TOL_PFEAS', 1e-8);
        % cvx_solver_settings('MSK_DPAR_INTPNT_TOL_DFEAS', 1e-8);
        % cvx_solver_settings('MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-6);
        
        % Variable declarations
        variable W(M, M, numClusters) Hermitian semidefinite
        variable A_n(numClusters) nonnegative
        variable B_n(numClusters) nonnegative
        variable A_f(numClusters) nonnegative
        variable B_f(numClusters) nonnegative
        variable A_c_n(numClusters) nonnegative
        variable B_c_n(numClusters) nonnegative
        variable R_n(numClusters) nonnegative
        variable R_f(numClusters) nonnegative
        variable R_c_n(numClusters) nonnegative

        % Objective function
        maximize(sum(para.weights_n * R_n + para.weights_f * R_f + para.weights_c * R_c_n))

        subject to
            % Power constraint
            sum_power = 0;
            for k = 1:numClusters
                sum_power = sum_power + trace(W(:, :, k));
            end
            sum_power <= P_max;

            for c = 1:numClusters
                % Rate constraints using precomputed Taylor approximations
                R_f(c) <= taylor_const_far(c, 1) - ...
                         taylor_const_far(c, 2) * (A_f(c) - A_f_prev(c)) - ...
                         taylor_const_far(c, 3) * (B_f(c) - B_f_prev(c));
                    
                R_n(c) <= taylor_const_near(c, 1) - ...
                         taylor_const_near(c, 2) * (A_n(c) - A_n_prev(c)) - ...
                         taylor_const_near(c, 3) * (B_n(c) - B_n_prev(c));

                R_c_n(c) <= taylor_const_back(c, 1) - ...
                           taylor_const_back(c, 2) * (A_c_n(c) - A_c_prev_n(c)) - ...
                           taylor_const_back(c, 3) * (B_c_n(c) - B_c_prev_n(c));

                % Minimum rate constraints
                R_f(c) >= R_f_min;
                R_n(c) >= R_n_min;
                R_c_n(c) >= R_c_min;

                % Interference calculations - CORRECTED VERSION
                inter_cluster_interference_near = 0;
                inter_cluster_interference_far = 0;
                inter_cluster_interference_near_b = 0;
                
                for j = 1:numClusters
                    if j ~= c
                        % Near user inter cluster interference  
                        inter_cluster_interference_near = inter_cluster_interference_near + ...
                            trace(W(:, :, j) * H_n_H_n{c});
                            
                        inter_cluster_interference_far = inter_cluster_interference_far + ...
                            trace(W(:, :, j) * H_f_H_f{c});
                            
                        inter_cluster_interference_near_b = inter_cluster_interference_near_b + ...
                            trace(W(:, :, j) * H_n_c_H_n_c{c});
                    end
                end

                % Main constraints
                inv_pos(A_n(c)) <= trace(W(:, :, c) * H_n_H_n{c}) * alpha_n(c);

                B_n(c) >= inter_cluster_interference_near + ...
                        trace(W(:, :, c) * H_n_c_H_n_c{c}) * eta_k + noise;
        
                inv_pos(A_f(c)) <= trace(W(:, :, c) * H_f_H_f{c}) * alpha_f(c);

                B_f(c) >= inter_cluster_interference_far + ...
                        trace(W(:, :, c) * H_f_H_f{c}) * alpha_n(c) + ...
                        trace(W(:, :, c) * H_f_c_H_f_c{c}) * eta_k + noise;

                inv_pos(A_c_n(c)) <= trace(W(:, :, c) * H_n_c_H_n_c{c}) * eta_k;

                B_c_n(c) >= inter_cluster_interference_near_b + noise;
            end
    cvx_end
    
    % Return results
    obj_prev = cvx_optval;
    A_n_opt = A_n;
    B_n_opt = B_n;
    A_f_opt = A_f;
    B_f_opt = B_f;
    A_c_n_opt = A_c_n;
    B_c_n_opt = B_c_n;
    W_opt = W;
    status = cvx_status;
    
    % Optional: Display results (comment out for maximum speed)
    if nargout == 0
        disp('R_n:'); disp(R_n);
        disp('R_f:'); disp(R_f);
        disp('R_c_n:'); disp(R_c_n);
    end
end













% 
% function [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = update(para,H_n, H_f, H_n_c, H_f_c,A_n_prev, ...
%     B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,alpha_f,alpha_n)
%     % Extract configuration parameters
% 
%     numClusters = para.K; % Number of clusters
%     M = para.M; % Number of BS antennas
%     R_n_min = para.R_min_n; % Minimum rate for near user
%     R_f_min = para.R_min_f; % Minimum rate for far user
%     R_c_min = para.R_c_min; % Minimum rate for backscatter user
%     eta_k = para.eta; % Backscatter coefficient
%     P_max = para.P_max; % Maximum transmit power
%     noise = para.noise;  % Noise scales with power
% 
% 
% 
%     cvx_begin quiet    sdp
%         cvx_solver mosek
% 
%         variable W(M,M,numClusters) Hermitian semidefinite
%         variable A_n(numClusters) nonnegative % Slack variable for near users
%         variable B_n(numClusters) nonnegative % Slack variable for near users
%         variable A_f(numClusters) nonnegative % Slack variable for far users
%         variable B_f(numClusters) nonnegative % Slack variable for far users
%         variable A_c_f(numClusters) nonnegative % Slack variable for backscatter devices at far user
%         variable B_c_f(numClusters) nonnegative % Slack variable for backscatter devices at far user
%         variable A_c_n(numClusters) nonnegative % Slack variable for backscatter devices at near user
%         variable B_c_n(numClusters) nonnegative % Slack variable for backscatter devices at near user
%         variable R_n(numClusters)  nonnegative% Slack variable for backscatter devices at near user
%         variable R_f(numClusters)  nonnegative% Slack variable for backscatter devices at near user
%         variable R_c_n(numClusters)  nonnegative% Slack variable for backscatter devices at near user
% 
%         % Objective function: Maximize weighted sum rate
%        maximize(sum(para.weights_n*R_n + para.weights_f*R_f + para.weights_c*R_c_n)) 
% 
%          subject to
% 
%             sum_power = 0;
%             for k = 1:numClusters
%                 sum_power = sum_power + real(trace(W(:,:,k)));
%             end
%             sum_power <= P_max;
% 
% 
%             for c = 1:numClusters
% 
%                     R_f(c) <= log2(1 + 1 ./ (A_f_prev(c) * B_f_prev(c))) -  ...
%                     (log2(exp(1)) * 1 ./ (A_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (A_f(c) - A_f_prev(c)) - ...
%                     (log2(exp(1)) * 1 ./ (B_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (B_f(c) - B_f_prev(c));
% 
% 
%                     R_n(c) <= log2(1 + 1 ./ (A_n_prev(c) * B_n_prev(c))) -  ...
%                     (log2(exp(1)) * 1 ./ (A_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (A_n(c) - A_n_prev(c)) - ...
%                     (log2(exp(1)) * 1 ./ (B_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (B_n(c) - B_n_prev(c));
% 
% 
% 
%                     R_c_n(c) <= log2(1 + 1 ./ (A_c_prev_n(c) * B_c_prev_n(c))) - ...
%                     (log2(exp(1)) *1 ./  (A_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (A_c_n(c) - A_c_prev_n(c)) - ...
%                     (log2(exp(1)) * 1 ./  (B_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (B_c_n(c) - B_c_prev_n(c));
% 
%                     R_f(c) >= R_f_min;
%                     R_n(c) >= R_n_min;
%                     R_c_n(c) >= R_c_min;
% 
%                     inter_cluster_interference_near = 0;
%                     inter_cluster_interference_far = 0;
%                     inter_cluster_interference_near_b=0;
%                     for j = 1:numClusters
%                         if j ~= c
%                             % Near user inter cluster interference  
%                             inter_cluster_interference_near = inter_cluster_interference_near + ...
%                                 real(trace(W(:,:,j) * H_n{c}' * H_n{c}));
% 
%                             inter_cluster_interference_far= inter_cluster_interference_far + ...
%                                 real(trace(W(:,:,j) * H_f{c}' * H_f{c})); 
% 
%                             inter_cluster_interference_near_b= inter_cluster_interference_near_b + ...
%                                 real(trace(W(:,:,j)* H_n_c{c}' * H_n_c{c}));          
%                         end
% 
%                     end
% 
%                     inv_pos(A_n(c)) <= real(trace(W(:,:,c) * H_n{c}' * H_n{c})) * alpha_n(c);
% 
%                     B_n(c) >=inter_cluster_interference_near + ...
%                             real(trace(W(:,:,c) * H_n_c{c}' * H_n_c{c})) * eta_k + noise;
% 
%                     inv_pos(A_f(c)) <= real(trace(W(:,:,c) * H_f{c}' * H_f{c})) * alpha_f(c);
% 
%                     B_f(c) >= inter_cluster_interference_far + ...
%                             real(trace(W(:,:,c) * H_f{c}' * H_f{c}))  * alpha_n(c) + ...
%                             real(trace(W(:,:,c) * H_f_c{c}' * H_f_c{c}))   * eta_k + noise;
% 
%                     inv_pos(A_c_n(c)) <= real(trace(W(:,:,c) * H_n_c{c}' * H_n_c{c})) * eta_k;
% 
%                     B_c_n(c) >= inter_cluster_interference_near_b + noise;
%             end
%     cvx_end
%     obj_prev = cvx_optval;
%     A_n_opt = A_n;
%     B_n_opt = B_n;
%     A_f_opt = A_f;
%     B_f_opt = B_f;
%     A_c_n_opt = A_c_n;
%     B_c_n_opt = B_c_n;
%     W_opt = W;
%     status = cvx_status;
%     disp(R_n);
%     disp(R_f);
%     disp(R_c_n);
% 
% 
% end
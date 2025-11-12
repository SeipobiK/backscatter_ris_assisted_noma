function [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, status] = feasible_2(para, H_n, H_f, H_n_c, ...
        H_f_c, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n,alpha_f,alpha_n)
    % Extract configuration parameters
    numClusters = para.K;
    M = para.M;
    R_n_min = para.R_min_n;
    R_f_min = para.R_min_f;
    R_c_min = para.R_c_min;
    eta_k = para.eta;
    P_max = para.P_max;

    cvx_begin quiet sdp
        cvx_solver mosek

        
        % Variable declarations
        variable W(M, M, numClusters) Hermitian semidefinite
        variable A_n(numClusters) nonnegative
        variable B_n(numClusters) nonnegative
        variable A_f(numClusters) nonnegative
        variable B_f(numClusters) nonnegative
        variable A_c_f(numClusters) nonnegative
        variable B_c_f(numClusters) nonnegative
        variable A_c_n(numClusters) nonnegative
        variable B_c_n(numClusters) nonnegative
        variable delta_g nonnegative
        variable R_n(numClusters) nonnegative
        variable R_f(numClusters) nonnegative
        variable R_c_n(numClusters) nonnegative
       
        
        minimize(delta_g)
        
        subject to
            % Power constraint
            sum_power = 0;
            for k = 1:numClusters
                sum_power = sum_power + real(trace(W(:, :, k)));
            end
            sum_power - P_max <= delta_g;
            
            for c = 1:numClusters
                
                    taylor_approx_far(c) = log2(1 + 1 ./ (A_f_prev(c) * B_f_prev(c))) -  ...
                    (log2(exp(1)) * 1 ./ (A_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (A_f(c) - A_f_prev(c)) - ...
                    (log2(exp(1)) * 1 ./ (B_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (B_f(c) - B_f_prev(c));

                    
                    taylor_approx_n(c) = log2(1 + 1 ./ (A_n_prev(c) * B_n_prev(c))) -  ...
                    (log2(exp(1)) * 1 ./ (A_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (A_n(c) - A_n_prev(c)) - ...
                    (log2(exp(1)) * 1 ./ (B_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (B_n(c) - B_n_prev(c));


                
                % Rate constraints
                R_f(c) - taylor_approx_far(c) <= delta_g;
                R_n(c) - taylor_approx_n(c) <= delta_g;
                % R_c_n(c) - taylor_approx_backscatter_n(c) <= delta_g;
                
                R_f_min - R_f(c) <= delta_g;
                R_n_min - R_n(c) <= delta_g;
                % R_c_min - R_c_n(c) <= delta_g;
                
                % Inter-cluster interference calculations
                inter_cluster_interference_near = 0;
                inter_cluster_interference_far = 0;
                test_int=0;
                testa_b=0;
                
                for j = 1:numClusters
                    if j ~= c
                        inter_cluster_interference_near = inter_cluster_interference_near + ...
                            real(trace(W(:, :, j) * H_n{c}' * H_n{c}));
                        test_int = test_int + trace(H_n{c}' * H_n{c});
                            
                        
                        inter_cluster_interference_far = inter_cluster_interference_far + ...
                            real(trace(W(:, :, j) * H_f{c}' * H_f{c}));
                        
                        % inter_cluster_interference_near_b = inter_cluster_interference_near_b + ...
                        %     real(trace(W(:, :, j) * H_n_c{c}' * H__c{c}));

                            testa_b = testa_b + trace(H_n_c{c}' * H_n_c{c});
                    end
                end
                
                % Create separate constraints for problematic ones
                % Constraint 1

            disp(['Near Interference :', num2str(test_int)]);
            disp(['far BST Interference :', num2str(testa_b)]);

            disp(['Near channel :',num2str(trace(H_n{c}' * H_n{c}))]);
            disp(['Far channel :',num2str(trace(H_f{c}' * H_f{c}))]);
            disp(['BST Near channel :',num2str(trace(H_n_c{c}' * H_n_c{c}))]);
            disp(['BST Far channel :',num2str(trace(H_f_c{c}' * H_f_c{c}))]);
            disp(['Noise :',num2str(para.noise)]);
            disp(['alpha_f :',num2str(alpha_f(c))]);
            disp(['alpha_f :',num2str(alpha_n(c))]);
            
               inv_pos(A_n(c)) - real(trace(W(:, :, c) * H_n{c}' * H_n{c})) * alpha_n(c) <= delta_g;
                
                % Constraint 2
                inter_cluster_interference_near  + para.noise - B_n(c) <= delta_g;
                
                % Constraint 3
                inv_pos(A_f(c)) - real(trace(W(:, :, c) * H_f{c}' * H_f{c})) * alpha_f(c) <= delta_g;
                
                % Constraint 4
                inter_cluster_interference_far + ...
                    real(trace(W(:, :, c) * H_f{c}' * H_f{c})) * alpha_n(c) - B_f(c) <= delta_g;
            end

            
            
            % PSD constraints
            for k = 1:numClusters
                W(:, :, k) + delta_g * eye(M) == hermitian_semidefinite(M);
            end
    cvx_end
 
    % Return optimized results
    obj_prev = cvx_optval;
    A_n_opt = A_n;
    B_n_opt = B_n;
    A_f_opt = A_f;
    B_f_opt = B_f;
    A_c_n_opt = 0;
    B_c_n_opt = 0;
    W_opt = W;
    status = cvx_status;
end
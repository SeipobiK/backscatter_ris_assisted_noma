function [W_opt, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, obj_history,converged,cvx_status] = ...
    find_feasible_ndris_active(para,Theta,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n,alpha_f,alpha_n, max_iter, mc,J_t,J_r)
    K=para.K;
    N=para.N;
    numClusters=para.K;


    % Initialize
    obj_history = zeros(max_iter, 1);
    converged = false;

    H_n = cell(1, K); H_f = cell(1, K);
    H_nc = cell(1, K); H_fc = cell(1, K);

    for c=1:K
            % disp('========================Active Feasible BFFFF======================================');

             H_n{c}  = g_1_all{c}'*J_r*Theta*J_t*G_all;
             H_f{c}  = g_2_all{c}'*J_r*Theta*J_t*G_all;

             H_nc{c}  = g_b_all{c}'*J_r*Theta*J_t*f1_all{c}*G_all;
             H_fc{c}  = g_b_all{c}'*J_r*Theta*J_t*f2_all{c}*G_all;  

           % disp('========================= END ======================================');



    end      

    for m = 1:max_iter
        % Update beamforming and Taylor parameters
        [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
            feasible_2(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n,alpha_f,alpha_n);

        if ~strcmp(cvx_status, 'Solved')
            warning('Update failed at MC %d iteration %d', mc, m);
            break;
        end

        % Update variables
        A_n_prev = A_n; B_n_prev = B_n;
        A_f_prev = A_f; B_f_prev = B_f;
        A_c_prev_n = A_cn; B_c_prev_n = B_cn;
        % W_init = W_opt;
        obj_history(m) = obj_curr;
        % disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_prev')]);
        % disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_prev')]);
        % disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_prev')]);
        % disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_prev')]);
        % disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_prev_n')]);
        % disp(['Iteration: ', num2str(m), ' pac far: ', num2str(alpha_f')]);
        % disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_prev_n')]);
        
        % % Display progress
        % disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.2e', obj_curr)]);
  
        % Check convergence
        if m > 1 && abs(obj_history(m)) < 1e-8
            converged = true;
            obj_history = obj_history(1:m);  % Trim unused entries
             % disp(['    Objective at Convergence: ', sprintf('%.2e', abs(obj_history(m)))]);
             converged = true;
             break;
        end

        % if m==24 && ~converged
        %     gjlkg
        %     break;
        % 
        % end
    end
end
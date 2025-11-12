function [W_opt, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, obj_history,obj_history_mc, converged,cvx_status] = ...
    active_BF_ndrs(para, H_n, ...
                     H_f,H_n_c,H_f_c, A_n_prev, B_n_prev, ...
                     A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n,alpha_f,alpha_n, max_iter)

    % Initialize
    obj_history = zeros(max_iter, 1);
    obj_history_mc = zeros(max_iter, 1);
    converged = false;   

    for m = 1:max_iter
        % Update beamforming and Taylor parameters
        [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
            update(para, H_n, H_f, H_n_c, H_f_c, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n,alpha_f,alpha_n);

        if ~strcmp(cvx_status, 'Solved')
            % disp('failedddd');
            % disp(cvx_status);
            break;
        end

        % Update variables
        A_n_prev = A_n; B_n_prev = B_n;
        A_f_prev = A_f; B_f_prev = B_f;
        A_c_prev_n = A_cn; B_c_prev_n = B_cn;
        % W_init = W_opt;
        obj_history(m) = obj_curr;
        obj_history_mc(m) = obj_curr;  % Store WSR for this iteration
        % 
        % disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_prev')]);
        % disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_prev')]);
        % disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_prev')]);
        % disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_prev')]);
        % disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_prev_n')]);
        % disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_prev_n')]);

        % Display progress
        % disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_curr)]);
        % disp(['    WSR: ', num2str(WSR)]);
        if m > 1
            % disp(['    Change: ', sprintf('%.10f', abs(obj_history(m) - obj_history(m-1)))]);
            % disp(['    Change(Calculated): ', sprintf('%.10f', obj_history_mc(m) - obj_history_mc(m-1))]);
        end
        for k = 1:size(W_opt, 3)
            % disp(['    Rank(W_', num2str(k), '): ', num2str(rank(W_opt(:,:,k)))]);
                    current_eig = eig(W_opt(:,:,k));
                    sorted_eig = sort(current_eig, 'descend');  
                    % disp(cvx_status)
                    % disp(['    Eigenvalues of W_', num2str(k), ': ', num2str(sorted_eig')]);
        end

        % Check convergenc
        if m > 1 && abs(obj_history(m) - obj_history(m-1)) < 1e-5
            converged = true;
            obj_history = obj_history(1:m);
            break;
        end
    end
end
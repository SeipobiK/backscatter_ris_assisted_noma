% function [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_history,obj_history_mc,cvx_status] =passive_bf_ndris_opt(para,w_k,G_all, g_1_all,...
%     g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n,alpha_f,alpha_n , max_iter,mc,J_t,J_r)
% 
% 
% 
%     % Initialize
%     obj_history = zeros(max_iter, 1);
%     obj_history_mc = zeros(max_iter, para.MC_MAX);
%     converged = false;
%     N=para.N;
%     max_iter=20;
% 
%      % Solve the relaxed problem
%     disp(['value of A_F :', num2str(B_f_prev')]);
% 
%     [V_opt_init,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_curr,status] = passive_relaxed_ndris(para,w_k,G_all, g_1_all,...
%      g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,alpha_f,alpha_n,J_t,J_r);
% 
%             A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
%             A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
%             A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
% 
% 
%     % % Extract V_0 and relaxation parameter_0
%     step_size = zeros(max_iter,1);
%     relax_parameter=zeros(max_iter,1);
%     max_eigenvalue_V_opt=zeros(max_iter,1);
%     max_eigVector_V_opt = zeros(N,max_iter);
%     [V_max, lambda_max] = max_eigVect(V_opt_init);
%     max_eigenvalue_V_opt(1)=lambda_max;
%     max_eigVector_V_opt(:,1)=V_max;
%     U_opt = zeros(N,N,max_iter);
%     U_opt(:,:,1) = V_opt_init;
%     obj_history(1) = obj_curr;
%     initial_ratio = max_eigenvalue_V_opt(1)/trace(U_opt(:,:,1));
%     step_size(1) = (1-initial_ratio)*0.7;
% 
% 
% 
%     relax_parameter(1)=min(1, max_eigenvalue_V_opt(1)/trace(U_opt(:,:,1)) + step_size(1));
% 
%     disp(relax_parameter(1));
% 
%     for m = 2:max_iter
% 
%         if ~strcmp(status, 'Solved')
%             V_opt=V_opt_init;
%             warning('Update failed at MC %d iteration %d', mc, m);
%             break;
%         end
% 
% 
% 
%          [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = passive_ndris(para,w_k,G_all, g_1_all,...
%             g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,max_eigVector_V_opt(:,m-1),relax_parameter(m-1),alpha_f,alpha_n,J_t,J_r);
% 
% 
%         if strcmp(cvx_status, 'Solved')
%             % Update variables
%             A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
%             A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
%             A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
% 
%             obj_history(m) = obj_prev;
%             % WSR = calculate_WSR(para, w_k, G_all, g_1_all, g_2_all, ...
%             %                         g_b_all, f1_all, f2_all, para.alpha_k_n, ...
%             %                         para.alpha_k_f, Theta);
% 
%             obj_history_mc(m,mc)=obj_prev;
% 
%             % Extract Max eigen Value and Max eigen Vector 
%             [V_max_, eig_max] = max_eigVect(V_opt);
%             U_opt(:,:,m)=V_opt;
%             max_eigenvalue_V_opt(m)=eig_max;
%             max_eigVector_V_opt(:,m)=V_max_;
%             current_eig_ = eig(V_opt);
%             sorted_eig_ = sort(current_eig_, 'descend');
%             step_size(m)=step_size(1);
%             % % 4. Calculate current ratio
%             current_ratio = eig_max/trace(V_opt);
% 
%             % % Display progress
%              disp(['Iteration: ', num2str(m-1), ' | Objective: ', sprintf('%.10f', obj_history(m-1))]);
%              if m > 1
%                  disp(['Change: ', sprintf('%.10f', (obj_history(m) - obj_history(m-1)))]);
%                  disp(current_ratio+step_size(1));
%              end
%              disp(['    Rank(V_', num2str(m), '): ', num2str(rank(V_opt))]);
%         else
% 
%             U_opt(:,:,m) = U_opt(:,:,m-1);
%             max_eigenvalue_V_opt(m) = max_eigenvalue_V_opt(m-1);
%             obj_history(m) = obj_history(m-1);
%             step_size(m)=step_size(m-1)/2;
%             disp('failed');
%             % disp(step_size(m));
%             % obj_history = obj_history(1:m);
%             % relax_parameter = relax_parameter(1:m);
%             if step_size(m)<1e-3
%                obj_history=NaN;
%                break;
%             end
% 
%         end
%         % disp(max_eigenvalue_V_opt(1)/trace(U_opt(:,:,m)));
%         current_ratio = max_eigenvalue_V_opt(m)/trace(U_opt(:,:,m));
% 
%             if  obj_history(m) - obj_history(m-1)<0
%                    relax_parameter(m) = min(1, current_ratio + 0.71*(1-current_ratio)); % 10% step   
% 
%                    disp(['Relax  parametr :',num2str(relax_parameter(m))]);                   
%             else 
% 
%                   relax_parameter(m) = min(1, current_ratio + 0.75*(1-current_ratio)); % 10% step
%                   disp(['Relax  parametr :',num2str(relax_parameter(m))]);  
%             end
% 
%         % Check convergence
%         if m > 1 && abs(obj_history(m) - obj_history(m-1)) < 1e-4 && strcmp(cvx_status, 'Solved') && abs(1-relax_parameter(m)) <= 1e-4
%              converged = true;
%              % obj_history = obj_history(1:m); 
%                          break;    
%         end
%     end
% 
% end


function [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_history,obj_history_mc,converged] =passive_bf_ndris_opt(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n,alpha_f,alpha_n , max_iter,mc,J_t,J_r)



    % Initialize
    obj_history = zeros(max_iter, 1);
    obj_history_mc = zeros(max_iter, para.MC_MAX);
    converged = false;
    N=para.N;
    max_iter=25;

     % Solve the relaxed problem
    % disp(['value of A_F :', num2str(B_f_prev')]);

    [V_opt_init,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_curr,status] = passive_relaxed_ndris(para,w_k,G_all, g_1_all,...
     g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,alpha_f,alpha_n,J_t,J_r);

            A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
            A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
            A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 


    % % Extract V_0 and relaxation parameter_0
    step_size = zeros(max_iter,1);
    relax_parameter=zeros(max_iter,1);
    max_eigenvalue_V_opt=zeros(max_iter,1);
    max_eigVector_V_opt = zeros(N,max_iter);
    [V_max, lambda_max] = max_eigVect(V_opt_init);
    max_eigenvalue_V_opt(1)=lambda_max;
    max_eigVector_V_opt(:,1)=V_max;
    U_opt = zeros(N,N,max_iter);
    U_opt(:,:,1) = V_opt_init;
    obj_history(1) = obj_curr;
    initial_ratio = max_eigenvalue_V_opt(1)/trace(U_opt(:,:,1));
    step_size(1) = (1-initial_ratio)*0.7;

    relax_parameter(1)=min(1, max_eigenvalue_V_opt(1)/trace(U_opt(:,:,1)) + step_size(1));

    for m = 2:max_iter

        if ~strcmp(status, 'Solved')
            V_opt=V_opt_init;
            warning('Update failed at MC %d iteration %d', mc, m);
            break;
        end

         [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = passive_ndris(para,w_k,G_all, g_1_all,...
            g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,max_eigVector_V_opt(:,m-1),relax_parameter(m-1),alpha_f,alpha_n,J_t,J_r);


        if strcmp(cvx_status, 'Solved')
            % Update variables
            A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
            A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
            A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 

            obj_history(m) = obj_prev;
            obj_history_mc(m,mc)=obj_prev;

            % Extract Max eigen Value and Max eigen Vector 
            [V_max_, eig_max] = max_eigVect(V_opt);
            U_opt(:,:,m)=V_opt;
            max_eigenvalue_V_opt(m)=eig_max;
            max_eigVector_V_opt(:,m)=V_max_;
            current_eig_ = eig(V_opt);
            sorted_eig_ = sort(current_eig_, 'descend');
            step_size(m)=step_size(1);
            % % 4. Calculate current ratio
            current_ratio = eig_max/trace(V_opt);

            % % Display progress
            %  disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_history(m))]);
             if m > 1
                %  disp(['Change: ', sprintf('%.10f', (obj_history(m) - obj_history(m-1)))]);
                 % disp(current_ratio+step_size(m));
                 % disp(step_size(m));
             end
             % disp(['    Rank(V_', num2str(m), '): ', num2str(rank(V_opt))]);
             % disp(sorted_eig_);

             step_size(m)=0.5*(1-current_ratio);

        else
            U_opt(:,:,m) = U_opt(:,:,m-1);
            max_eigenvalue_V_opt(m) = max_eigenvalue_V_opt(m-1);
            obj_history(m) = obj_history(m-1);
            step_size(m)=step_size(m)/2;
            % disp('failed due to step size :');
            if step_size(m)<1e-6
               obj_history=NaN;
               break;
            end

        end
        % disp(max_eigenvalue_V_opt(1)/trace(U_opt(:,:,m)));
        % current_ratio = max_eigenvalue_V_opt(m)/trace(U_opt(:,:,m));


        % relax_parameter(m) = min(1, current_ratio + step_size(m)); % 10% step
        current_ratio = max_eigenvalue_V_opt(m)/trace(U_opt(:,:,m));



        relax_parameter(m) = min(1, current_ratio + 0.65*(1-current_ratio)); % 10% step
        % disp(['Relax  parametr :',num2str(relax_parameter(m))]);  
        % disp(['Step size :',num2str(step_size(m))]);
      
        % Check convergence
        if m > 1 && abs(obj_history(m) - obj_history(m-1)) < 1e-4 && strcmp(cvx_status, 'Solved') && abs(1-relax_parameter(m)) <= 1e-5
             converged = true;
             obj_history = obj_history(1:m); 
            break;
        end
    end
end
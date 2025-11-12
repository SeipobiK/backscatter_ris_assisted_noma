function [g_1_all, g_2_all, g_b_all, f1_all, f2_all, decoding_order] = ensure_decoding_order(para, Theta, w_k, G_all_matrix, g_local, f_local, J_t, J_r)
    K = para.K;
    noise=para.noise;
    eta_k=para.eta;
    g_1_all = cell(1, K);
    g_2_all = cell(1, K);
    g_b_all = cell(1, K);
    f1_all = cell(1, K);
    f2_all = cell(1, K);
    decoding_order = cell(1, K);

    scal=para.scall;
    
    for c = 1:K
        % Calculate equivalent combined channel gains with proper normalization
        H_n_combined = g_local(:, c, 1)' * J_r * Theta * J_t * G_all_matrix * w_k(:, c);
        H_f_combined = g_local(:, c, 2)' * J_r * Theta * J_t * G_all_matrix * w_k(:, c);
        H_n_combined_c = g_local(:, c, 2)' * J_r * Theta * J_t * G_all_matrix* f_local(c, 1)  * w_k(:, c);
        H_f_combined_c = g_local(:, c, 2)' * J_r * Theta * J_t * G_all_matrix * f_local(c, 2) * w_k(:, c);


        
        % Calculate inter-cluster interference for proper gain comparison
        inter_cluster_near = 0;
        inter_cluster_far = 0;
        
        for j = 1:K
            if j ~= c
                inter_cluster_near = inter_cluster_near + ...
                    abs(g_local(:, c, 1)' * J_r * Theta * J_t * G_all_matrix * w_k(:, j))^2;
                inter_cluster_far = inter_cluster_far + ...
                    abs(g_local(:, c, 2)' * J_r * Theta * J_t * G_all_matrix * w_k(:, j))^2;
            end
        end
        
        % Equivalent combined channel gains (following paper's Eq. 11)
        gamma_n = abs(H_n_combined)^2/(abs(H_n_combined_c)^2*eta_k + inter_cluster_near  + noise);
        gamma_f = abs(H_f_combined)^2/(abs(H_f_combined_c)^2 *eta_k + inter_cluster_far + noise);
        
        % Ensure near user has better channel gain (decodes later in SIC)
        if gamma_n >= gamma_f
            % Near user is stronger (good) - decode later
            g_1_all{c} = g_local(:, c, 1)*scal;  % Near user
            g_2_all{c} = g_local(:, c, 2)*scal;  % Far user  
            f1_all{c} = f_local(c, 1);
            f2_all{c} = f_local(c, 2);
            decoding_order{c} = 'near_strong';
            % disp(['Cluster ', num2str(c), ': Near user stronger - Gamma_n=', ...
            %       num2str(gamma_n), ', Gamma_f=', num2str(gamma_f)]);
        else
            % Far user is stronger (bad) - swap assignment
            g_1_all{c} = g_local(:, c, 2)*scal;  % Now: g_1 = far user (weaker)
            g_2_all{c} = g_local(:, c, 1)*scal;  % Now: g_2 = near user (stronger)
            f1_all{c} = f_local(c, 2);
            f2_all{c} = f_local(c, 1);
            decoding_order{c} = 'far_strong_swapped';

            % disp(['Cluster ', num2str(c), ': SWAPPED - Gamma_n=', ...
            %       num2str(gamma_n), ', Gamma_f=', num2str(gamma_f)]);
        end
        
        g_b_all{c} = g_local(:, c, 3)*scal;  % Backscatter user remains same
    end
end
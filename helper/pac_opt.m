function [alpha_n, alpha_f] = pac_opt(para,w_k,G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, Theta)
    numClusters = para.K; % Number of clusters
    alpha_f = zeros(numClusters,1);
    alpha_n = zeros(numClusters,1);

    noise = para.noise;  % Noise power

    for c = 1:numClusters
        % Inter-cluster interference
        inter_cluster_near = 0;
        inter_cluster_far = 0;

        for j = 1:numClusters
            if j ~= c
                inter_cluster_near = inter_cluster_near + abs(g_1_all{c}'*Theta*G_all*w_k(:, j))^2;
                inter_cluster_far  = inter_cluster_far  + abs(g_2_all{c}'*Theta*G_all*w_k(:, j))^2;
            end
        end

        % Signal power for current cluster
        A_n = abs(g_1_all{c}'*Theta*G_all*w_k(:, c))^2;
        A_f = abs(g_2_all{c}'*Theta*G_all*w_k(:, c))^2;

        B_n = inter_cluster_near + noise;
        B_f = inter_cluster_far + noise;

        % Compute SINR-based fraction
        gain_f = A_f / B_f;  % or A_n/B_n if using near user? check your formula
        alpha_f(c) = (para.R_min_f / (1 + para.R_min_f)) * (1 + 1/gain_f);
        alpha_f_temp = (para.R_min_f / (1 + para.R_min_f)) * (1 + 1/gain_f);


        % Ensure alpha_f in [0,1]
                % Ensure alpha_f in (0,0.9]
        if alpha_f_temp <= 0
            alpha_f(c) = 0.4;      % small positive minimum
        elseif alpha_f_temp > 0.9
            alpha_f(c) = 0.9;       % upper bound
        else
            alpha_f(c) = alpha_f_temp;
        end

        % Ensure alpha_n + alpha_f = 1
        alpha_n(c) = 1 - alpha_f(c);
    end
end

function [H, g, f] = generate_channel(para, BS_array, RIS_array)
    % Rician factor
    epsilon = para.rician;

    BS_loc = para.BS_loc;
    userloc = para.userloc;   % 3 x Nclusters x 3 (d,elev,az)
    Nclusters = para.K;       % number of clusters
    Kusers = para.K_u;        % users per cluster (3)
    
    % ============================
    % 1. BS -> RIS channel
    % ============================
    H_NLOS = 1/sqrt(2) * (randn(para.N, para.M) + 1i*randn(para.N, para.M));
    a_BR = steering_vector(BS_array, -BS_loc(2), -BS_loc(3));
    a_RB = steering_vector(RIS_array, BS_loc(2), BS_loc(3));
    H_LOS = a_RB * a_BR.';
   
    path_loss = sqrt(10.^(-para.pathloss(BS_loc(1))/10));
    H = path_loss * ( sqrt(epsilon/(epsilon+1)) * H_LOS + sqrt(1/(epsilon+1)) * H_NLOS );

    % ============================
    % 2. RIS -> Users channel
    % ============================
    g = zeros(para.N, Nclusters, Kusers);
    g_NLOS = 1/sqrt(2) * (randn(para.N, Nclusters, Kusers) + 1i*randn(para.N, Nclusters, Kusers));
    g_LOS = zeros(para.N, Nclusters, Kusers);

    for c = 1:Nclusters
        for k = 1:Kusers
            g_LOS(:,c,k) = steering_vector(RIS_array, userloc(k,c,2), userloc(k,c,3));
            d = userloc(k,c,1); 
            pl = sqrt(10.^(-para.pathloss(d)/10));
            g(:,c,k) = pl * ( sqrt(epsilon/(epsilon+1)) * g_LOS(:,c,k) + sqrt(1/(epsilon+1)) * g_NLOS(:,c,k) );
        end
    end

        % ============================
        % 3. BD -> Users channel (f)
        % Using stored distances from para.distances
        % ============================
        
        f = zeros(Nclusters, 2); % complex channel: BD->User1 and BD->User2
        
        for c = 1:Nclusters
            % Access distances from your struct
            dist_user1 = para.BSTdist.(['Cluster' num2str(c)]).User1_BD;
            dist_user2 = para.BSTdist.(['Cluster' num2str(c)]).User2_BD;
            cluster_dists = [dist_user1, dist_user2];
            disp(cluster_dists);
        
            for k = 1:2
                dist = cluster_dists(k);  % use precomputed distance
                
                % Path loss (convert from dB to linear)
                pl = sqrt(10^(-para.pathloss(dist)/10));
                
                % Small-scale Rayleigh fading
                f(c, k) = pl * (randn + 1i*randn) / sqrt(2);
            end
        end
end

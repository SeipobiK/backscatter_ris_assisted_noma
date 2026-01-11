function [angles, users,distances] = calculate_angles()
    % ============================
    % Fixed BS and RIS positions
    % ============================
BS = [0, 10, 20];
RIS = [0, 30, 20];

% ============================
% Hard-coded cluster centers
% ============================
cluster_centers = [0, 25, 0;
                   0, 35, 0;
                   5, 35, 0]; % new cluster center

radii = [5,5,5]; % radius for each cluster

% ============================
% Hard-coded users per cluster (3 users each)
% ============================

% ----- Cluster 1 (center: [0,25,0], radius: 5m) -----
users{1} = [ 0, 28, 0;   
             0, 26, 0;   
             2, 27, 0];  

% ----- Cluster 2 (center: [0,35,0], radius: 5m) -----
% users{2} = [ -1, 31, 0;   % User 1: North-East
%              -1, 40, 0;   % User 2: South-West
%              0, 31, 0];   % User 3: North-East
% 
% ----- Cluster 3 (center: [5,30,0], radius: 5m) -----
users{2} = [ 3, 30, 0;   % User 1: North-East   
             6, 30, 0;   % User 2: South-West
             4, 29, 0];  % User 3: South-East



% ----- Define cluster centers and radii -----
cluster_centers = [0 25 0;   % Cluster 1 center
                   0 35 0;   % Cluster 2 center
                   5 30 0];  % Cluster 3 center
    radius = 5;                  % Cluster radius (meters)
    
    num_clusters = size(cluster_centers,1);
    num_users_per_cluster = 3;
    
    % ----- Generate user positions -----
    % users = cell(num_clusters,1);
    % ----- Cluster 1 (center: [0,25,0], away direction +Y) -----
% users{1} = [ -1, 20, 0;     % left
%               0,  20, 0;     % center
%               1,  20, 0];    % right
% 
% % ----- Cluster 2 (center: [0,35,0], away direction +Y) -----
% users{2} = [ -1, 40, 0;
%               0,  40, 0;
%               1,  40, 0];
% 
% % ----- Cluster 3 (center: [5,30,0], away direction +X+Y) -----
% users{3} = [ 10, 31, 0;
%               9, 32, 0;
%               7, 26, 0];

    Nclusters = length(users);
    disp(Nclusters);
    Kusers = size(users{1},1);

    % ============================
    % 1. AoD at BS (BS -> RIS)
    % ============================
    delta_BS_RIS = RIS - BS;
    angles.BS_AoD.azimuth = atan2d(delta_BS_RIS(2), delta_BS_RIS(1));
    angles.BS_AoD.elevation = atan2d(delta_BS_RIS(3), sqrt(delta_BS_RIS(1)^2 + delta_BS_RIS(2)^2));
    angles.BS_AoD.distance = norm(delta_BS_RIS);

    % ============================
    % 2. AoA at RIS
    % ============================
    angles.RIS_AoA.azimuth = mod(angles.BS_AoD.azimuth + 180, 360);
    angles.RIS_AoA.elevation = -angles.BS_AoD.elevation;
    angles.RIS_AoA.distance = angles.BS_AoD.distance;

    % ============================
    % 3. AoD at RIS -> Users
    % ============================
    for c = 1:Nclusters
        for k = 1:Kusers
            user = users{c}(k,:);
            delta_RIS_user = user - RIS;

            azimuth = atan2d(delta_RIS_user(2), delta_RIS_user(1));
            elevation = atan2d(delta_RIS_user(3), sqrt(delta_RIS_user(1)^2 + delta_RIS_user(2)^2));

            angles.RIS_AoD.clusters(c).users(k).user = sprintf('C%d-U%d', c, k);
            angles.RIS_AoD.clusters(c).users(k).azimuth = azimuth;
            angles.RIS_AoD.clusters(c).users(k).elevation = elevation;
            angles.RIS_AoD.clusters(c).users(k).distance = norm(delta_RIS_user);
        end
    end
    
    % Initialize struct array
    distances = struct();
   
    
    for c = 1:length(users)
        BD = users{c}(3, :); % BD position
        clusterDistances = struct();
        
        % Distance between each user and BD
        for u = 1:2
            user_pos = users{c}(u, :);
            d_user_BD = norm(user_pos - BD);
            clusterDistances.(['User' num2str(u) '_BD']) = d_user_BD;
        end
        
        % Optional: distance between users
        d_user1_user2 = norm(users{c}(1,:) - users{c}(2,:));
        clusterDistances.User1_User2 = d_user1_user2;
    
        % Store in main struct
        distances.(['Cluster' num2str(c)]) = clusterDistances;
    end

%     % 
%     % % % % % % % % % ============================
%     % % % % % % % % 4. Visualization
%     % % % % % % % ============================
    figure; hold on; grid on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('BS, RIS, and User Clusters with Circles');

    % Plot BS and RIS
    scatter3(BS(1),BS(2),BS(3),100,'r','filled','^'); text(BS(1),BS(2),BS(3)+1,'BS');
    scatter3(RIS(1),RIS(2),RIS(3),100,'b','filled','s'); text(RIS(1),RIS(2),RIS(3)+1,'RIS');

    colors = lines(Nclusters);
    for c = 1:Nclusters
        % Cluster center
        center = cluster_centers(c,:);
        scatter3(center(1), center(2), center(3), 80, colors(c,:),'d','filled');
        text(center(1), center(2), center(3)+0.5,sprintf('Cluster %d', c),'Color','k','FontSize',10,'FontWeight','bold');

        % Draw cluster circle in XY-plane
        theta = linspace(0, 2*pi, 100);
        x_circle = center(1) + radii(c) * cos(theta);
        y_circle = center(2) + radii(c) * sin(theta);
        z_circle = center(3) * ones(size(theta));
        plot3(x_circle, y_circle, z_circle, '--','Color', colors(c,:), 'LineWidth',1.2);

        % Users
        cluster_users = users{c};
        scatter3(cluster_users(:,1), cluster_users(:,2), cluster_users(:,3), 50, colors(c,:),'filled','o');

        % Label users
        for k = 1:Kusers
            text(cluster_users(k,1), cluster_users(k,2), cluster_users(k,3)+0.5, ...
                sprintf('C%d-U%d', c, k), 'Color', colors(c,:), 'FontSize', 8, 'FontWeight','bold');
        end
    end
    legend({'BS','RIS','Cluster Centers','Cluster Circles','Users'});
    view(45,25);
    hold off;
end

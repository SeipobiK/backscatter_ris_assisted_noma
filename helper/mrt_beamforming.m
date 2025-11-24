function [w_k] = mrt_beamforming(para, H_users)
% MRT Beamforming for K-Cluster NOMA System
% Inputs:
%   - para: struct containing system parameters
%   - H_users: cell array of size {K clusters} containing each cluster's channel [M x K_u]
% Outputs:
%   - w_k: Beamforming matrix [M x K] (one beam per cluster)

M = para.M;      % Number of BS antennas
K = para.K;      % Number of clusters
P_max = para.P_max;

% Optional: power allocation per cluster (uniform here)
alpha = ones(1,K)/K;  % sum(alpha) = 1
% alpha=[0.4,0.6];

% disp(alpha);





w_k = zeros(M, K);

for c = 1:K
    Hc = H_users{c};        % M x K_u (users in cluster c)
    
    % Combine user channels in the cluster (mean, or strongest user)
    h_eff = mean(Hc, 2);    % effective channel vector
    
    % Normalize
    h_eff = h_eff / norm(h_eff);
    
    % Scale by allocated power
    w_k(:,c) = sqrt(P_max * alpha(c)) * h_eff;
end

% Verify total power constraint
total_power = sum(diag(w_k' * w_k));
if total_power > P_max * 1.01
    warning('Power constraint violated. Normalizing...');
    w_k = w_k * sqrt(P_max / total_power);
end

end

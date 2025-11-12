function w_k = init_active_BF(H_eff, P_max, method)
%INIT_ACTIVE_BF Initialize active beamforming matrix (N x K)
%
%   INPUTS:
%       H_eff  : K x N effective channel matrix 
%                (each row = effective channel for one user, size 1xN)
%       P_max  : total BS transmit power
%       method : 'MRT', 'ZF', 'RZF', or 'RANDOM'
%
%   OUTPUT:
%       w_k    : N x K beamforming matrix (each column = beam for cluster)

    [K, N] = size(H_eff);   % K users, N antennas
    w_k = zeros(N, K);

    switch upper(method)
        case 'MRT'
            % MRT per user
            for c = 1:K
                h = H_eff(c,:).'; 
                if norm(h) > 0
                    w_k(:,c) = h / norm(h);
                end
            end

        case 'ZF'
            if K > N
                error('ZF infeasible: K > N');
            end
            W_zf = H_eff' / (H_eff * H_eff');   % N x K
            for c = 1:K
                if norm(W_zf(:,c)) > 0
                    W_zf(:,c) = W_zf(:,c) / norm(W_zf(:,c));
                end
            end
            w_k = W_zf;

        case 'RZF'
            sigma2 = 1e-9; % small reg param (replace with para.noise if available)
            alpha = sigma2 * K / P_max;
            A = H_eff' * H_eff + alpha * eye(N);
            W_rzf = A \ H_eff';
            for c = 1:K
                if norm(W_rzf(:,c)) > 0
                    W_rzf(:,c) = W_rzf(:,c) / norm(W_rzf(:,c));
                end
            end
            w_k = W_rzf;

        case 'RANDOM'
            W_rand = (randn(N,K) + 1j*randn(N,K))/sqrt(2);
            for c = 1:K
                W_rand(:,c) = W_rand(:,c) / norm(W_rand(:,c));
            end
            w_k = W_rand;

        otherwise
            error('Unknown method. Use MRT, ZF, RZF, or RANDOM.');
    end

    % Scale to satisfy total power P_max
    total_power = sum(vecnorm(w_k,2,1).^2);
    if total_power > 0
        scale = sqrt(P_max / total_power);
        w_k = w_k * scale;
    end
end

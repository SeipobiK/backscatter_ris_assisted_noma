function [Jt, Jr] = get_permutation_matrices(G, H)
% Compute permutation matrices Jt and Jr for non-diagonal RIS
%
% Input:
%   G : N x M BS-RIS channel matrix
%   H : K x N RIS-Users channel matrix
%
% Output:
%   Jt : N x N permutation matrix (sorts g)
%   Jr : N x N permutation matrix (sorts h)

    [N, M] = size(G);
    [K, N_check] = size(H);
    assert(N == N_check, 'RIS dimension mismatch between G and H');

    % --- Average channel magnitudes
    g_avg = mean(abs(G), 2);    % N x 1
    h_avg = mean(abs(H), 1).';  % N x 1 (transpose!)

    % --- Sort indices in ascending order
    [~, idx_g] = sort(g_avg, 'ascend');
    [~, idx_h] = sort(h_avg, 'ascend');

    % --- Build permutation matrices
    Jt = eye(N);
    Jt = Jt(idx_g, :);   % sorts g

    Jr = eye(N);
    Jr = Jr(idx_h, :);   % sorts h

end

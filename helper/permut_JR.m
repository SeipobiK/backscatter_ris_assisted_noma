function [J_r] = permut_JR(g)
    N=10;
    a = abs(g);
    J_r = zeros(N, N);
    [~, idx_g] = sort(a, 'ascend');
    for i = 1:N
        J_r(idx_g(i), i) = 1;
    end
    % J_r=eye(N);
end
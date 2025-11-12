function [J_t] = permut_JT(g)
    N=10;
    a = abs(g);
    J_t = zeros(N, N);
    [~, idx_g] = sort(a, 'ascend');
    for i = 1:N
        J_t(i,idx_g(i)) = 1;
    end
    J_t=eye(N);
end
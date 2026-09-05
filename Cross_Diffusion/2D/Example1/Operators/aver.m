function [ax, ay] = aver(M, N)
    ax = Axs(M, N); ay = Ays(M, N);
end

function g = Axs(M, N)
    s = ones(1, N * (M - 1)); s = cat(2, s, s);
    i = cat(2, 1 : N * (M - 1), 1 : N * (M - 1));
    j = zeros(1, N * (M - 1));
    for n = 1 : N
        for m = 1 : M - 1
            j(m + (M - 1) * (n - 1)) = m + M * (n - 1);
        end
    end
    j = cat(2, j, j + 1);
    g = sparse(i, j, s)/2;
end

function g = Ays(M, N)
    s = ones(1, M * (N - 1)); s = cat(2, s, s);
    i = cat(2, 1 : M * (N - 1), 1 : M * (N - 1));
    j = 1 : M * (N - 1);
    j = cat(2, j, j + M);
    g = sparse(i, j, s)/2;
end
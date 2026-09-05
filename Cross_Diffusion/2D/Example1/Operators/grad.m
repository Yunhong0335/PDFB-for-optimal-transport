function [gx, gy] = grad(M, N, dx)
    gx = Gxs(M, N, dx); gy = Gys(M, N, dx);
end

function g = Gxs(M, N, dx)
    s = ones(1, N * (M - 1)); s = cat(2, -s, s);
    i = cat(2, 1 : N * (M - 1), 1 : N * (M - 1));
    j = zeros(1, N * (M - 1));
    for n = 1 : N
        for m = 1 : M - 1
            j(m + (M - 1) * (n - 1)) = m + M * (n - 1);
        end
    end
    j = cat(2, j, j + 1);
    g = sparse(i, j, s)/dx;
end

function g = Gys(M, N, dx)
    s = ones(1, M * (N - 1)); s = cat(2, -s, s);
    i = cat(2, 1 : M * (N - 1), 1 : M * (N - 1));
    j = 1 : M * (N - 1);
    j = cat(2, j, j + M);
    g = sparse(i, j, s)/dx;
end
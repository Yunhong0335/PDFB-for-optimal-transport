function dm = grad(N, dx)
    s = ones(1, N - 1); s = cat(2, -s, s);
    i = cat(2, 1 : N - 1, 1 : N - 1);
    j = cat(2, 1 : N - 1, 2 : N);
    dm = sparse(i, j, s)/dx;
end
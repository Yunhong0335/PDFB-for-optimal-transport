function dm = grad(N, dx)
    s = ones(1, N); s = cat(2, -s, s);
    i = cat(2, 1 : N, N, 1 : N - 1);
    j = cat(2, 1 : N, 1 : N);
    dm = sparse(i, j, s)/dx;
end
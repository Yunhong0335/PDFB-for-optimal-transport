function Am = aver(N)
    s = ones(1, 2 * (N - 1));
    i = cat(2, 1 : N - 1, 1 : N - 1);
    j = cat(2, 1 : N - 1, 2 : N);
    Am = sparse(i, j, s)/2;
end
function m = mobility(u)
    m11 = u(:, 1);
    m12 = zeros(size(u(:, 1)));
    m22 = u(:, 2);
    m = [m11, m12, m12, m22];
end
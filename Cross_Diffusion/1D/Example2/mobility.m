function m = mobility(u)
    rho = 2e-4;
    
    m11 = u(:, 1).^3/3;
    m12 = u(:, 1).^2 .* u(:, 2)/2;
    m22 = u(:, 1) .* u(:, 2).^2 + rho * u(:, 2);
    m = [m11, m12, m12, m22];
end
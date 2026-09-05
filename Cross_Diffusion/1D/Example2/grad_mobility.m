function [m1, m2] = grad_mobility(u)
    rho = 2e-4;
    
    m11 = u(:, 1).^2;
    m12 = u(:, 1).* u(:, 2);
    m22 = u(:, 2).^2;
    m1 = [m11, m12, m12, m22];
    
    m11 = zeros(size(u(:, 1)));
    m12 = u(:, 1).^2/2;
    m22 = 2 * u(:, 1).* u(:, 2) + rho;
    m2 = [m11, m12, m12, m22];
end
function [m1, m2] = grad_mobility(u)
    m11 = ones(size(u(:, 1)));
    m12 = zeros(size(u(:, 1)));
    m22 = zeros(size(u(:, 1)));
    m1 = [m11, m12, m12, m22];
    
    m11 = zeros(size(u(:, 1)));
    m12 = zeros(size(u(:, 1)));
    m22 = ones(size(u(:, 1)));
    m2 = [m11, m12, m12, m22];
end
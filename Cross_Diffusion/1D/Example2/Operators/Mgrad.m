function G = Mgrad(u, Q)
    G = zeros(size(u));
    [M1, M2] = grad_mobility(u);
    G(:, 1) = M1(:, 1) .* Q(:, 1) + M1(:, 2) .* Q(:, 2) + M1(:, 3) .* Q(:, 3) + M1(:, 4) .* Q(:, 4);
    G(:, 2) = M2(:, 1) .* Q(:, 1) + M2(:, 2) .* Q(:, 2) + M2(:, 3) .* Q(:, 3) + M2(:, 4) .* Q(:, 4);
end
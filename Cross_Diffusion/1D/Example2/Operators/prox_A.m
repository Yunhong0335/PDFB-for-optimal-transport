function [p, m] = prox_A(p0, m0, pn, K)
    [~, n] = size(p0);
    p = p0; m = m0;
    for ell = 1 : n
        [p(:, ell), m(:, ell)] = proj(p0(:, ell), m0(:, ell), pn(:, ell), K);
    end
end

function [p, m] = proj(p0, m0, pn, K)
    Mp = length(p0); Mm = length(m0);
    x0 = cat(1, p0, m0);
    y0 = pn;
    y = (K * K')\(K * x0 - y0);
    x = x0 - K' * y;
    p = x(1 : Mp);
    m = x(Mp + 1 : Mp + Mm);
end
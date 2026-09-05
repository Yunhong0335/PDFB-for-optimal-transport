function E = grad_e(p, ep, Kx, Ky, Hes)
    alpha = 0.4;
    E = ep^2 * Bound(p, Kx, Ky, alpha) + p .* (p.^2 - 1) - Hes * p;
end

function E3 = Bound(p, Kx, Ky, alpha)
    p1 = Kx * p; p2 = Ky * p;
    pnorm = sqrt(p1.^2 + p2.^2);
    n1 = zeros(size(p1)); n2 = zeros(size(p2));
    for i = 1 : length(pnorm)
        if pnorm(i) > 0
            n1(i) = p1(i)/pnorm(i);
            n2(i) = p2(i)/pnorm(i);
        end
    end
    
    g = V(n1, n2, alpha);
    p1 = g.^2 .* p1; p2 = g.^2 .* p2;
    
    [g1, g2] = dV(n1, n2, alpha);
    q1 = (1 - n1.^2) .* g1 - n1 .* n2 .* g2;
    q2 = (1 - n2.^2) .* g2 - n1 .* n2 .* g1;
    q1 = g .* pnorm .* q1;
    q2 = g .* pnorm .* q2;
    
    p1 = p1 + q1; p2 = p2 + q2;
    E3 = Kx' * p1 + Ky' * p2;
end

function g = V(n1, n2, alpha)
%     g2 = 32 * n2.^6 - 48 * n2.^4 + 18 * n2.^2 - 1;
%     g2 = 2 * n2.^2 - 1;
    g2 = 4 * n2.^3 - 3 * n2;
    g = 1 + alpha * g2;
end

function [g1, g2] = dV(n1, n2, alpha)
    g2 = 12 * n2.^2 - 3;
    
    g1 = zeros(size(n1));
    g2 = alpha * g2;
end
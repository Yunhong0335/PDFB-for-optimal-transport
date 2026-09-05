function E = ENE(p, Kx, Ky, vare, alpha)
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
    
    E = 1/4 * (1 - p.^2).^2 + vare^2 * g.^2 .* pnorm.^2;
    E = sum(E);
end

function g = V(n1, n2, alpha)
%     g2 = 32 * n2.^6 - 48 * n2.^4 + 18 * n2.^2 - 1;
%     g2 = 2 * n2.^2 - 1;
    g2 = 4 * n2.^3 - 3 * n2;
    g = 1 + alpha * g2;
end
function [phi, p1, p2] = prox_H(phi, p1, p2)
    spsi = p1.^2 + p2.^2;
    ind = phi + spsi/2 > 0;
    itermax = 100;
    
    lam = 1e3 * ones(size(phi));
    lam = lam(ind);
    if ~isempty(ind)
%         length(lam)
        for ell = 1 : itermax
            f = (1 + lam).^2 .* (phi(ind) - lam) + spsi(ind)/2;
%             max(f)
            if max(abs(f)) < 1e-10
                break
            end
            df = 2 * (1 + lam) .* (phi(ind) - lam) - (1 + lam).^2;
            grad = f./df;
            lam = lam - grad;
        end
    end
    phi(ind) = phi(ind) - lam;
    p1(ind) = p1(ind)./(1 + lam);
    p2(ind) = p2(ind)./(1 + lam);
end
function eta = prox_E(eta0, sig, dt)
    % p - p0 + sig * dt * e(p/dt)
    i = 1 : length(eta0); s = ones(size(eta0));
    I = sparse(i, i, s);
    
    itermax = 100;
    eta = zeros(size(eta0));
    for ell = 1 : itermax
        F = (eta - eta0) + sig * dt * exp(eta/dt);
%         max(abs(F))
        if max(abs(F)) < 1e-12
%             ell
            break
        end
        GF = I + sig * sparse(i, i, exp(eta/dt));
        eta = eta - GF\F;
    end
end
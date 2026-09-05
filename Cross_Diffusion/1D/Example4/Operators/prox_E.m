function eta = prox_E(eta0, sig, dt)
    % p - p0 + sig * e(p/dt)
    i = 1 : length(eta0);
    
    itermax = 100;
    eta = eta0;
    for ell = 1 : itermax
        F = (eta - eta0) + sig * exp(eta/dt);
%         max(abs(F))
        if max(abs(F)) < 1e-12
%             ell
            break
        end
        GF = 1 + sig * exp(eta/dt)/dt;
        eta = eta - F./GF;
    end
end
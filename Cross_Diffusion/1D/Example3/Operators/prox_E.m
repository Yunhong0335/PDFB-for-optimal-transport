function eta = prox_E(eta0, tau, epi)
    % p - p0 + tau * epi^2/p^3 - tau * epi^8/p^9
    itermax = 100;
    eta = eta0;
    for ell = 1 : itermax
        F = (eta.^10 - eta0 .* eta.^9) + tau * epi^2 * eta.^6 - tau * epi^8;
        if max(abs(F)) < 1e-12
            break
        end
        GF = 10 * eta.^9 - 9 * eta0 .* eta.^8 + tau * 6 * epi^2 * eta.^5;
        eta = eta - F./GF;
    end
end
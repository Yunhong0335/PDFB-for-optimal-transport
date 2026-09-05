function reference = IEsolver(M, T)
    R = 5; dx = 2 * R/M; t = 1; dt = t/T;
    space = linspace(-R + dx/2, R - dx/2, M)';

    % equation parameters
    p0 = exp(-(space - 1/2).^2/2)/sqrt(2 * pi);
    e0 = exp(-(space + 1/2).^2/2)/sqrt(2 * pi);

    % recorders
    itermax = 1000;
    rho = zeros(M, T + 1);
    eta = zeros(M, T + 1);
    rho(:, 1) = p0;
    eta(:, 1) = e0;

    % matrix initialize
    Gmat = grad(M, dx); Amat = aver(M);
    i = 1 : M; s = ones(M, 1); I = sparse(i,i,s);
    tic
    for tout = 1 : T
        p0 = rho(:, tout);
        e0 = eta(:, tout);

        p = p0; e = e0;
        for i = 1 : itermax
            Re = e; Rp = p;
            F = (p - p0) + dt * Gmat' * ((Amat * (2 * p + Re)) .* (Gmat * p) + (Amat * p) .* (Gmat * Re));
            G = (e - e0) + dt * Gmat' * ((Amat * (2 * e + Rp)) .* (Gmat * e) + (Amat * e) .* (Gmat * Rp));

            DF = eye(M) + dt * (2 * Gmat' * diag(Gmat * p) * Amat + Gmat' * diag(Amat * (2 * p + Re)) * Gmat + Gmat' * diag(Gmat * Re) * Amat);
            GF = eye(M) + dt * (2 * Gmat' * diag(Gmat * e) * Amat + Gmat' * diag(Amat * (2 * e + Rp)) * Gmat + Gmat' * diag(Gmat * Rp) * Amat);

            p = p - DF\F; e = e - GF\G;
            if max(F) < 1e-15 && max(G) < 1e-15
%                 [max(F), max(G)]
                break
            end
        end
        rho(:, tout + 1) = p;
        eta(:, tout + 1) = e;
        tout
    end
    reference = [rho(:, end), eta(:, end)];
    toc
end
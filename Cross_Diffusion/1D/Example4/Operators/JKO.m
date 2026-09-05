function sol = JKO(M, T)
    R = 5; dx = 2*R/M;
    space = linspace(-R + dx/2, R - dx/2, M)';
    t = 1; dt = t/T;

    % matrix initialize
    Dmat = -grad(M, dx)'; Amat = aver(M);
    i = 1 : M; s = ones(M, 1); I = sparse(i,i,s);
    K = cat(2, I, Dmat);

    mu1 = 1/sqrt(2 * pi) * exp(-(space - 1/2).^2/2);
    mu2 = 1/sqrt(2 * pi) * exp(-(space + 1/2).^2/2);
    u0 = cat(2, mu1, mu2);
    Rec1 = zeros(M, T + 1); Rec1(:, 1) = u0(:, 1);
    Rec2 = zeros(M, T + 1); Rec2(:, 1) = u0(:, 2);

    tic
    for tout = 1 : T
        % initialization
        u = u0; m = zeros(M - 1, 2);
        Q = zeros(M, 4); q = zeros(M, 2);
        nu = zeros(M, 2);

        % optimization
        itermax = 20000;
        L = 1; tau = 2; sig = L^2/tau;

        zu = u; zm = m;
        for ell = 1 : itermax
            Ru = u; Rm = m;
            
            Q = Q + sig * Mfun((zu + u0)/2);
            q = q + sig * Amat' * zm;
            [Q, q] = prox_C(Q, q);
            nu = nu + sig * zu;
            for i = 1 : 2
                nu(:, i) = prox_E(nu(:, i), sig, dt);
            end

            u = u - tau * (nu + Mgrad((u + u0)/2, Q)/2);
            m = m - tau * (Amat * q);
            [u, m] = prox_A(u, m, u0, K);

            zu = 2 * u - Ru;
            zm = 2 * m - Rm;

            % recorder update
            res = norm(u - Ru)/norm(Ru)
            if res < 1e-8
                res
                break
            end
        end
        u0 = u;
        Rec1(:, tout + 1) = u(:, 1);
        Rec2(:, tout + 1) = u(:, 2);
        tout
    end
    sol = [Rec1(:, end), Rec2(:, end)];
    toc
end
clear;
addpath('Operators')

% spatiotemporal size
M = 1024; R = 1; dx = 2*R/M;
space = linspace(-R + dx/2, R - dx/2, M)';
t = 0.8/24; T = 10; dt = t/T;
time = linspace(0, t, T + 1);

c = 2e-2;

% matrix initialize
Dmat = -grad(M, dx)'; Amat = aver(M);
i = 1 : M; s = ones(M, 1); I = sparse(i,i,s);
K = cat(2, I, Dmat);
Hes = c * (Dmat * Dmat');

mu1 = max(0, 1/4 - space.^2);
mu2 = 0.9/2 * max(0, (1 - 1e-4) - tanh(50 * (abs(space) - 0.1)));
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
    L = 1; tau = 5; sig = L^2/tau;
    
    zu = u; zm = m;
    for ell = 1 : itermax
        Ru = u; Rm = m;
        
        Mob = Mfun((zu + u0)/2);
        Q = Q + sig * Mfun((zu + u0)/2);
        q = q + sig * Amat' * zm;
        [Q, q] = prox_C(Q, q);
        
        nu = nu + sig * zu;
        ex = (I + dt/sig * Hes)\(nu(:, 1)/sig);
        nu(:, 1) = nu(:, 1) - sig * ex;
        nu(:, 2) = prox_E(nu(:, 2), sig, dt);
        
        u = u - tau * (nu + Mgrad((u + u0)/2, Q)/2);
        m = m - tau * Amat * q;
        for i = 1 : 2
            [u(:, i), m(:, i)] = prox_A_new(u(:, i), m(:, i), u0(:, i), [0, inf], K);
        end
        
        zu = 2 * u - Ru;
        zm = 2 * m - Rm;
        
        % recorder update
        res = norm(u - Ru)/norm(Ru)
        if res < 1e-6
            break
        end
    end
    u0 = u;
    Rec1(:, tout + 1) = u(:, 1);
    Rec2(:, tout + 1) = u(:, 2);
end
toc

rmpath('Operators')

figure
hold on
for ell = 1 : 5 : T + 1
    plot(space, Rec1(:, ell), 'Linewidth', 2, 'Color', [0, 1 - ell/(T + 1), ell/(T + 1)])
end
for ell = 1 : 5 : T + 1
    plot(space, Rec2(:, ell), 'Linewidth', 2, 'Color', [1, 1 - ell/(T + 1), ell/(T + 1)])
end
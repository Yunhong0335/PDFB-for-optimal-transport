clear;
addpath('Operators')

% spatiotemporal size
M = 400; R = 1; dx = 2*R/M;
space = linspace(-R + dx/2, R - dx/2, M)';
t = 0.5; T = 5000; dt = t/T;
time = linspace(0, t, T + 1);
nuc = 1.2; vare = 0.01;

% matrix initialize
Dmat = -grad(M, dx)'; Amat = aver(M);
i = 1 : M; s = ones(M, 1); I = sparse(i,i,s);
K = cat(2, I, Dmat);
Hes = Dmat * Dmat';

mu1 = 1 - (1 + cos(pi * space))/4;
mu2 = (1 + vare) * ones(size(space));
u0 = cat(2, mu1, mu2);
Rec1 = zeros(M, T + 1); Rec1(:, 1) = u0(:, 1);
Rec2 = zeros(M, T + 1); Rec2(:, 1) = u0(:, 2);

tic
for tout = 1 : T
    % initialization
    u = u0; m = zeros(M - 1, 2);
    Q = zeros(M, 4); q = zeros(M, 2);
    nu = zeros(M, 1);
    
    % optimization
    itermax = 20000;
    L = 1; tau = 0.02; sig = L^2/tau;
    
    zu = u; zm = m;
    for ell = 1 : itermax
        Ru = u; Rm = m;
        
        Mob = Mfun((zu + u0)/2);
        Q = Q + sig * Mfun((zu + u0)/2);
        q = q + sig * Amat' * zm;
        [Q, q] = prox_C_new(Q, q, Mob);
        
        nu = nu + sig * (zu(:, 2) - zu(:, 1));
        ex = prox_E(nu/sig, dt/sig, vare);
        nu = real(nu - sig * ex);
        
        u = u - tau * Mgrad((u + u0)/2, Q)/2 - tau * cat(2, dt * nuc * Hes * u(:, 1) - nu, dt * Hes * u(:, 2) + nu);
        m = m - tau * Amat * q;
        [u, m] = prox_A(u, m, u0, K);
        
        zu = 2 * u - Ru;
        zm = 2 * m - Rm;
        
        % recorder update
        res = norm(u - Ru)/norm(Ru)
        if res < 1e-9
            break
        end
    end
    [tout, ell]
    
    u0 = u;
    Rec1(:, tout + 1) = u(:, 1);
    Rec2(:, tout + 1) = u(:, 2);
end
toc
rmpath('Operators')

T = 5000;
figure
energy = zeros(1, T + 1);
for ell = 1 : T + 1
    D1 = nuc/2 * (Dmat' * Rec1(:, ell)).^2;
    D1 = sum(D1);
    
    D2 = 1/2 * (Dmat' * Rec2(:, ell)).^2;
    D2 = sum(D2);
    
    D3 = Rec2(:, ell) - Rec1(:, ell);
    D3 = vare^8./(8*D3.^8) - vare^2./(2*D3.^2);
    P = sum(D3);
    energy(ell) = (D1 + D2 + P) * dx;
end


plot(time, energy, 'Linewidth', 2, 'Color', [0, 0, 1])

T = 1500;
frame = [1, 10, 301, 1501];
figure
hold on
for ell = 0 : 3
    plot(space, Rec1(:, frame(ell + 1)), 'Linewidth', 2, 'Color', [0, 1 - ell/3, ell/3])
end
for ell = 0 : 3
    plot(space, Rec2(:, frame(ell + 1)), 'Linewidth', 2, 'Color', [1, 1 - ell/3, ell/3])
end
ylim([0.3, 1.2])
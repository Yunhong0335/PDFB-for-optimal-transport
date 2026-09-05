clear;
addpath('Operators')

% spatiotemporal size
M = 400; R = 4; dx = 2 * R/M;
t = 10; T = 100; dt = t/T;

time = linspace(0, t, T + 1);
space = linspace(-R + dx/2, R - dx/2, M)';

% equation parameters
p0 = sum(exp(-1/2 * max((space.^2 - 1), 0))) * ones(size(space)) * dx/8;

% recorders
itermax = 3000; conver = zeros(itermax, T);
rho = zeros(M, T + 1); rho(:, 1) = p0;

% matrix initialize
Dmat = -grad(M, dx)'; Amat = aver(M);
i = 1 : M; s = ones(M, 1); I = sparse(i,i,s);
% Pcond = I + Dmat * Dmat';
% Pcond = ichol(Pcond, struct('michol', 'on'));
Kp = cat(2, I, Dmat);

tic
for tout = 1 : T
    p0 = rho(:, tout);
    
    % grid initialization
    m = zeros(M - 1, 1); p = p0; xi = zeros(M, 1);
    phi = zeros(M, 1); psi = zeros(M, 1);
    
    % primal dual optimization
    L = 1;
    tau = 0.2; sig = L^2/tau;
    
    Tm = m; Tp = p;
    for ell = 1 : itermax
        Rm = m; Rp = p;
        
        % phi, psi
        Ap = (p0 + p)/2;
        
        Td = Ap .* (1 - Ap) - (1/2 - Ap) .* p;
        phi = phi + sig * (1/2 - Ap) .* Tp + sig * Td;
        psi = psi + sig * Amat' * Tm;
        [phi, psi] = prox_H(phi, psi);
        
        xi = xi + sig * Tp;
        xi = prox_E(xi, sig, dt);
        
        % m, rho
        m = m - tau * Amat * psi;
        
        Ap = (p0 + p)/2;
        gradp = phi .* (1/2 - Ap) + dt * grad_e(p, space);
        p = p - tau * gradp - tau * xi;
        
        % proxA
        [p, m] = prox_A(p, m, p0, [0, 1], Kp);
        
        Tm = 2 * m - Rm;
        Ap = (p0 + p)/2;
        Tgradp = phi .* (1/2 - Ap) + dt * grad_e(p, space);
        Tp = 2 * p - Rp - tau * (Tgradp - gradp);
        
        % recorder update
        res = sqrt(sum((p - Rp).^2)/sum(Rp.^2));
        conver(ell, tout) = res;
        if res < 1e-7
            break
        end
    end
    ell
    rho(:, tout + 1) = p;
end
toc

rmpath('Operators')

figure
hold on
for ell = 1 : 10: T + 1
    plot(space, rho(:, ell), 'Linewidth', 2, 'Color', [0, 1 - ell/(T + 1), ell/(T + 1)])
end
pinf = exp(-1/2 * max((space.^2 - 1), 0));
plot(space, pinf, '--', 'Linewidth', 3, 'Color', [0.75, 0.75, 0.75])
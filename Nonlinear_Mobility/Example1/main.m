clear;

addpath('Operators')

% spatiotemporal size
M = 200; R = 1; dx = 2 * R/M;
t = 1e-2; T = 20; dt = t/T;

time = linspace(0, t, T + 1);
space = linspace(-R + dx/2, R - dx/2, M)';

% equation parameters
t0 = 1e-3;
p0 = t0^(-1/3) * max(0, (3/16)^(1/3) - t0^(-2/3) * space.^2/12);

% recorders
itermax = 10000; conver = zeros(itermax, T);
pos = zeros(T + 1, 1); vol = zeros(T + 1, 1);

rho = zeros(M, T + 1); rho(:, 1) = p0;

energy = zeros(1, T + 1);
energy(1) = sum(p0.^2/2) * dx;
pos(1) = min(p0);
vol(1) = sum(p0) * dx;
    
% matrix initialize
Dmat = -grad(M, dx)'; Amat = aver(M);
i = 1 : M; s = ones(M, 1); I = sparse(i,i,s);
Kp = cat(2, I, Dmat);

tic
for tout = 1 : T
    p0 = rho(:, tout);
    
    % grid initialization
    m = zeros(M - 1, 1); p = zeros(M, 1); 
    phi = zeros(M, 1); psi = zeros(M, 1);
    
    % primal dual optimization
    L = 1;
    tau = 1; sig = L^2/tau;
    
    Tm = m; Tp = p;
    for ell = 1 : itermax
        Rm = m; Rp = p;
        
        Am = Amat' * Tm;
        Ap = (p0 + Tp)/2;
        
        % phi, psi
        phi = phi + sig * Ap;
        psi = psi + sig * Am;
        [phi, psi] = prox_H(phi, psi);
        
        % m, rho
        m = m - tau * (Amat * psi);
        p = p - tau * phi/2 - tau * dt * grad_e(p);
        
        % proxA
        [p, m] = prox_A(p, m, p0, [0, inf], Kp);

        Tm = m + (m - Rm);
        Tp = p + (p - Rp) - dt * tau * (grad_e(p) - grad_e(Rp));

        % recorder update
        res = sqrt(sum((p - Rp).^2)/sum(Rp.^2))
        conver(ell, tout) = res;
        if res < 1e-5
            break
        end
    end
%     ell
    pos(tout + 1) = min(p);
    vol(tout + 1) = sum(p) * dx;
    rho(:, tout + 1) = p;
    energy(tout + 1) = sum(p.^2/2) * dx;
end
toc

rmpath('Operators')

figure
hold on
for ell = 1 : 4: T + 1
    plot(space, rho(:, ell), 'Linewidth', 2, 'Color', [0, 1 - ell/(T + 1), ell/(T + 1)])
end
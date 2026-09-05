clear;
addpath('Operators')

% spatiotemporal size
M = 200; R = 1; dx = 2 * R/M;
t = 0.1; T = 100; dt = t/T;

time = linspace(0, t, T + 1);
space = linspace(-R + dx/2, R - dx/2, M)';

% equation parameters
p0 = 0.8 - cos(pi * space) + 0.25 * cos(2 * pi * space);

% recorders
itermax = 1000;
pos = zeros(T, 1); vol = zeros(T, 1);
nos = zeros(T, 1);

rho = zeros(M, T + 1); rho(:, 1) = p0;

% matrix initialize
Gmat = grad(M, dx); Amat = aver(M);
i = 1 : M; s = ones(M, 1); I = sparse(i,i,s);
Bes = I + Gmat' * Gmat;
Kp = cat(1, cat(2, I, I), cat(2, I, -Gmat' * Gmat));

tic
for tout = 1 : T
    p0 = rho(:, tout);
    
    % grid initialization
    m = zeros(M, 1); p = p0; 
    phi = zeros(M, 1); psi = zeros(M, 1);
    
    % primal dual optimization
    L = 1;
    tau = 0.03; sig = L^2/tau;
    
    Tm = m; Tp = p;
    for ell = 1 : itermax
        Rm = m; Rp = p;
        
        Am = Amat' * Tm;
        Ap = (p0 + Tp)/2;
        
        % phi, psi
        phi = phi + sig * Ap.^3;
        psi = psi + sig * Am;
        [phi, psi] = prox_H(phi, psi);
        
        % m, rho
        m = m - tau * (Amat * psi);
        p = p - tau * 3 * phi .* ((p0 + p)/2).^2/2 - tau * dt * grad_e(p, Bes);

        % proxA            
        [p, Ca] = prox_A(p, p0 + Gmat' * m, [0, inf], Kp, Bes);
        m = m + Gmat * Ca;
        
        Tm = m + (m - Rm);
        Tp = p + (p - Rp) + dt * tau * (grad_e(p, Bes) - grad_e(Rp, Bes));

        % recorder update
        res = sqrt(sum((p - Rp).^2)/sum(Rp.^2));
        res_m = sqrt(sum((m - Rm).^2)/M);
        [res, res_m]
        if res < 1e-7 && res_m < 1e-5
            break
        end
    end
    pos(tout) = min(p0);
    nos(tout) = max(p0);
    vol(tout) = sum(p0) * dx;
    rho(:, tout + 1) = p;
end
toc

rmpath('Operators')

figure
hold on
for ell = 1 : 20 : T + 1
    plot(space, rho(:, ell), 'Linewidth', 2, 'Color', [0, 1 - ell/(T + 1), ell/(T + 1)])
end
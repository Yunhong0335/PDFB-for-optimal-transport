clear;

addpath('Operators')

% initialize
t = 1; T = 1000; dt = t/T;
time = linspace(0, t, T + 1);

R = 0.5; M = 128; dx = R/M;
X = linspace(-R + dx/2, R - dx/2, M);

H = 0.5; N = 128; dh = H/N;
Y = linspace(-H + dh/2, H - dh/2, N);

% vectorize
vare = 0.1/4; F = ones(M, N);
for m = 1 : M
    for n = 1 : N
        F(m, n) = -tanh((sqrt(X(m)^2 + Y(n)^2) - 0.3)/vare);
    end
end

p0 = zeros(M * N, 1);
for m = 1 : M
    for n = 1 : N
        p0(m + M * (n - 1)) = F(m, n);
    end
end

% matrix initialize
[Dx, Dy] = Grad(M, N, dx, dh); Dx = -Dx'; Dy = -Dy';
[Ax, Ay] = Aver(M, N);
Kx = -Ax' * Dx'; Ky = -Ay' * Dy';

i = 1 : M * N; s = ones(M * N, 1); I = sparse(i,i,s);
Kp = cat(2, I, Dx, Dy);
Hes = Dx * Dx' + Dy * Dy';

% recorders
itermax = 3000;
rho = zeros(M * N, T + 1); rho(:, 1) = p0;
conver = zeros(itermax, T);
energy = zeros(T + 1, 1);

vare = 0.01; alpha = 0.4;
energy(1) = ENE(p0, Kx, Ky, vare, alpha) * dx^2;

for tout = 1 : T
    p0 = rho(:, tout);
    
    % grid initialization
    m1 = zeros((M - 1) * N, 1);
    m2 = zeros(M * (N - 1), 1);
    
    p = p0;
    
    phi = zeros(M * N, 1);
    nu = zeros(M * N, 1);
    
    p1 = zeros(M * N, 1);
    p2 = zeros(M * N, 1);
    
    % primal dual optimization
    L = 1; tau = 2; sig = L^2/tau;
    
    beta = 1e-4;
    
    Tm1 = m1; Tm2 = m2; Tp = p;
    for ell = 1 : itermax
        Rm1 = m1; Rm2 = m2;
        Rp = p; Rp1 = p1; Rp2 = p2;
        
        % m, rho
        m1 = Tm1 - tau * Ax * p1;
        m2 = Tm2 - tau * Ay * p2;
        
        Ap = (p + p0)/2;
        gradp = -Ap .* phi + nu + dt * grad_e(p, vare, Kx, Ky, alpha * vare^2 * Hes);
        
        p = Tp - tau * gradp;
        
        % proxA
        [p, m1, m2] = prox_A(p, m1, m2, p0, [-1, 1], Kp);
        
        % phi, psi
        Am1 = Ax' * m1; Am2 = Ay' * m2;
        Ap = (p0 + p)/2;
        
        phi = phi + sig * (1 - Ap.^2);
        p1 = p1 + sig * Am1;
        p2 = p2 + sig * Am2;
        
        [phi, p1, p2] = prox_H(phi, p1, p2);
                
        nu = nu + sig * p;
        eta = prox_E2(nu/sig, beta * vare^2 * dt/sig * (Hes * Hes) + alpha * vare^2 * dt/sig * Hes);
        nu = real(nu - sig * eta);
        
        Tm1 = m1 - tau * Ax * (p1 - Rp1);
        Tm2 = m2 - tau * Ay * (p2 - Rp2);
        
        Ap = (p0 + p)/2;
        Tgradp = -Ap .* phi + nu + dt * grad_e(p, vare, Kx, Ky, alpha * vare^2 * Hes);
        Tp = p - tau * (Tgradp - gradp);

        % recorder update
        res = sqrt(sum((p - Rp).^2)/sum(Rp.^2));
%         conver(ell, tout) = res;
        if res < 1e-5
            break
        end
    end
    ell
    rho(:, tout + 1) = p;
    energy(tout + 1) = ENE(p, Kx, Ky, vare, alpha) * dx^2;
end

rmpath('Operators')

M = 128; N = 128;
F = zeros(M, N);
for m = 1 : M
    for n = 1 : N
    	F(m, n) = rho(m + M * (n - 1), end);
    end
end
Heatmap(F)
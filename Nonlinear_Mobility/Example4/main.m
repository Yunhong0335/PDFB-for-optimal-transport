clear;
addpath('Operators')

% initialize
M = 64; N = 64; dx = 1/M;

X = dx * (0 : M - 1) + dx/2;
Y = dx * (0 : N - 1) + dx/2;
[X, Y] = ndgrid(X, Y);

% img = rand(M, N);
load('initial')

t = 0.02; T = 20; dt = t/T;
time = linspace(0, t, T + 1);

% grid initilize
% F = img;
% p0 = zeros(M * N, 1);
% for m = 1 : M
%     for n = 1 : N
%         p0(m + M * (n - 1)) = F(m, n);
%     end
% end

p0 = initial;

% recorders
itermax = 1700; % pos = zeros(T, 1); vol = zeros(T, 1);
rho = zeros(M * N, T + 1); rho(:, 1) = p0;

conver = zeros(itermax, T);

% matrix initialize
[Gxmat, Gymat] = Grad(M, N, dx);
[Axmat, Aymat] = Aver(M, N);

i = 1 : M * N;
s = ones(M * N, 1);
Bes = sparse(i,i,s) + Gxmat' * Gxmat + Gymat' * Gymat;
Hes = Gxmat' * Gxmat + Gymat' * Gymat;
Kp = cat(1, cat(2, sparse(i,i,s), sparse(i,i,s)), cat(2, sparse(i,i,s), -Gxmat' * Gxmat - Gymat' * Gymat));

iter_number = zeros(1, T);

tic
for tcont = 1 : T
    p0 = rho(:, tcont);
    
    p = p0;
    m1 = zeros((M - 1) * N, 1); m2 = zeros(M * (N - 1), 1);
    phi = zeros(M * N, 1);
    p1 = zeros(M * N, 1); p2 = zeros(M * N, 1);
    
    % primal dual optimization
    L = 1; tau = 20; sig = L^2/tau;
    
    Tm1 = m1; Tm2 = m2; Tp = p;
    for ell = 1 : itermax
        Rm1 = m1; Rm2 = m2; Rp = p;
        
        Am1 = Axmat' * Tm1; Am2 = Aymat' * Tm2;
        Ap = (p0 + Tp)/2;
        
        % phi, psi
        phi = phi + sig * (1 - Ap.^2);
        p1 = p1 + sig * Am1;
        p2 = p2 + sig * Am2;
        [phi, p1, p2] = prox_H(phi, p1, p2);
        
        % m, rho
        m1 = m1 - tau * Axmat * p1;
        m2 = m2 - tau * Aymat * p2;
        p = p + tau * phi .* Ap - tau * dt * grad_e(p, Hes);
        
        % proxA
        [p, Ca] = prox_A(p, p0 + Gxmat' * m1 + Gymat' * m2, [-1, 1], Kp, Bes);
        m1 = m1 + Gxmat * Ca;
        m2 = m2 + Gymat * Ca;
        
        Tm1 = 2 * m1 - Rm1;
        Tm2 = 2 * m2 - Rm2;
        Tp = 2 * p - Rp - tau * (dt * grad_e(p, Hes) - phi .* (p + p0)/2) + tau * (dt * grad_e(Rp, Hes) - phi .* (Rp + p0)/2);

        % recorder update
        res = sqrt((sum((p - Rp).^2) + sum((m1 - Rm1).^2) + sum((m2 - Rm2).^2))/(sum(Rp.^2) + sum(Rm1.^2) + sum(Rm2.^2)));
        conver(ell, tcont) = res;
        res
        if res < 1e-5
            break
        end
    end
    iter_number(tcont) = ell;
    rho(:, tcont + 1) = p;
end
toc

rmpath('Operators')

display_p = zeros(M, N, T + 1);
for tcont = 1 : T + 1
    for m = 1 : M
        for n = 1 : N
            display_p(m, n, tcont) = rho(m + M * (n - 1), tcont);
        end
    end
end
Heatmap(display_p(:, :, end))
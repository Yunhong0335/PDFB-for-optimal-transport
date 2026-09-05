clear;
addpath('Operators')

% spatiotemporal size
M = 64; N = 64; dx = 10/M;

X = dx * (-M/2 : M/2 - 1) + dx/2;
Y = dx * (-N/2 : N/2 - 1) + dx/2;
[X, Y] = ndgrid(X, Y);

t = 1; T = 10; dt = t/T;
time = linspace(0, t, T + 1);

% matrix initialize
[Gxmat, Gymat] = grad(M, N, dx);
[Axmat, Aymat] = aver(M, N);
i = 1 : M * N; s = ones(M * N, 1); I = sparse(i,i,s);
K = cat(2, I, -Gxmat', -Gymat');

smu1 = 1/sqrt(2*pi) * exp(-(X - 1/2).^2/2 - Y.^2/2);
smu2 = 1/sqrt(2*pi) * exp(-(X + 1/2).^2/2 - Y.^2/2);
mu1 = zeros(M * N, 1);
mu2 = zeros(M * N, 1);
for ell = 1 : N
    mu1(1 + (ell - 1) * M : ell * M) = smu1(:, ell);
    mu2(1 + (ell - 1) * M : ell * M) = smu2(:, ell);
end
u0 = cat(2, mu1, mu2);
Rec1 = zeros(M * N, T + 1); Rec1(:, 1) = u0(:, 1);
Rec2 = zeros(M * N, T + 1); Rec2(:, 1) = u0(:, 2);

tic
for tout = 1 : T
    % grid initialization
    u = u0;
    m = zeros((M - 1) * N, 2); v = zeros(M * (N - 1), 2);
    phi = zeros(M * N, 4); p1 = zeros(M * N, 2); p2 = zeros(M * N, 2);
    nu = zeros(M * N, 2);
    
    % primal dual optimization
    itermax = 20000;
    L = 1; tau = 0.2; sig = L^2/tau;
    
    Tu = u; Tm = m; Tv = v;
    for ell = 1 : itermax
        Ru = u; Rm = m; Rv = v;
        
        phi = phi + sig * Mfun((u0 + Tu)/2);
        p1 = p1 + sig * Axmat' * Tm;
        p2 = p2 + sig * Aymat' * Tv;
        [phi, p1, p2] = prox_C(phi, p1, p2);
        
        nu = nu + sig * Tu;
        nu(:, 1) = prox_E(nu(:, 1), sig, dt);
        nu(:, 2) = prox_E(nu(:, 2), sig, dt);
        
        % m, rho
        u = u - tau * Mgrad((u + u0)/2, phi)/2 - tau * nu;
        m = m - tau * (Axmat * p1);
        v = v - tau * (Aymat * p2);
        
        % proxA
        [u(:, 1), m(:, 1), v(:, 1)] = prox_A(u(:, 1), m(:, 1), v(:, 1), u0(:, 1), [0, inf], K);
        [u(:, 2), m(:, 2), v(:, 2)] = prox_A(u(:, 2), m(:, 2), v(:, 2), u0(:, 2), [0, inf], K);
        
        Tu = 2 * u - Ru;
        Tm = 2 * m - Rm;
        Tv = 2 * v - Rv;
        
        % recorder update
        res = norm(u - Ru)/norm(Ru)
        if res < 1e-5
            break
        end
    end
    u0 = u;
    Rec1(:, tout + 1) = u(:, 1);
    Rec2(:, tout + 1) = u(:, 2);
end
toc
rmpath('Operators')

display_1 = zeros(M, N);
display_2 = zeros(M, N);
for m = 1 : M
    for n = 1 : N
        display_1(m, n) = Rec1(m + M * (n - 1), end);
        display_2(m, n) = Rec2(m + M * (n - 1), end);
    end
end
Heatmap(display_1, [0, 0.2])
Heatmap(display_2, [0, 0.2])
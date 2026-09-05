clear;
addpath('Operators')

% spatiotemporal size
M = 400; R = 1; dx = 2*R/M;
space = linspace(-R + dx/2, R - dx/2, M)';
t = 3; T = 30; dt = t/T;
time = linspace(0, t, T + 1);

% matrix initialize
Dmat = -grad(M, dx)'; Amat = aver(M);
i = 1 : M; s = ones(M, 1); I = sparse(i,i,s);
INull = sparse(M, M); DNull = sparse(M, M - 1);

K1 = cat(2, I, INull, Dmat, DNull);
K2 = cat(2, INull, I, DNull, Dmat);
K = cat(1, K1, K2); Kp = cat(2, I, Dmat);

Kold = cat(2, I, Dmat);

D = 0.2; C1 = 4; C2 = 2;
f = 2/5 * (1 - (4 * space/3).^2);
mu1 = max(0, f .* (1 - cos(8 * pi * space)/2));
mu2 = max(0, f .* (1 + cos(8 * pi * space)/2));
u0 = cat(2, mu1, mu2);

Rec1 = zeros(M, T + 1); Rec1(:, 1) = u0(:, 1);
Rec2 = zeros(M, T + 1); Rec2(:, 1) = u0(:, 2);

T_initial = 1;
tic
for tout = T_initial : T
    % initialization
    u = u0; m = zeros(M - 1, 2);
    Q = zeros(M, 4); q = zeros(M, 2);
    xi = zeros(M, 1);
    
    % optimization
    itermax = 30000;
    L = 1; tau = 0.1; sig = L^2/tau;
    
    zu = u; zm = m;
    for ell = 1 : itermax
        Ru = u; Rm = m;
        
        Q = Q + sig * Mfun((zu + u0)/2);
        q = q + sig * Amat' * zm;
        [Q, q] = prox_C_Newton(Q, q);
        
        s = u(:, 1) + u(:, 2);
        u = u - tau * Mgrad((u + u0)/2, Q)/2 - dt * tau * cat(2, D * s + C1 * space.^2/2, D * s + C2 * space.^2/2);
        m = m - tau * Amat * q;
        [u, m] = prox_A_new(u, m, u0, K);
%         for i = 1 : 2
%             [u(:, i), m(:, i)] = prox_A_old(u(:, i), m(:, i), u0(:, i), [0, inf], Kold);
%         end

        zu = 2 * u - Ru;
        zm = 2 * m - Rm;
        
        % recorder update
        res = norm(u - Ru)/norm(Ru)
        if res < 1e-8
            break
        end
    end
    u0 = u;
    Rec1(:, tout + 1) = u(:, 1);
    Rec2(:, tout + 1) = u(:, 2);
end
toc

rmpath('Operators')

energy = zeros(1, T + 1);
for ell = 1 : T + 1
    energy(ell) = sum(D/2 * (Rec1(:, ell) + Rec2(:, ell)).^2 + C1 * space.^2/2 .* Rec1(:, ell) + C2 * space.^2/2 .* Rec2(:, ell)) * dx;
end

figure
hold on
plot(space, Rec1(:, 1), 'Linewidth', 2, 'Color', [0, 0, 1])
plot(space, Rec1(:, 1) + Rec2(:, 1), 'Linewidth', 2, 'Color', [0, 1, 0])
ylim([0, 1])

figure
hold on
plot(space, Rec1(:, end), 'Linewidth', 2, 'Color', [0, 0, 1])
plot(space, Rec1(:, end) + Rec2(:, end), 'Linewidth', 2, 'Color', [0, 1, 0])
ylim([0, 1])

figure
hold on
plot(space, Rec1(:, tout) + Rec2(:, tout), 'Linewidth', 2, 'Color', [0, 1, 0])
plot(space, Rec1(:, tout), 'Linewidth', 2, 'Color', [0, 0, 1])

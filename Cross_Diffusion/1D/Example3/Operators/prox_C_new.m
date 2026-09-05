function [Q, q] = prox_C_new(Q0, q0, M0)
    [Lam, U] = Eig_de(M0);
    Lam(:, 1) = Q0(:, 1) .* U(:, 1) .* U(:, 1) + 2 * Q0(:, 2) .* U(:, 1) .* U(:, 3) + Q0(:, 4) .* U(:, 3) .* U(:, 3);
    Lam(:, 2) = Q0(:, 1) .* U(:, 2) .* U(:, 2) + 2 * Q0(:, 2) .* U(:, 2) .* U(:, 4) + Q0(:, 4) .* U(:, 4) .* U(:, 4);
    
    q1 = U(:, 1) .* q0(:, 1) + U(:, 3) .* q0(:, 2);
    q2 = U(:, 2) .* q0(:, 1) + U(:, 4) .* q0(:, 2);
    q = cat(2, q1, q2);
    
    for ell = 1 : 2
        [Lam(:, ell), q(:, ell)] = prox_H(Lam(:, ell), q(:, ell));
    end
    
    Q(:, 1) = Lam(:, 1) .* U(:, 1) .* U(:, 1) + Lam(:, 2) .* U(:, 2) .* U(:, 2);
    Q(:, 2) = Lam(:, 1) .* U(:, 1) .* U(:, 3) + Lam(:, 2) .* U(:, 2) .* U(:, 4);
    Q(:, 4) = Lam(:, 1) .* U(:, 3) .* U(:, 3) + Lam(:, 2) .* U(:, 4) .* U(:, 4);
    Q(:, 3) = Q(:, 2);
    
    q1 = U(:, 1) .* q(:, 1) + U(:, 2) .* q(:, 2);
    q2 = U(:, 3) .* q(:, 1) + U(:, 4) .* q(:, 2);
    q = cat(2, q1, q2);
end

function [L, U] = Eig_de(Lam)
    [M,~] = size(Lam);
    Det = Lam(:, 1) .* Lam(:, 4) - Lam(:, 2).^2;
    Trace = Lam(:, 1) + Lam(:, 4);
    
    D = sqrt(abs(0.25 * Trace.^2 - Det));
    Lmax = 0.5 * Trace + D;
    Lmin = 0.5 * Trace - D;
    
    Smax = zeros(M, 1); Smin = zeros(M, 1);
    U11 = zeros(M, 1); U12 = zeros(M, 1);
    U21 = zeros(M, 1); U22 = zeros(M, 1);
    
    for i = 1 : M
        Smax_i = Lmax(i); Smin_i = Lmin(i);
        
        S11 = Lam(i, 1); S12 = Lam(i, 2); S22 = Lam(i, 4);
        if S12 == 0
            if S11 >= S22
                v11 = 1; v21 = 0; v12 = 0; v22 = 1;
            else
                v11 = 0; v21 = 1; v12 = 1; v22 = 0;
            end
        else
            v11 = Smax_i - S22; v21 = S12;
            L1 = sqrt(v11^2 + v21^2);
            v11 = v11/L1; v21 = v21/L1;
                
            v12 = Smin_i - S22; v22 = S12;
            L2 = sqrt(v12^2 + v22^2);
            v12 = v12/L2; v22 = v22/L2;
        end    
        Smax(i) = Smax_i; Smin(i) = Smin_i;
        U11(i) = v11; U12(i) = v12;
        U21(i) = v21; U22(i) = v22;
    end
    L = cat(2, Smax, Smin);
    U = cat(2, U11, U12, U21, U22);
end

function [phi, psi] = prox_H(phi, psi)
    spsi = psi.^2;
    ind = phi + spsi/2 > 0;
    itermax = 100;
    
    lam = 1e3 * ones(size(phi));
    lam = lam(ind);
    if ~isempty(ind)
%         length(lam)
        for ell = 1 : itermax
            f = (1 + lam).^2 .* (phi(ind) - lam) + spsi(ind)/2;
%             max(f)
            if max(abs(f)) < 1e-10
                break
            end
            df = 2 * (1 + lam) .* (phi(ind) - lam) - (1 + lam).^2;
            grad = f./df;
            lam = lam - grad;
        end
    end
    phi(ind) = phi(ind) - lam;
    psi(ind) = psi(ind)./(1 + lam);
end
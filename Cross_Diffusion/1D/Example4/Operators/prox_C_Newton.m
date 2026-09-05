function [Q, q] = prox_C_Newton(Q0, q0)
    Itermax = 100; Res = 1e-8;
    Q = Q0; q = q0; [M, N] = size(q0);
    for i = 1 : M
        Q0_i = [Q0(i, 1), Q0(i, 2); Q0(i, 3), Q0(i, 4)];
        q_i = q(i, :)'; q0_i = q0(i, :)';
        for ell = 1 : Itermax
            Y = Q0_i + q_i * q_i'/2; [U, Lam] = eig_de(Y);
            Y_plus = U * max(0, Lam) * U';
            
            grad_f = q_i - q0_i + Y_plus * q_i;
            if norm(grad_f) < Res
                break
            end
            
            A = zeros(2);
            for i1 = 1 : 2
                for i2 = 1 : 2
                    if abs(Lam(i1) - Lam(i2)) < 1e-12;
                        A(i1, i2) = (Lam(i1) > 0);
                    else
                        A(i1, i2) = (max(0, Lam(i1)) - max(0, Lam(i2)))/(Lam(i1) - Lam(i2));
                    end
                end
            end
            Z1 = [2 * q_i(1), q_i(2); q_i(2), 0]/2;
            Z2 = [0, q_i(1); q_i(1), 2 * q_i(2)]/2;
            D1 = U * (A .* (U' * Z1 * U)) * U' * q_i;
            D2 = U * (A .* (U' * Z2 * U)) * U' * q_i;
            
            H = eye(2) + Y_plus + cat(2, D1, D2);
            s = - H\grad_f;
            q_i = q_i + s;
        end
%         ell
%         norm(grad_f)
        Y = Q0_i + q_i * q_i'/2; [U, Lam] = eig_de(Y);
        Y_plus = U * min(0, Lam) * U';
        Q_i = Y_plus - q_i * q_i'/2;
        
        q(i, :) = q_i';
        Q(i, 1) = Q_i(1, 1); Q(i, 2) = Q_i(1, 2);
        Q(i, 3) = Q_i(2, 1); Q(i, 4) = Q_i(2, 2);
    end
end

function [U, L] = eig_de(Lam)
    Det = Lam(1, 1) .* Lam(2, 2) - Lam(1, 2).^2;
    Trace = Lam(1, 1) + Lam(2, 2);
    
    D = sqrt(abs(0.25 * Trace.^2 - Det));
    Lmax = 0.5 * Trace + D;
    Lmin = 0.5 * Trace - D;
    
    S11 = Lam(1, 1); S12 = Lam(1, 2); S22 = Lam(2, 2);
    if S12 == 0
        if S11 >= S22
            v11 = 1; v21 = 0; v12 = 0; v22 = 1;
        else
            v11 = 0; v21 = 1; v12 = 1; v22 = 0;
        end
    else
        v11 = Lmax - S22; v21 = S12;
        L1 = sqrt(v11^2 + v21^2);
        v11 = v11/L1; v21 = v21/L1;
        
        v12 = Lmin - S22; v22 = S12;
        L2 = sqrt(v12^2 + v22^2);
        v12 = v12/L2; v22 = v22/L2;
    end
    L = [Lmax, 0; 0, Lmin];
    U = [v11, v12; v21, v22];
end
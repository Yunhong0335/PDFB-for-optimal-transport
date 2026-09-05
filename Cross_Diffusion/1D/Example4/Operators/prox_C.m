function [Q, q] = prox_C(Q0, q0)
    gam = 1; itermax = 1000;
    Lam = zeros(size(Q0));
    Q = Q0; q = q0;
    for ell = 1 : itermax
        RLam = Lam;
        
        q11 = q(:, 1) .* q(:, 1); q12 = q(:, 2) .* q(:, 1);
        q22 = q(:, 2) .* q(:, 2);
        Res = Q + cat(2, q11, q12, q12, q22)/2;
        
        Lam = Lam + gam * Res;
        Lam = proj_SP(Lam);
        
        if sqrt(sum(sum((RLam - Lam).^2)))/sqrt(sum(sum(RLam.^2))) < 1e-10
%             ell
            break
        end
        
        det = (1 + Lam(:, 1)) .* (1 + Lam(:, 4)) - Lam(:, 2).^2;
        q1 = (1 + Lam(:, 4)) .* q0(:, 1) - Lam(:, 2) .* q0(:, 2);
        q2 = (1 + Lam(:, 1)) .* q0(:, 2) - Lam(:, 3) .* q0(:, 1);
        
        Q = Q0 - Lam; q = cat(2, q1./det, q2./det);
    end
%     max(abs(Q(:, 1)))
end

function Lam = proj_SP(Lam)
    [M,~] = size(Lam);
    Det = Lam(:, 1) .* Lam(:, 4) - Lam(:, 2).^2;
    Trace = Lam(:, 1) + Lam(:, 4);
    
    D = sqrt(abs(0.25 * Trace.^2 - Det));
    Lmax = 0.5 * Trace + D;
    Lmin = 0.5 * Trace - D;
    
    for i = 1 : M
        Smax_i = Lmax(i); Smin_i = Lmin(i);
        if Smin_i < 0
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
            
            Smax_i = max(Smax_i, 0); Smin_i = max(Smin_i, 0);
            
            Lam(i, 1) = Smax_i * v11 * v11 + Smin_i * v12 * v12;
            Lam(i, 2) = Smax_i * v11 * v21 + Smin_i * v12 * v22;
            Lam(i, 4) = Smax_i * v21 * v21 + Smin_i * v22 * v22;
        end
    end
    
    Lam(:, 3) = Lam(:, 2);
end
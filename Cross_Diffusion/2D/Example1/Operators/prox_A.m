function [p, m1, m2] = prox_A(p0, m1, m2, pn, U, Kp, Pcond)
    c = 1; itermax = 10;
    
    M1 = length(m1); M2 = length(m2);
    
    m0 = cat(1, m1, m2);
    Mp = length(p0);
    Mm = length(m0);
    
    x0 = cat(1, p0, m0);
    y0 = pn;

    y = (Kp * Kp')\(Kp * x0 - y0);
%     [y, ~,error, r,~] = pcg(Kp * Kp', (Kp * x0 - y0), 1e-6, 200, Pcond, Pcond');
%     r
    
    x = x0 - Kp' * y;
    p = x(1 : Mp);
    m = x(Mp + 1 : Mp + Mm);
    
    lam = zeros(Mp, 1); Ninx = 0; inx = [];
    for ell = 1 : itermax
        RNinx = Ninx;
        Rinx = inx;
        
        ninx = find(lam + c * (p - U(1)) < 0); Nn = length(ninx);
        pinx = find(lam + c * (p - U(2)) > 0); Np = length(pinx);
        if Nn + Np > 0
            
            inx = cat(1, ninx, pinx);
            Ninx = length(inx);
            PA = sparse(1 : Ninx, inx, ones(Ninx, 1), Ninx, Mp);
            
            if Nn > 0 && Np > 0
                psi = cat(1, U(1) * ones(Nn, 1), U(2) * ones(Np, 1));
            elseif Nn > 0
                psi = U(1) * ones(Nn, 1);
            else
                psi = U(2) * ones(Np, 1);
            end
            
            G = cat(2, PA, sparse(Ninx, Mm));
            G = cat(1, Kp, G);
            
            y0 = cat(1, pn, psi);
            y = (G * G')\(G * x0 - y0);
            
%             Null = sparse(Mp, Ninx);
%             I = sparse(1 : Ninx, 1 : Ninx, ones(Ninx, 1));
%             P = cat(1, cat(2, Pcond, Null), cat(2, Null', I));
%             [y, ~, error, r,~] = pcg(G * G', (G * x0 - y0), 1e-9, 300, P, P');
% %             r
% %             error
            
            x = x0 - G' * y;
            p = x(1 : Mp);
            m = x(Mp + 1 : Mp + Mm);
            lam = PA' * y(Mp + 1 : end);
        end
        if Ninx == RNinx
            if inx == Rinx
%                 ell
                break
            end
        end
    end
    m1 = m(1 : M1);
    m2 = m(M1 + 1 : M1 + M2);
end
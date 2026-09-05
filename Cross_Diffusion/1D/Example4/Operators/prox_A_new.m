function [p, m] = prox_A_new(p0, m0, pn, K)
    c = 1; itermax = 30;
    [Mp,~] = size(p0); [Mm,~] = size(m0);
    pn = cat(1, pn(:, 1), pn(:, 2));
    
    x0 = cat(1, p0(:, 1), p0(:, 2), m0(:, 1), m0(:, 2));
    y0 = pn;
    
    y = (K * K')\(K * x0 - y0);
    x = x0 - K' * y;
    
    p1 = x(1 : Mp); p2 = x(Mp + 1 : 2 * Mp);    
    m1 = x(2 * Mp + 1 : 2 * Mp + Mm);
    m2 = x(end - Mm + 1 : end);
    
    lam0 = zeros(Mp, 1); Ninx0 = 0; inx0 = [];
    lam1 = zeros(Mp, 1); Ninx1 = 0; inx1 = [];
    lam2 = zeros(Mp, 1); Ninx2 = 0; inx2 = [];
    for ell = 1 : itermax
        RNinx0 = Ninx0; Rinx0 = inx0;
        inx0 = find(lam0 + c * (p1 - 0) < 0);
        Ninx0 = length(inx0);
        
        RNinx1 = Ninx1; Rinx1 = inx1;
        inx1 = find(lam1 + c * (p2 - 0) < 0);
        Ninx1 = length(inx1);
        
        RNinx2 = Ninx2; Rinx2 = inx2;
        inx2 = find(lam2 + c * (p1 + p2 - 1) > 0);
        Ninx2 = length(inx2);
        
        Flag0 = 0; Flag1 = 0; Flag2 = 0;
        if Ninx0 > 0
            PA0 = sparse(1 : Ninx0, inx0, ones(Ninx0, 1), Ninx0, Mp);
            G1 = cat(2, PA0, sparse(Ninx0, Mp), sparse(Ninx0, 2 * Mm));
            Flag0 = 1;
        end
        if Ninx1 > 0
            PA1 = sparse(1 : Ninx1, inx1, ones(Ninx1, 1), Ninx1, Mp);
            G2 = cat(2, sparse(Ninx1, Mp), PA1, sparse(Ninx1, 2 * Mm));
            Flag1 = 1;
        end
        if Ninx2 > 0
            PA2 = sparse(1 : Ninx2, inx2, ones(Ninx2, 1), Ninx2, Mp);
            G3 = cat(2, PA2, PA2, sparse(Ninx2, 2 * Mm));
            Flag2 = 1;
        end
        
        G = K;
        psi = zeros(0, 1);
        if Flag0 == 1 && Flag1 == 1 && Flag2 == 1
            psi0 = zeros(Ninx0, 1);
            psi1 = zeros(Ninx1, 1);
            psi2 = ones(Ninx2, 1);
            psi = cat(1, psi0, psi1, psi2);
            G = cat(1, K, G1, G2, G3);
        
        elseif Flag0 == 1 && Flag1 == 1 && Flag2 == 0
            psi0 = zeros(Ninx0, 1);
            psi1 = zeros(Ninx1, 1);
            psi = cat(1, psi0, psi1);
            G = cat(1, K, G1, G2);
        elseif Flag0 == 1 && Flag1 == 0 && Flag2 == 1
            psi0 = zeros(Ninx0, 1);
            psi2 = ones(Ninx2, 1);
            psi = cat(1, psi0, psi2);
            G = cat(1, K, G1, G3);
        elseif Flag0 == 0 && Flag1 == 1 && Flag2 == 1
            psi1 = zeros(Ninx1, 1);
            psi2 = ones(Ninx2, 1);
            psi = cat(1, psi1, psi2);
            G = cat(1, K, G2, G3);
        
        elseif Flag0 == 1 && Flag1 == 0 && Flag2 == 0
            psi = zeros(Ninx0, 1);
            G = cat(1, K, G1);
        elseif Flag0 == 0 && Flag1 == 1 && Flag2 == 0
            psi = zeros(Ninx1, 1);
            G = cat(1, K, G2);
        elseif Flag0 == 0 && Flag1 == 0 && Flag2 == 1
            psi = ones(Ninx2, 1);
            G = cat(1, K, G3);
        end
        y0 = cat(1, pn, psi);
        y = (G * G')\(G * x0 - y0);
        x = x0 - G' * y;
        
        p1 = x(1 : Mp);
        p2 = x(Mp + 1 : 2 * Mp);
            
        m1 = x(2 * Mp + 1 : 2 * Mp + Mm);
        m2 = x(end - Mm + 1 : end);
        
        if Flag0 == 1
            lam0 = PA0' * y(2 * Mp + 1 : 2 * Mp + Ninx0);
        end
        if Flag1 == 1
            lam1 = PA1' * y(2 * Mp + 1 + Ninx0 : 2 * Mp + Ninx0 + Ninx1);
        end
        if Flag2 == 1
            lam2 = PA2' * y(end - Ninx2 + 1 : end);
        end
            
        % convergence
        flag0 = 0; flag1 = 0; flag2 = 0;
        if Flag0 == 1
            if Ninx0 == RNinx0
                if inx0 == Rinx0
                    flag0 = 1;
                end
            end
        elseif Flag0 == 1
            flag0 = 1;
        end
        if Flag1 == 1
            if Ninx1 == RNinx1
                if inx1 == Rinx1
                    flag1 = 1;
                end
            end
        elseif Flag1 == 0
            flag1 = 1;
        end
        if Flag2 == 1
            if Ninx2 == RNinx2
                if inx2 == Rinx2
                    flag2 = 1;
                end
            end
        elseif Flag2 == 0
            flag2 = 1;
        end
        if (flag0 * flag1 * flag2) == 1
%             ell
            break
        end
    end
    p = cat(2, p1, p2); m = cat(2, m1, m2);
end
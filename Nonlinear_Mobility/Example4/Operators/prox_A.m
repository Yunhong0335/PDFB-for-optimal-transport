function [p, ca] = prox_A(p0, m0, U, Kp, Bes)
    c = 1; itermax = 30;
    
    M = length(p0);
    ca = Bes\(p0 - m0);
    p = p0 - ca;
    
    lam = zeros(M, 1);
    for ell = 1 : itermax
        Rlam = lam;
        ninx = find(lam + c * (p - U(1)) < 0); Nn = length(ninx);
        pinx = find(lam + c * (p - U(2)) > 0); Np = length(pinx);
        if Nn + Np > 0
            inx = cat(1, ninx, pinx);
            Ninx = length(inx);
            A = sparse(1 : Ninx, inx, ones(Ninx, 1), Ninx, M);
            if Nn > 0 && Np > 0
                psi = cat(1, U(1) * ones(Nn, 1), U(2) * ones(Np, 1));
            elseif Nn > 0
                psi = U(1) * ones(Nn, 1);
            else
                psi = U(2) * ones(Np, 1);
            end
            EA = cat(2, A, sparse(Ninx, M));
            
            K = cat(1, cat(2, Kp, EA'), cat(2, EA, sparse(Ninx, Ninx)));

            sol = K\cat(1, p0, m0, psi);
            p = sol(1 : M);
            ca = sol(M + 1 : 2 * M);
            lam = A' * sol(2 * M + 1 : end);
        end
        if lam == Rlam
%             ell
            break
        end
    end
end
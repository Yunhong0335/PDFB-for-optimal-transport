function E = grad_e(p, Hes)
    E = (0.018)^2 * Hes * p + p .* (p.^2 - 1);
end
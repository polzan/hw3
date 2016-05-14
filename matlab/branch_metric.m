function bm = branch_metric(sigma_j, sigma_i, rho_k, f)
    if any(sigma_j(2:length(sigma_j)) ~= sigma_i(1:length(sigma_i)-1))
        bm = Inf;
    else
        bm = abs(rho_k - f(sigma_j, sigma_i))^2;
    end
end

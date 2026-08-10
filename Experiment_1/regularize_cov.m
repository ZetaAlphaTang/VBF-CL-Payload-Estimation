function P = regularize_cov(P, eps_val)
% Make covariance matrix symmetric positive definite approximately

P = (P + P') / 2;

min_eig = min(real(eig(P)));

if min_eig < eps_val
    P = P + (eps_val - min_eig) * eye(size(P));
end

P = (P + P') / 2;
end
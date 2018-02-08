%calculates mean of (Z(i)|X(i))
function cond_mu = cond_mean(x, A, cov, mu)
    cond_mu = A'*inv(A*A' + cov)*(x - mu);
end

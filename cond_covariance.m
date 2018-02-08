%calculates covariance of (Z(i)|X(i))
function cond_cov = cond_covariance(k, A, cov)

    cond_cov = eye(k) - A'*inv(A*A' + cov)*A;
end

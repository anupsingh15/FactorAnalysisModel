function l = likelihood1(x,mu,cov,A)
l =0; 
for i =1:size(x,2)
    l = l+log(((2*pi)^(size(x,1)/2)*det(A*A'+cov)^0.5)^-1 * exp(-0.5*(x(:,i)-mu)'*inv(A*A'+cov)*(x(:,i)-mu)));
 end
end
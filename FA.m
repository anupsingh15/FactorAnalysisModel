function [mu covv A L] = FA(X)

t= 0; % to count no. of iterations before convergence
[n m] = size(X); %dimensions of data X matrix = no.of dimensions of each training example * no.of training examples
k = input('Enter the no. of factors you want: '); 
L = [];
%Initialization
Z = rand([k m]);
A = rand([n k]); %factor loading matrix
covv = cov(X');
%update mean vector i.e mu (update has to be done just once)
mu = mean(X')';
l0 = likelihood1(X, mu, covv, A);
l1 = 10^5;

%EM iteration
while(abs(l1-l0) >= 10^-4)

    % E-STEP
    % in this step it is just defining the function cond_mean and
    % cond_covariance which will be called in M step 
    
    % M-STEP
    a = 0;
    %update covariance matrix i.e cov
    for i = 1:m
        phi = a + X(:,i)*X(:,i)' - X(:,i)*cond_mean(X(:,i), A, covv, mu)'*A' - A*cond_mean(X(:,i), A, covv, mu)*X(:,i)' + A*((cond_mean(X(:,i), A, covv, mu))*(cond_mean(X(:,i), A, covv, mu))' + cond_covariance(k, A, covv))*A';
    end
    
    a = 0;
    b = 0;
    %update A
    for i = 1:m
        a = a + (X(:,i)-mu)*(cond_mean(X(:,i), A, covv, mu))';
        b = b + (cond_mean(X(:,i), A, covv, mu))*(cond_mean(X(:,i), A, covv, mu))' + cond_covariance(k, A, covv);
    end
    A = a*inv(b);
     
    %updated diagonal covariance matrix
    covv =diag(diag(phi))/m;
    
    l0 = l1;
    %likelihood value at each iteration
    l1 = likelihood1(X, mu, covv, A);
    t = t+1;
    L = [L;[t l1]]; %storing value of likelihood at each iteration
end

%e =mvnrnd(zeros(n,1),covv,m);
%x =mu +A*Z(:,11) + e(11,:)'
end

        
    
    
    
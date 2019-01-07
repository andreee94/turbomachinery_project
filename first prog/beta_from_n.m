function [ beta, u] = beta_from_n( n, Dm, Tin, cp, k, eta)
 
    if(nargin < 6)
        eta = 1;
    end    

    u = 2 * pi * n / 60 * Dm / 2;
    
  
    f = @(beta) -u^2 + cp * Tin * (beta^-k - 1) * eta;
    beta = secants(f, 0.7, 0.9);
    beta = 1/beta;
end


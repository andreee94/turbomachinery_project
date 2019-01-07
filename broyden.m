function [xv,it] = broyden( f, x0, toll, iter)
    
    if(nargin < 4)
        iter = 1000;
    end
    if(nargin < 3)
        toll = 1e-6;
    end
    
    % Broyden's method for solving a system of n non-linear equations
    % in n variables.
    %
    % Example call: [xv,it] = broyden(x,f,n,tol)
    % Requires an initial approximation column vector x. tol is required
    % accuracy. User must define function f, for example see page 115.
    % xv is the solution vector, parameter it is number of iterations
    % taken. WARNING. Method may fail, for example, if initial estimates
    % are poor.
    %
    n = length(x0);
    fr = zeros(n,1); 
    it = 0;
    xv = x0;
    %Set initial Br
    Br = eye(n);
    fr = feval(f, xv);
    while norm(fr) > toll && it <= iter
        it = it + 1;
        pr = -Br * fr;
        tau = 1;
        xv1 = xv + tau * pr;
        xv = xv1;
        oldfr = fr; 
        fr = feval(f,xv);
        %Update approximation to Jacobian using Broydens formula
        y = fr-oldfr;
        oldBr = Br;
        oyp = oldBr * y-pr;
        pB = pr' * oldBr;
        for i = 1:n
            for j = 1:n
                M(i,j) = oyp(i) * pB(j);
            end;
        end;
        Br = oldBr-M./(pr' * oldBr * y);
    end;
    
    xv = xv(1, :);
    
    
    if it >= iter
       disp('Broyden algorithm do not converge in all iterations') 
    end
    
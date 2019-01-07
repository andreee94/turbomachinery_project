function [ x ] = secanti( f,x0,x1,toll,iter)
    
    if(nargin < 5)
        iter = 1000;
    end
    if(nargin < 4)
        toll = 1e-6;
    end
    
    xold = x0;
    fold = f(xold);
    x = x1;
    
    d = toll+1;
    
    for ii = 1:iter
        if abs(d) < toll
            break;
        end
        fx = f(x);
        d = -(x-xold)/(fx-fold)*fx;
        fold = fx;
        xold = x;
        x = x + d;
    end
    
    
end


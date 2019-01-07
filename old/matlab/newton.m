function [ x ] = newton( f,df,x0,toll,iter)
    
    if(nargin < 5)
        iter = 100;
    end
    if(nargin < 4)
        toll = 1e-6;
    end
    
    x = x0;
    d = toll+1;
    
    for ii = 1:iter
        if abs(d) < toll
            break;
        end
        d = -f(x)/df(x);
        x = x + d;
    end
    
    
end


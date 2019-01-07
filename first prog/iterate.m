function [res, iter, err] = iterate(gl, f, p0)
    err = inf;
    old = p0;
    %first guess
    p = p0;
    iter = 0;
    while err > gl.toll  && iter <= gl.maxiter
        p = f(p);
        
        if length(p0) > 1
            err = abs(old(1) - p(1));
        else
            err = abs(old - p);
        end
        
        old = p;
        
        iter = iter + 1;
    end
    
    res = p;
end


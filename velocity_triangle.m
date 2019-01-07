function [ res ] = velocity_triangle( u, vA, deltavT, alpha1 )
    
    
    res.u1 = u;
    res.u2 = u;
    res.u = u;
    
    res.v1A = vA;
    res.v1T = res.v1A * tan(alpha1);
    res.v1 = norm([res.v1A, res.v1T]);
    res.alpha1 = alpha1;
    res.alpha1_deg = rad2deg(res.alpha1);
    
    res.w1T = res.v1T - u;
    res.w1A = vA;
    res.w1 = norm([res.w1A, res.w1T]);
    res.beta1 = atan(res.w1T/res.w1A);
    res.beta1_deg = rad2deg(res.beta1);
    
    res.v2T = res.v1T - abs(deltavT);
    res.v2A = vA;
    res.v2 = norm([res.v2A, res.v2T]);
    res.alpha2 = atan(res.v2T/res.v2A);
    res.alpha2_deg = rad2deg(res.alpha2);
    
    res.w2T = res.v2T - u;
    res.w2A = vA;
    res.w2 = norm([res.w2A, res.w2T]);
    res.beta2 = atan(res.w2T/res.w2A);
    res.beta2_deg = rad2deg(res.beta2);
    
    res.wvmax = max([res.v1, res.v2, res.w1, res.w2]);
    res.deltaBeta = res.beta1 - res.beta2;
    res.deltaBetaDeg = rad2deg(res.deltaBeta);
    
end


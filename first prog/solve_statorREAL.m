function [ P ] = solve_statorREAL(gl, iterating, P0, triangle, Yp, beta)%, k1, cp)
    
    P = init_point();
    
    P.hT = P0.hT; % fixed
    P.v = triangle.v1;
    P.h = P.hT - P.v^2 / 2; % fixed
    
    pT1 = @(s) XSteam('p_hs', P.hT / 1000, s);
    p1 = @(s) XSteam('p_hs', P.h / 1000, s);
    
    % Yp  = (pT0 - pT1) / (pT1 - p1);
    
    f = @(s) Yp  - (P0.pT - pT1(s)) / (pT1(s) - p1(s));
    
    P.s = secants(f, P0.s, P0.s + gl.secant_entropy_delta0);
    
    P.pT = XSteam('p_hs', P.hT / 1000, P.s);
    P.p = XSteam('p_hs', P.h / 1000, P.s);
    
    P.rho = XSteam('rho_ph', P.p, P.h / 1000);
    P.mu = XSteam('my_ph', P.p, P.h / 1000);
    
    P.hTrel = P.h + triangle.w1^2 / 2;
    P.pTrel = XSteam('p_hs', P.hTrel / 1000, P.s);
    
    %P.v = triangle.v1;
    
    if iterating
        P.TT = XSteam('T_ph', P.pT, P.hT / 1000) + 273.15;
        P.T = XSteam('T_ph', P.p, P.h / 1000) + 273.15;
        
        P.TTrel = XSteam('T_ph', P.pT, P.hTrel / 1000) + 273.15;
        
        P.a = XSteam('w_pT', P.p, P.T - 273.15);
        P.M = triangle.v1 / P.a;
        P.Mrel = triangle.w1 / P.a;
        
        P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
        
        P.hIS = XSteam('h_ps', P.p, P0.s) * 1000;
    end
end

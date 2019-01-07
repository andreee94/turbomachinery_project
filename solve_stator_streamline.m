function [ rho, s ] = solve_stator_streamline( gl, P0, v1sq, Yp )
    
    P = init_point();
    
    P.hT = P0.hT; % fix
    %P.v = triangle.v1;
    P.h = P.hT - v1sq / 2; % fix
    
    pT1 = @(s) XSteam('p_hs', P.hT / 1000, s);
    p1 = @(s) XSteam('p_hs', P.h / 1000, s);
    
    % Yp  = (pT0 - pT1) / (pT1 - p1);
    
    f = @(s) Yp  - (P0.pT - pT1(s)) / (pT1(s) - p1(s));
    
    P.s = secants(f, P0.s, P0.s + gl.secant_entropy_delta0);
    
    %P.pT = XSteam('p_hs', P.hT / 1000, P.s);
    P.p = XSteam('p_hs', P.h / 1000, P.s);
    
    P.rho = XSteam('rho_ph', P.p, P.h / 1000);
    
    rho = P.rho;
    s = P.s;
    
end


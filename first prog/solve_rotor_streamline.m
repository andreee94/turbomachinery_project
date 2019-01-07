function [ rho, P ] = solve_rotor_streamline( gl, P1, w2sq, Yp )
    
    P = init_point();
    
    P.hTrel = P1.hTrel; % fixed
    P.h = P.hTrel - w2sq / 2; % fixed
    
    pT2rel = @(s) XSteam('p_hs', P.hTrel / 1000, s);
    p2 = @(s) XSteam('p_hs', P.h / 1000, s);
    
    % Yp  = (pT1rel - pT2rel) / (pT2rel - p2);
    f = @(s) Yp  - (P1.pTrel - pT2rel(s)) / (pT2rel(s) - p2(s));
    
    P.s = secants(f, P1.s, P1.s + gl.secant_entropy_delta0);
    
    P.p = XSteam('p_hs', P.h / 1000, P.s);
    
    P.rho = XSteam('rho_ph', P.p, P.h / 1000);
    
    rho = P.rho;
end


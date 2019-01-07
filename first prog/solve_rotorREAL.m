function [ P ] = solve_rotorREAL(gl, iterating, P1, l, triangle, Yp, beta)%, k1, cp)
    
    P = init_point();
    
    P.hTrel = P1.hTrel; % fix
    P.v = triangle.v2;
    P.w = triangle.w2;
    P.h = P.hTrel - triangle.w2^2 / 2; % fix
    
    pT2rel = @(s) XSteam('p_hs', P.hTrel / 1000, s);
    p2 = @(s) XSteam('p_hs', P.h / 1000, s);
    
    % Yp  = (pT1rel - pT2rel) / (pT2rel - p2);
    
    f = @(s) Yp  - (P1.pTrel - pT2rel(s)) / (pT2rel(s) - p2(s));
    
    P.s = secants(f, P1.s, P1.s + gl.secant_entropy_delta0);
    
    P.pTrel = XSteam('p_hs', P.hTrel / 1000, P.s) + 273.15;
    P.hT = P.h + triangle.v2^2 / 2;
    P.pT = XSteam('p_hs', P.hT / 1000, P.s);
    P.p = XSteam('p_hs', P.h / 1000, P.s);
    
    P.TT = XSteam('T_ph', P.pT, P.hT / 1000) + 273.15;
    P.T = XSteam('T_ph', P.p, P.h / 1000) + 273.15;
    
    P.rho = XSteam('rho_ph', P.p, P.h / 1000);
    P.mu = XSteam('my_ph', P.p, P.h / 1000);
    P.T = XSteam('T_ph', P.p, P.h / 1000) + 273.15;
    P.s = XSteam('s_ph', P.p, P.h / 1000);
    P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
    
    P.TTrel = XSteam('T_ph', P.pT, P.hTrel / 1000) + 273.15;
    
    
    %P.TT = P.T + P.v^2 / 2 / P.cp;
    %P.pT = P.p + P.rho * P.v^2 / 2 / 1e5;
    % check also hT2 = hT1 - |l|;
    %disp( P.hT - (P1.hT - abs(l)))
    %P.hTrel = P.h + triangle.w2^2 / 2;
    %P.TTrel = P.T + triangle.w2^2 / 2 / P.cp;
    %P.pTrel =  P.p + P.rho * triangle.w2^2 / 2 / 1e5;
    %P.sT = XSteam('s_pT', P.pT, P.TT - 273.15);
    
    %P.v = triangle.v1;
    P.a = XSteam('w_pT', P.p, P.T - 273.15);
    P.M = triangle.v1 / P.a;
    P.Mrel = triangle.w1 / P.a;
    
    % as gamma we use a mean value between inlet and outlet
    P.hIS = XSteam('h_ps', P.p, P1.s) * 1000;
end
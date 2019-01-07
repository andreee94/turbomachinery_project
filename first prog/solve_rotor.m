function [ P ] = solve_rotor(gl, P1, l, triangle, Yp, beta)%, k1, cp)
    
    P = init_point();
    
    P.hTrel = P1.hTrel; % fix
    P.v = triangle.v2;
    P.w = triangle.w2;
    P.h = P.hTrel - triangle.w2^2 / 2; % fix
    
    % first guess of p
    P.p = P1.p / (sqrt(beta));
    
    err = inf;
    old = gl.pT0;
    ii = 0;
    
    p2 = P.p;
    
    while err > gl.toll
        
        P.rho = XSteam('rho_ph', p2, P.h / 1000);
        
        % definition of pT2rel = @(p2) p2 + P.rho * w^2 / 2;
        % Yp  = (pT1rel - pT2rel) / (pT2rel - p2);
        p2 = P1.pTrel - (1 + Yp) * P.rho * P.w^2 / 2 / 1e5;
        
        err = abs(p2 - old);
        old = p2;
        
        ii = ii + 1;
    end
    
    P.p = p2;
    P.rho = XSteam('rho_ph', P.p, P.h / 1000);
    P.mu = XSteam('my_ph', P.p, P.h / 1000);
    P.T = XSteam('T_ph', P.p, P.h / 1000) + 273.15;
    P.s = XSteam('s_ph', P.p, P.h / 1000);
    P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
    
    P.TT = P.T + P.v^2 / 2 / P.cp;
    P.pT = P.p + P.rho * P.v^2 / 2 / 1e5;
    P.hT = P.h + triangle.v2^2 / 2;
    % check also hT2 = hT1 - |l|;
    disp( P.hT - (P1.hT - abs(l)))
    P.hTrel = P.h + triangle.w2^2 / 2;
    P.TTrel = P.T + triangle.w2^2 / 2 / P.cp;
    P.pTrel =  P.p + P.rho * triangle.w2^2 / 2 / 1e5;
    P.sT = XSteam('s_pT', P.pT, P.TT - 273.15);
  
    %P.v = triangle.v1;
    P.a = XSteam('w_pT', P.p, P.T - 273.15);
    P.M = triangle.v1 / P.a;
    P.Mrel = triangle.w1 / P.a;
    
    % as gamma we use a mean value between inlet and outlet
    P.hIS = XSteam('h_ps', P.p, P1.s) * 1000;
end
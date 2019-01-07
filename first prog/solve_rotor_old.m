function [ P ] = solve_rotor_old(gl, P1, l, triangle, Yp, beta)%, k1, cp)
    
    P = init_point();
    
    %h = @(T) XSteam('h_pT', p0, T - 273.15) * 1000;
    
    %if nargin < 4
    %    cp = @(T) XSteam('Cp_pT', p0, T - 273.15) * 1000;
    %    cp0 = cp(TT0);
    %end
    
    P.hT = P1.hT - abs(l); % exact
    P.hTrel = P1.hTrel; % exact
    P.h = P.hTrel - triangle.w2 ^ 2 / 2; % exact
    % FIRST GUESS INITIALIZZATION
    P.p = P1.p / (sqrt(beta)) ;
    P.T = XSteam('T_ph', P.p, P.h / 1000) + 273.15;
    P.TT = P.T + triangle.v2^2 / 2 / gl.cp0;
    P.TTrel = P.T + triangle.w2^2 / 2 / gl.cp0;
    %P.TTrel = P.hTrel / gl.cp0;
    %P.TT = P.hT / gl.cp0;
    %P.T = P.h / gl.cp0;
    
    
    err = inf;
    %toll = 1e-6;
    old = gl.pT0;
    ii = 0;
    
    %first guess
    p2 = P.p;
    
    while err > gl.toll
        
        P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
        P.cv = XSteam('Cv_pT', P.p, P.T - 273.15) * 1000;
        P.gamma = P.cp / P.cv;
        P.k = 1 - P.gamma^-1;
        
        P.T = XSteam('T_ph', P.p, P.h / 1000) + 273.15;
        %P.T = P.h / gl.cp0;				
        %P.T = (P.hTrel - 0.5 * triangle.w2^2) / P.cp;
        
        P.TTrel = P.T + triangle.w2^2 / 2 / P.cp;
        %P.TTrel = P.hTrel / P.cp;		
        
        % TODO take constant coef out from pT1 function handle
        % TT / T = cpTT / cpT = hT / h
        % convenient since I know exactly the enthalpies
        pT2rel = @(p2) p2 * (P.TTrel / P.T) ^ (P.k^-1);
        %pT2rel = @(p2) p2 * (P.TTrel / P.T) ^ (P.k^-1);
        
        f = @(p2) Yp - (P1.pTrel - pT2rel(p2)) ./ (pT2rel(p2) - p2);
        
        %figure
        %plot(1:1:200, f(1:1:200))
        
        % rough approx of P1
        
        p2 = secants(f, p2 - 0.1, p2);
        
        err = abs(p2 - old);
        
        old = p2;
        ii = ii + 1;
        %f(p1)
    end
    
    
    P.p = p2;
    
    P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
    P.cv = XSteam('Cv_pT', P.p, P.T - 273.15) * 1000;
    P.gamma = P.cp / P.cv;
    P.k = 1 - P.gamma^-1;
    
    P.T = XSteam('T_ph', P.p, P.h / 1000) + 273.15;
    P.TT = P.T + triangle.v2^2 / 2 / P.cp;
    P.TTrel = P.T + triangle.w2^2 / 2 / P.cp;
    %P.T = P.h / P.cp; %(P.hTrel - 0.5 * triangle.w2^2) / P.cp;
    %P.TTrel = P.hTrel / P.cp;
    %P.TT = P.hT / P.cp;
    
    P.pT =  p2 * (P.TT / P.T) ^ (P.k^-1);
    P.pTrel =  p2 * (P.TTrel / P.T) ^ (P.k^-1);
    
    %P.h = XSteam('h_pT', P.p, P.T - 273.15) * 1000;
    P.mu = XSteam('my_pT', P.p, P.T - 273.15);
    P.rho = XSteam('rho_pT', P.p, P.T - 273.15);
    
    P.v = triangle.v2;
    P.w = triangle.w2;
    P.M = triangle.v2 / sqrt(P.gamma * gl.R * P.T);
    P.Mrel = triangle.w2 / sqrt(P.gamma * gl.R * P.T);
    
    %P.s = XSteam('s_pT', P.p, P.T - 273.15);
    P.s = XSteam('s_ph', P.p, P.h/1000);
    P.hIS = XSteam('h_ps', P.p, P1.s) * 1000;
    
    %P.hIS = P.cp * P1.T * (P.p / P1.p) ^ (0.5*P.k + 0.5*P1.k);
    %P.TIS = P.hIS / P.cp;
    
    
    
end
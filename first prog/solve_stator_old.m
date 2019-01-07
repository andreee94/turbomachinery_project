function [ P ] = solve_stator_old(gl, P0, triangle, Yp, beta)%, k1, cp)
    
    P = init_point();
    
    %h = @(T) XSteam('h_pT', p0, T - 273.15) * 1000;
    
    %if nargin < 4
    %    cp = @(T) XSteam('Cp_pT', p0, T - 273.15) * 1000;
    %    cp0 = cp(TT0);
    %end
    
    P.hT = P0.hT; % fix
    P.v = triangle.v1;
    P.h = P.hT - P.v^2 / 2; % fix
    
    % FIRST GUESS INITIALIZZATION
    %P.TT = P.hT / gl.cp0;
    P.p = P0.p / (sqrt(beta));
    P.T = XSteam('T_ph', P.p, P.h / 1000) + 273.15;
    P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
    P.TT = P.T + P.v^2 / 2 / P.cp;
    %P.T = P.h / gl.cp0;
    
    
    err = inf;
    %toll = 1e-6;
    old = gl.pT0;
    ii = 0;
    
    %first guess
    p1 = P.p;
    
    while err > gl.toll
        
        P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
        P.cv = XSteam('Cv_pT', P.p, P.T - 273.15) * 1000;
        P.gamma = P.cp / P.cv;
        P.k = 1 - P.gamma^-1;%(gamma-1)/gamma
        
        %P.T = P.h / P.cp;%(P.hT - 0.5 * triangle.v1^2) / P.cp;
        %P.TT = P.hT / P.cp;
        
        %P.TT = P.T + 0.5 * triangle.v1^2 / P.cp;
        
        % TODO take constant coef out from pT1 function handle
        % TT / T = cpTT / cpT = hT / h
        % convenient since I know exactly the enthalpies
        pT1 = @(p1) p1 * (P.TT / P.T) ^ (1 / P.k);
        %pT1 = @(p1) p1 * (P.TT / P.T) ^ (P.k^-1);
        
        f = @(p1) Yp - (P0.pT - pT1(p1)) / (pT1(p1) - p1);
        
        %figure
        %plot(1:1:200, f(1:1:200))
        
        % rough approx of P1
        
        p1 = secants(f, p1 - 0.1, p1);
        
        err = abs(p1 - old);
        
        old = p1;
        ii = ii + 1;
        %f(p1)
    end
    
    
    P.p = p1;
    P.pT =  p1 * (P.TT / P.T) ^ (P.k^-1);
    
    P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
    P.cv = XSteam('Cv_pT', P.p, P.T - 273.15) * 1000;
    P.gamma = P.cp / P.cv;
    P.k = 1 - P.gamma^-1;
    
    P.T = XSteam('T_ph', P.p, P.h / 1000) + 273.15;
    P.TT = P.T + P.v^2 / 2 / P.cp;
    %P.T = P.h / P.cp;%(P.hT - 0.5 * triangle.v1^2) / P.cp;
    %P.TT = P.hT / P.cp;
    
    %P.T = (P.hT - 0.5 * triangle.v1^2) / P.cp;
    %P.TT = P.T + 0.5 * triangle.v1^2 / P.cp;
    
    %P.h = XSteam('h_pT', P.p, P.T - 273.15) * 1000;
    %P.h_other = P.hT - 0.5 * triangle.v1^2;
    P.mu = XSteam('my_pT', P.p, P.T - 273.15);
    P.rho = XSteam('rho_pT', P.p, P.T - 273.15);
    
    P.hTrel = P.h + triangle.w1^2 / 2;
    P.TTrel = P.T + triangle.w1^2 / 2 / P.cp;
    P.pTrel =  p1 * (P.TTrel / P.T) ^ (P.k^-1);
  
    %P.v = triangle.v1;
    P.M = triangle.v1 / sqrt(P.gamma * gl.R * P.T);
    P.Mrel = triangle.w1 / sqrt(P.gamma * gl.R * P.T);
    
    %P.s = XSteam('s_pT', P.p, P.T - 273.15);
    P.s = XSteam('s_ph', P.p, P.h/1000);
    
    % as gamma we use a mean value between inlet and outlet
    P.hIS = XSteam('h_ps', P.p, P0.s) * 1000;
    %P.cp * P0.T * (P.p / P0.p) ^ (0.5*P.k + 0.5*P0.k);
    %P.TIS = P.hIS / P.cp;
    
    
end


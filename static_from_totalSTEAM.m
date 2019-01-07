function [ P ] = static_from_totalSTEAM(gl, pT0, TT0, v0, iterating)
    
    P = init_point();
    
    % I want to find v0 so that b/Dm = 0.025
    % but m = v0 * rho0(v0) * pi * (b/Dm)^2 * Dm;
    % This equation is implicit in v0
    
    P.TT = TT0;
    P.pT = pT0;
    %first  guess
    P.v = v0;
    
    P.hT = XSteam('h_pT', P.pT, P.TT - 273.15) * 1000;
    %static entropy = total entropy
    P.s = XSteam('s_pT', P.pT, P.TT - 273.15);
    P.h = P.hT - v0^2/2;
    P.p = XSteam('p_hs', P.h / 1000, P.s);
    P.T = XSteam('T_hs', P.h / 1000, P.s) + 273.15;
    
    P.rho = XSteam('rho_ph', P.p, P.h / 1000);
    P.mu = XSteam('my_pT',  P.p, P.T - 273.15);
    %P.s = XSteam('s_ph', P.p, P.h / 1000);
    
    %P.hT = P.h + v0^2 / 2;
    %P.sT = XSteam('s_pT', P.pT, P.TT - 273.15);
    
    P.a = XSteam('w_pT', P.p, P.T - 273.15);
    P.M = v0 / P.a;
    
    P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
    P.cv = XSteam('Cv_pT', P.p, P.T - 273.15) * 1000;
    %P.gamma = P.cp / P.cv;
    %P.k = 1 - P.gamma^-1;
end


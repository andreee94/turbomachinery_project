function [ P ] = pressure_temperature_from_total(gl, pT0, TT0, v0)
    
    P = init_point();
    
    P.TT = TT0;
    P.pT = pT0;
    P.v = v0;
    
    % guess
    P.p = P.pT;
    P.T = P.TT;
    
    p0 = P.p;
    T0 = P.T;
    
    err = inf;
    old1 = p0;
    old2 = T0;
    
    ii=0;
    
    while err > gl.toll
        
        P.cp = XSteam('Cp_pT', p0, T0 - 273.15) * 1000;
        P.rho = XSteam('rho_pT', p0, T0 - 273.15);
        
        T0 = P.TT - v0^2 / 2 / P.cp;
        p0 = P.pT - P.rho * v0^2 / 2 / 1e5;
        
        err = max([abs(old1 - p0), abs(old2 - T0)]);
        old1 = p0;
        old2 = T0;
        ii = ii + 1;
    end
    
    P.p = p0;
    P.T = T0;
    P.rho = XSteam('rho_pT', P.p, P.T - 273.15);
    P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
    P.h = XSteam('h_pT', P.p, P.T - 273.15) * 1000;
    P.mu = XSteam('my_pT',  P.p, P.T - 273.15);
    P.s = XSteam('s_ph', P.p, P.h / 1000);
    
    P.hT = P.h + v0^2 / 2;
    P.sT = XSteam('s_pT', P.pT, P.TT - 273.15);
    
    P.a = XSteam('w_pT', P.p, P.T - 273.15);
    P.M = v0 / P.a;
    
    P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
    P.cv = XSteam('Cv_pT', P.p, P.T - 273.15) * 1000;
    P.gamma = P.cp / P.cv;
    P.k = 1 - P.gamma^-1;
end


function [ P ] = pressure_temperature_from_total_old(gl, pT0, TT0, v0, cp0)
    
    P = init_point();
    
    P.TT = TT0;
    P.pT = pT0;
    P.v = v0;
    
    %toll = 1e-14;
    
    err = inf;
    P.p = P.pT;
    old = P.p;
    
    ii=0;
    
    while err > gl.toll
        
        [P.T, P.h] = static_temperature(P.TT, P.v, P.p, cp0);% it is an approximation since we need also p0 to get cp0
        P.cp = XSteam('Cp_pT', P.p, P.T - 273.15) * 1000;
        P.cv = XSteam('Cv_pT', P.p, P.T - 273.15) * 1000;
        P.gamma = P.cp / P.cv;
        P.k = 1 - P.gamma^-1;
        % pT/p = (TT/T)^(1/K)
        P.p = P.pT * (P.TT / P.T) ^ (1 / P.k);
        err = abs(old - P.p);
        old = P.p;
        ii = ii+1;
        
        if gl.maxiter > ii
            break;
        end
    end
    
    P.rho = XSteam('rho_pT', P.p, P.T - 273.15);
    % speed of sound
    P.a = XSteam('w_pT', P.p, P.T - 273.15);
    P.M = P.v / P.a;
    %P.hT = P.cp * P.T + P.v^2 / 2;%XSteam('h_pT', P.p, P.TT - 273.15) * 1000;
    
    P.h = XSteam('h_pT', P.p, P.T - 273.15) * 1000;
    P.hT = P.h +  P.v^2/2;
    
    %P.hT = P.cp * P.TT;
    %P.h = P.hT - P.v^2/2;
    %P.h = XSteam('h_pT', P.p, P.T - 273.15) * 1000;
    
    P.mu = XSteam('my_pT', P.p, P.T - 273.15);
    %P.htable = XSteam('h_pT', P.p, P.T - 273.15) * 1000;
    
    %P.s = XSteam('s_pT', P.p, P.T - 273.15);
    P.s = XSteam('s_ph', P.p, P.h/1000);
    
    %ii
end


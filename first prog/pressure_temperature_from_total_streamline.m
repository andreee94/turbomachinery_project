function [ rho ] = pressure_temperature_from_total_streamline(gl, pT0, TT0, v0)
    
    P = init_point();
    
    P.TT = TT0;
    P.pT = pT0;
    P.v = v0;
    
    P.hT = XSteam('h_pT', P.pT, P.TT - 273.15) * 1000;
    %static entropy = total entropy
    P.s = XSteam('s_pT', P.pT, P.TT - 273.15);
    P.h = P.hT - v0^2/2;
    P.p = XSteam('p_hs', P.h / 1000, P.s);
    
    P.rho = XSteam('rho_ph', P.p, P.h / 1000);
    
    rho = P.rho;
end
function [P] = init_point()
    
    P.p = [];
    P.pT = [];
    P.pTrel = [];
    
    P.T  = [];
    P.TT  = [];
    P.TTrel = [];
    P.TIS = [];
    
    P.h  = [];
    P.hT = [];
    P.hTrel = [];
    P.hIS = 0;%[];
    
    P.cp = [];
    P.cv = [];
    P.gamma = [];
    P.k = [];
    
    P.rho = [];
    P.mu = [];
    
    P.v = [];
    P.w = [];
    P.a = [];
    P.M = [];
    P.Mrel = [];
    
    P.s = [];
    P.sT = [];
    
    %% non thermodynamical properties
    P.b = [];
    P.c = [];
    P.nBlades = [];
    P.chord_su_b = [];
    
    P.span = [];
end


function [ profile ] = create_profile(camber,T)
    
    %camber = deg2rad(camber);
    m = 0;
    p = 0.4;
    
    f = @(x) camber_da_profilo(x,p)-camber;
    
    m = secanti(f,0.03,0.1);
    
    %camber_da_profilo(m,p)
    %camber
    
    a0 = 0.2969;
    a1 = -0.126;
    a2 = -0.3516;
    a3 = 0.2842;
    a4 = -0.1036;%oppure 0.1036
    
    %metà dello spessore
    yT = @(x) T/0.2*(a0*sqrt(x) + a1*x + a2*x.^2 + a3*x.^3 + a4*x.^4);
    
    yC_finoP = @(x) m/p^2*(2*p*x-x.^2);
    yC_dopoP = @(x) m/(1-p)^2*(1-2*p+2*p*x-x.^2);
    
    yU_f = @(x,yC) yC(x) + yT(x)*cos(camber);
    yL_f = @(x,yC) yC(x) - yT(x)*cos(camber);
    
    %iterare x con beta tra 0 e pi
    beta = linspace(0, pi, 100);
    x = 0.5*(1-cos(beta));
    
    xU = x - yT(x) * sin(camber);
    xL = x + yT(x) * sin(camber);
    
    yU = [yU_f(x(x<p),yC_finoP),...
        yU_f(x(x>=p),yC_dopoP) ];
    
    yL = [yL_f(x(x<p),yC_finoP),...
        yL_f(x(x>=p),yC_dopoP) ];
    
    yC = [yC_finoP(x(x<p)),...
        yC_dopoP(x(x>=p)) ];
    
    profile.x = x;
    profile.xU = xU;
    profile.xL = xL;
    profile.yU = yU;
    profile.yL = yL;
    profile.yC = yC;
    profile.m = m;
    profile.p = p;
    profile.T = T;
    
    
end


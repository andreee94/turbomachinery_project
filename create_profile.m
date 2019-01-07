function [ profile ] = create_profile(camber,T)
    
    %camber = deg2rad(camber);
    m = 0;
    p = 0.4;
    
    %camber = abs(camber);
    
    %fold = @(x) camber_da_profilo(x,p)-camber;
    
    %f = @(x) 4*x^2 + 2 * x / tan(abs(camber)) + p * (1-p);
    
    %m = secants(f, 0.25 / tan(abs(camber)), 0.25 / tan(abs(camber)) + 0.1);
    
    %m = sign(camber) * m;
    
    a = 4 * tan(abs(camber)) ;
    b = 2;
    c = p * (p - 1) * tan(abs(camber));
    
    [ms] = roots([a b c]);
    
    m = ms(ms>=0);
    if length(m) > 1
        m = m(1);
    end
    
    if abs(m) > 10
        k = k;
    end
    
    m = sign(camber) * m;
    
    
    %camber_da_profilo(m,p)
    %camber
    
    %     a0 = 0.2969;
    %     a1 = -0.126;
    %     a2 = -0.3516;
    %     a3 = 0.2842;
    %     a4 = -0.1036;%oppure 0.1036
    
    %metà dello spessore
    %yT = @(x) T/0.2*(a0*sqrt(x) + a1*x + a2*x.^2 + a3*x.^3 + a4*x.^4);
    
    yC_finoP = @(x) m/p^2*(2*p*x-x.^2);
    yC_dopoP = @(x) m/(1-p)^2*(1-2*p+2*p*x-x.^2);
    
    pendenzaC_finoP = @(x) 2*m/p^2*(p-x);
    pendenzaC_dopoP = @(x) 2*m/(1-p)/(1-p)*(p-x);
    
    curvaturaC_finoP = @(x) -2*m/p^2;
    curvaturaC_dopoP = @(x) -2*m/(1-p)/(1-p);
    
    %yU_f = @(x,yC) yC(x) + yT(x)*cos(camber);
    %yL_f = @(x,yC) yC(x) - yT(x)*cos(camber);
    
    %iterare x con beta tra 0 e pi
    beta = linspace(0, pi, 100);
    x = 0.5*(1-cos(beta));
    
    %xUP = x - yT(x) * sin(camber);
    %xDOWN = x + yT(x) * sin(camber);
    %
    %yUP = [yU_f(x(x<p),yC_finoP),...
    %    yU_f(x(x>=p),yC_dopoP) ];
    %
    %     yDOWN = [yL_f(x(x<p),yC_finoP),...
    %         yL_f(x(x>=p),yC_dopoP) ];
    
    yC = [yC_finoP(x(x<p)),...
        yC_dopoP(x(x>=p)) ];
    
    mC = [pendenzaC_finoP(x(x<p)), pendenzaC_dopoP(x(x>=p)) ];
    
    rC = [ (1 + pendenzaC_finoP(x(x<p) .^ 2)).^1.5 ./ abs(curvaturaC_finoP(x(x<p))), (1 + pendenzaC_dopoP(x(x>=p) .^ 2)).^1.5 ./ abs(curvaturaC_dopoP(x(x>=p))) ];
    
    abc1 = [x(end) yC(end);
        x(end - 1) yC(end-1);
        x(end - 2) yC(end - 2)];
    
    
    profile.e = fit_circle_through_3_points(abc1);
    
    %pendenzaC_dopoP(x(end))
    
    %     disp('-------------------------')
    %     disp(rad2deg(camber))
    %     disp(rad2deg(atan(mC(1)) - atan(mC(end))))
    
    profile.x = x;
    %     profile.xUp = xUP;
    %     profile.xDOWN = xDOWN;
    %     profile.yUP = yUP;
    %     profile.yDOWN = yDOWN;
    profile.yC = yC;
    profile.mC = mC;
    profile.rC = rC;
    profile.m = m;
    profile.p = p;
    profile.T = T;
    
%     figure
%     hold on
%     grid on
%     plot(profile.x, profile.yC)
%     plot(profile.x, profile.mC)
%     plot(profile.x, profile.rC)
    
end

function [ camber, camber_deg ] = camber_da_profilo(m, p)
    
    % m = freccia massima percentuale di corda
    % p = posizione freccia massima
    
    camber = atan(2*m./p)+atan(2*m./(1-p));
    
    camber_deg = rad2deg(camber);
    
end





syms m p tan_phi

A = 2*m/p;

B = 2*m/(p+1);

p_m = solve(tan_phi == (A+B)/(1-A*B),m);

p_m_func = matlabFunction(p_m(1));

%plot(linspace(0,1,100), p_m_func(linspace(0,1,100), tan(deltabeta_geo)));

m = linspace(0,1,100);
p = p_m_func(m, tan(deltabeta_geo));

plot(camber_da_profilo(m,p))
%plot(atan(2*m./p)+ atan(2*m./(p+1)));
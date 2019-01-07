
global m_P beta_TOT R Cp gamma p0 T0 ws Ds eta l graph v_M my_warning

% indice politropica
n = @(T_OUT, T_IN, beta) log(beta)/(log(beta)-log(T_OUT/T_IN));

%isobara = @(s,rho) T0 * exp(s/Cp)-R*ln(rho);
isobaraTS = @(T,p) (Cp-R) * log(T) + R * log(R*T/p);
adiabaticaTS = @(T,p_IN,rho_IN,n) (Cp-R) * log(T) - R/(n-1) * log(T*R*rho_IN/p_IN) - R*log(rho_IN);

adiabaticaPV = @(p,p_IN,rho_IN,n) (p/p_IN).^(1/n)*rho_IN;

%piano PV
figure(ff); ff = ff+1;
subplot(1,2,1)
for stadio = ST
    array_p = linspace(stadio.p_IN,stadio.p_OUT,100);
    politropica = n(stadio.T_OUT,stadio.T_IN,stadio.beta);
    
    plot(1./adiabaticaPV(array_p,stadio.p_IN,stadio.rho_IN,politropica),array_p,'k-')
    hold on
    plot(1./adiabaticaPV(array_p,stadio.p_IN,stadio.rho_IN,gamma),array_p,'b--')
    plot(1/stadio.rho_IN, stadio.p_IN,'ro')
    plot(1/stadio.rho_OUT, stadio.p_OUT,'ro')
end

array_p = linspace(ST(1).p_IN,ST(end).p_OUT,100);
politropica = n(ST(end).T_OUT, ST(1).T_IN, ST(end).p_OUT/ST(1).p_IN );
plot(1./adiabaticaPV(array_p,ST(1).p_IN,ST(1).rho_IN,politropica),array_p,'m-')

xlabel('Volume specifico')
ylabel('Pressione')


%piano TS
%figure(ff); ff = ff+1;
subplot(1,2,2)
for stadio = ST
    
    politropica = n(stadio.T_OUT,stadio.T_IN,stadio.beta);
    array_T = linspace(stadio.T_IN,stadio.T_OUT,100);
    
    plot(isobaraTS(array_T,stadio.p_IN),array_T,'k--')
    hold on
    plot(isobaraTS(array_T,stadio.p_OUT),array_T,'k--')
    plot(adiabaticaTS(array_T,stadio.p_IN,stadio.rho_IN,politropica),array_T,'k-')
    plot((Cp-R)*log(stadio.T_IN)-R*log(stadio.rho_IN), stadio.T_IN,'ro')
    plot((Cp-R)*log(stadio.T_OUT)-R*log(stadio.rho_OUT), stadio.T_OUT,'ro')
end

array_T = linspace(ST(1).T_IN, ST(end).T_OUT, 100);
politropica = n(ST(end).T_OUT, ST(1).T_IN, ST(end).p_OUT/ST(1).p_IN );
plot(adiabaticaTS(array_T,ST(1).p_IN,ST(1).rho_IN,politropica),array_T,'m-')

xlabel('Entropia specifica')
ylabel('Temperatura')

clear n
clear isobaraTS adiabaticaTS adiabaticaPV
clear politropica array_p array_T









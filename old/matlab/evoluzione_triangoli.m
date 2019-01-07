
T_mid = triangolo_velocita( beta1,beta2,v_M, ST(1).U);

ii = 1;
ST_ii = ST(ii);
d_array = linspace(ST_ii.D * (1 - ST_ii.bsuD), ST_ii.D * (1 + ST_ii.bsuD), 10);

deltabeta_evoluzione = zeros(length(d_array), 1);
beta1_evoluzione = zeros(length(d_array), 1);
beta2_evoluzione = zeros(length(d_array), 1);
l_evoluzione = zeros(length(d_array), 1);

zz = 0;
for d = d_array
    
    zz = zz + 1;
    
    U = ST_ii.w * d/2;
    delta_vT = ST_ii.l/U/v_M;
    
    %vortice libero
    v1T = T_mid.v1T*ST_ii.D/d;
    v2T = T_mid.v2T*ST_ii.D/d;
    
    %alpha1 = atan(v1T/v_M);
    %alpha2 = atan(v2T/v_M);
    
    beta1_evoluzione(zz) =  -atan(v2T/v_M);
    beta2_evoluzione(zz) =  -atan(v1T/v_M);
    % beta2 - beta1
    deltabeta_evoluzione(zz) = beta2_evoluzione(zz) - beta1_evoluzione(zz);
    % come verifica che l è costante
    l_evoluzione(zz) = U*(v2T-v1T);
    
end

% figure(ff); ff = ff + 1;
% plot(d_array, l_evoluzione-ST_ii.l)
% legend('Lavoro lungo la pala')

T_hub = triangolo_velocita( beta1_evoluzione(1),beta2_evoluzione(1),v_M, ST(1).U);
%T_mid = triangolo_velocita( beta1,beta2,v_M, ST(1).U);
T_tip = triangolo_velocita( beta1_evoluzione(end),beta2_evoluzione(end),v_M, ST(1).U);


figure(ff); ff = ff + 1;
plot(d_array, rad2deg(deltabeta_evoluzione))
hold on
plot(ST_ii.D, T_mid.deltaBeta_deg, 'o')
plot(ST_ii.D*(1-ST_ii.bsuD), T_hub.deltaBeta_deg, 'o')
plot(ST_ii.D*(1+ST_ii.bsuD), T_tip.deltaBeta_deg, 'o')
grid on
legend('Delta Beta', 'Diametro medio', 'Hub', 'Tip')

xlabel('Diametro')
ylabel('Delta Beta')


figure(ff); ff = ff + 1;
plot(d_array, rad2deg(beta1_evoluzione))
hold on
plot(d_array, rad2deg(beta2_evoluzione))
plot(ST_ii.D*(1-ST_ii.bsuD), rad2deg(T_hub.beta1), 'ro', 'LineWidth',2)
plot(ST_ii.D*(1+ST_ii.bsuD), rad2deg(T_tip.beta1), 'bo', 'LineWidth',2)
plot(ST_ii.D*(1-ST_ii.bsuD), rad2deg(T_hub.beta2), 'ro', 'LineWidth',2)
plot(ST_ii.D*(1+ST_ii.bsuD), rad2deg(T_tip.beta2), 'bo', 'LineWidth',2)
grid on
legend('Beta1', 'Beta2', 'Hub','Tip');%,'Hub', 'Tip')

xlabel('Diametro')
ylabel('Beta')

clear ii ST_ii zz 
%clear d_array
%clear deltabeta_evoluzione beta1_evoluzione beta2_evoluzione
clear U delta_vT V1T V2T






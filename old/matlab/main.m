
warning off

global m_P beta_TOT R Cp gamma p0 T0 ws Ds eta l graph v_M my_warning

graph = true;
my_warning = false;
ff = 1;

howell_init
balje_init

m_P = 65;%kg/s
beta_TOT = 13.5;
R = 287;%K/kgK
gamma = 1.4;
Cp = gamma/(gamma-1)*R;%K/kgK
p0 = 0.85*1e5;%pa
T0 = 268;%k
rho0 = p0/R/T0;
%da balje
% ws = 3.4;%4;
% Ds = 1./ws*2*500/sqrt(412*Cp*(1.265^(0.4/1.4)-1));
%  ws = 3.2;
%  Ds = 1.73;
% ETA OTTIMO
%ws = 2.8737;
%Ds = 2.5;
% BETA MEDIO OTTIMO
ws = 2.3828;
Ds = 2.0333;

%beta che va male
%ws = 4.2;
%Ds = 2.2;
% BETA TOT 13.5
%ws = 2.68072289156627;
%Ds = 1.88433734939759;
% ws = 3;
% Ds = 1.7857;
%ws = 2.8737;
%Ds = 1.8053;
%ws = 2.8895
%Ds = 1.8
eta = balje(ws,Ds);%da balje

l = @(eta, T_in, beta) Cp*T_in/eta*(beta^((gamma-1)/gamma)-1);


%% OTTIMIZZAZIONE BETA CON U MASSIMO POSSIBILE

uu_opt = zeros(500,1); jj = 1;
%beta_approx_array = linspace(1.19,1.21,500);
beta_approx_array = linspace(1.1,1.4,500);

for beta_approx = beta_approx_array;
    
    [ST,~,n_intermedio] = calcolo_stadi_beta_costante(beta_approx);
    
    uu_opt(jj) = ST(n_intermedio).U;
    jj = jj + 1;
    
end

%plot(beta_approx_array,uu_opt);

beta_approx = beta_approx_array(max(find(uu_opt<500)));
%beta_approx = 1.2;
%uu_opt(max(find(uu_opt<500)))

%% CALCOLO STADI CON BETA OTTIMIZZATO

[ST,N,n_intermedio] = calcolo_stadi_beta_costante(beta_approx);

ST_stima_beta_cost = ST;

% ST(n_intermedio)

%calcolo vmeridiana imponendo b/Dm nel primo stadio massimo

%v_m = m_P/(ST(1).rho_IN*pi*ST(1).D*0.4) % 0.4 = b/D max

%b_su_Dm_OUT =

% disp('L tot con formula totale')
% l_TOT = l(eta,T0,beta_TOT)
%
% disp('L tot sommando tutti i lavori')
% sum([ST.l])
%
% disp('L tot lavoro intermedio per N')
% ST(n_intermedio).l * N

jj = 1;
v_m = zeros(500,1);
alpha1_array = linspace(0.1,deg2rad(50),500);

for alpha1 = alpha1_array
    
    beta2 = -alpha1;
    
    %il meno perchè beta2<beta1
    deltaBeta = howell(howell_POLY,beta2);
    %deltaBeta = deg2rad(deltaBeta);
    
    v_m(jj) = ST(n_intermedio).l/ST(n_intermedio).U/(abs(tan(abs(beta2)+deltaBeta))-tan(abs(beta2)));
    %v_m(jj) = ST(n_intermedio).l/ST(n_intermedio).U/(tan(beta2+deltaBeta)-tan(beta2));
    jj = jj +1;
end

% v_m_MIN = m_P / pi / ST(n_intermedio).D^2 / rho0 / 0.4;
%
% T_OUT_IS_stima = T0 * beta_TOT^((gamma-1)/gamma);
%
% T_OUT_stima = T0+(T_OUT_IS_stima-T0)/0.75;
%
% rho_OUT_stima = p0*beta_TOT/R/T_OUT_stima;
%
% v_m_MAX = m_P / pi / ST(n_intermedio).D^2 / rho_OUT_stima / 0.03;

if graph
    figure(ff); ff = ff+1;
    plot(rad2deg(alpha1_array),v_m);
    xlabel('|\beta2|');
    ylabel('V_m');
end

ST0 = ST(1);
ST0.l = ST(n_intermedio).l;
ST0.U = ST(n_intermedio).U;
ST0.D = ST(n_intermedio).D;
ST0.w = ST(n_intermedio).w;

%ST = calcolo_stadi_l_costante(N,n_intermedio,ST(n_intermedio));

v_M = max(v_m);

ricalcola_stadi = false;
ST = calcolo_stadi_l_costante_da_inizio(N+2,ST0);

v_M_min = m_P / pi / ST(n_intermedio).D^2 / ST(1).rho_IN / 0.4;
v_M_max = m_P / pi / ST(n_intermedio).D^2 / ST(end).rho_IN / 0.03;

if v_M > v_M_max
    v_M = v_M_max;
    disp('Si limita v_M alla massima possibile per bsuD')
    disp('Si ricalcolano gli stadi')
    ricalcola_stadi = true;
end

if ricalcola_stadi
    ST = calcolo_stadi_l_costante_da_inizio(N+2,ST0);
end

%beta2_deg = rad2deg(-alpha1_array(v_M == v_m));
beta2_deg = rad2deg(-alpha1_array(abs(v_M-v_m) == min(abs(v_M-v_m))));
beta2 = deg2rad(beta2_deg);

[deltabeta, beta1] = howell(howell_POLY, deg2rad(beta2_deg));

deltabeta_deg = rad2deg(deltabeta);
beta1_deg = beta2_deg-deltabeta_deg;%rad2deg(beta1);
beta1 = deg2rad(beta1_deg);

if graph
    hold on
    plot(-beta2_deg,v_M,'ro');
    grid on
end

plot_point_on_balje


disp('----------------------------------------------------')
disp('INIZIO SIMULAZIONE:')

disp(['T0 = ', num2str(T0)])
disp(['p0 = ', num2str(p0)])
disp(['Portata massica = ', num2str(m_P)])
disp(['Beta richiesto = ', num2str(beta_TOT)])
disp(['ws = ', num2str(ws)])
disp(['Ds = ', num2str(Ds)])
disp(['eta = ', num2str(eta)])

disp(['Beta per stadi a beta costante = ', num2str(beta_approx)])
disp(['U = ', num2str(ST(1).U)])
disp(['Diametro medio = ', num2str(ST(1).D)])
disp(['Velocità di rotazione = ', num2str(ST(1).w)])
disp(['Velocità di rotazione = ', num2str(ST(1).rpm)])

disp('----------------------------------------------------')
disp('Triangoli di velocità:')

disp(['Beta1 = ', num2str(beta1_deg)])
disp(['Beta2 = ', num2str(beta2_deg)])
disp(['Delta beta = ', num2str(deltabeta_deg)])
disp(['Vm = ', num2str(v_M)])
disp(['Vm minima = ', num2str(v_M_min)])
disp(['Vm massima = ', num2str(v_M_max)])

if ST(end).bsuD < 0.03
    disp('v_m è troppo grande');
elseif ST(1).bsuD > 0.4
    disp('v_m è troppo piccolo');
else disp('v_m è OK')
end

disp('----------------------------------------------------')
disp(['Rapporto compressione = ', num2str(ST(end).p_OUT/p0)])

disp(['Rapporto compressione medio = ', num2str((ST(end).p_OUT/p0)^(1/(length(ST))))])

disp(['Temperatura finale = ', num2str(ST(end).T_OUT)])

disp(['Numero stadi = ', num2str(length(ST))])

L_punto = sum([ST.deltaH_IS])*m_P;
L_punto_entrante = length(ST)*ST(1).l*m_P;

disp(['Rendimento totale politropico = ', num2str(L_punto/L_punto_entrante)])

disp(['Rendimento totale isoentropico = ', num2str(l(1,T0,ST(end).p_OUT/p0)*m_P/L_punto_entrante)])

disp('----------------------------------------------------')
disp('Angoli geometrici e cinematici:')

evoluzione_triangoli

plot_triangolo_velocita

deltabeta_geo_init

plot_howell

plot_TS

figure(ff); ff = ff+1;
plot_profile(T_mid.profile,T_mid.beta1,4)

%%
figure(ff); ff = ff+1;
plot_profile(T_hub.profile,T_hub.beta1)
plot_profile(T_mid.profile,T_mid.beta1)
plot_profile(T_tip.profile,T_tip.beta1)

disp('FINE SIMULAZIONE')
disp('----------------------------------------------------')













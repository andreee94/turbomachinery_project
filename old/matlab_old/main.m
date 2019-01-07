
global m_P beta_TOT R Cp gamma p0 T0 ws Ds eta l 


m_P = 65;%kg/s
beta_TOT = 13.5;
R = 287;%K/kgK
gamma = 1.4;
Cp = gamma/(gamma-1)*R;%K/kgK
p0 = 0.85*1e5;%pa
T0 = 268;%k
%da balje
ws = 3;
Ds = 2.5;
eta = 0.8;%da balje

l = @(eta, T_in, beta) Cp*T_in/eta*(beta^((gamma-1)/gamma)-1);


%% OTTIMIZZAZIONE BETA CON U MASSIMO POSSIBILE

uu_opt = zeros(500,1); jj = 1;
beta_approx_array = linspace(1.19,1.21,500);

for beta_approx = beta_approx_array;
    
    ST = calcolo_stadi_beta_costante(beta_approx);
    
    uu_opt(jj) = ST(n_intermedio).U;
    jj = jj + 1;
    
end

%plot(beta_approx_array,uu_opt);

beta_approx = beta_approx_array(max(find(uu_opt<500)));

%% CALCOLO STADI CON BETA OTTIMIZZATO

[ST,N,n_intermedio] = calcolo_stadi_beta_costante(beta_approx);

ST(n_intermedio)

%calcolo vmeridiana imponendo b/Dm nel primo stadio massimo

v_m = m_P/(ST(1).rho_IN*pi*ST(1).D*0.4) % 0.4 = b/D max

%b_su_Dm_OUT = 

disp('L tot con formula totale')
l_TOT = l(eta,T0,beta_TOT)

disp('L tot sommando tutti i lavori')
sum([ST.l])

disp('L tot lavoro intermedio per N')
ST(n_intermedio).l * N

jj = 1;
v_m = zeros(500,1);
alpha1_array = linspace(0,50,500);

for alpha1 = alpha1_array

   beta2 = -alpha1;
   
   deltaBeta = 20;%howell(beta2);
   
   v_m(jj) = ST(n_intermedio).l/ST(n_intermedio).U/(tan(beta2*deltaBeta)-tan(beta2));
    
end


%plot(alpha1_array,v_m);

ST0 = ST(1);
ST0.l = ST(n_intermedio).l;

ST = calcolo_stadi_l_costante(N,n_intermedio,ST(n_intermedio));

plot([ST.ws],[ST.Ds])















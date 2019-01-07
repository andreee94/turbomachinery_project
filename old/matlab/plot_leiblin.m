
beta1_array = linspace(0,65,100);
tsuc_array = linspace(0,0.11,100);

figure(ff); ff = ff+1;

subplot(2,2,1)
plot(beta1_array,polyval(i0_10_POLY,beta1_array));
hold on
plot(abs(T_mid.beta1_deg),polyval(i0_10_POLY,abs(T_mid.beta1_deg)),'o');
grid on

legend('(i_0)_1_0','Condizione operativa')
xlabel('|\beta1|')
ylabel('(i_0)_1_0')


subplot(2,2,2)
plot(tsuc_array,polyval(k_i_t_POLY,tsuc_array));
hold on
plot(T_mid.GEO.t_su_c,polyval(k_i_t_POLY,T_mid.GEO.t_su_c),'o');
grid on

legend('(K_i)_t','Condizione operativa')
xlabel('t/c')
ylabel('(K_i)_t')



subplot(2,2,3)
plot(beta1_array,polyval(delta0_10_POLY,beta1_array));
hold on
plot(abs(T_mid.beta1_deg),polyval(delta0_10_POLY,abs(T_mid.beta1_deg)),'o');
grid on

legend('(\delta_0)_1_0','Condizione operativa')
xlabel('|\beta1|')
ylabel('(\delta_0)_1_0')


subplot(2,2,4)
plot(tsuc_array,polyval(k_delta_t_POLY,tsuc_array));
hold on
plot(T_mid.GEO.t_su_c,polyval(k_delta_t_POLY,T_mid.GEO.t_su_c),'o');
grid on

legend('(K_\delta)_t','Condizione operativa')
xlabel('t/c')
ylabel('(K_\delta)_t')

tightfig


clear beta1_array
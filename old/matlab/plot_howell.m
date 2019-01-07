
beta2_array = linspace(0,50,100);

figure(ff); ff = ff+1;
plot(beta2_array,polyval(howell_POLY,beta2_array));
hold on
plot(abs(beta2_deg),deltabeta_deg,'o');
grid on

legend('Howell','Condizione operativa')
xlabel('|\beta2|')
ylabel('\Delta\beta')

clear beta2_array
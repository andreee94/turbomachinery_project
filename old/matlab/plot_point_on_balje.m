
% if graph
figure(ff); ff = ff+1;
loglog(balje_data_80_X,balje_data_80_Y)%,'*');
hold on
axis equal
loglog(balje_data_75_X,balje_data_75_Y)%,'*');
loglog(balje_data_70_X,balje_data_70_Y)%,'*');
loglog(balje_data_60_X,balje_data_60_Y)%,'*');
loglog(balje_data_50_X,balje_data_50_Y)%,'*');
%loglog(balje_data_axial_X,balje_data_axial_Y)%,'*');
%loglog(balje_data_axial_X_false,balje_data_axial_Y_false)%,'*');

ws_ = [ST.ws];
Ds_ = [ST.Ds];

%ws_ = [iter_res_noempty.ws];
%Ds_ = [iter_res_noempty.Ds];

for ii = 1:length(ws_)
    plot(ws_(ii),Ds_(ii),'ro');
end

ws__ = linspace(1,20,30);
ds__14 = 1./ws__*2*ST(1).U/sqrt(T0*Cp*(1.4^(0.4/1.4)-1));
ds__13 = 1./ws__*2*500/sqrt(T0*Cp*(1.3^(0.4/1.4)-1));
ds__12 = 1./ws__*2*500/sqrt(T0*Cp*(1.2^(0.4/1.4)-1));
ds__11 = 1./ws__*2*500/sqrt(T0*Cp*(1.1^(0.4/1.4)-1));

plot(ws__,ds__14)
%plot(ws__,ds__13)
%plot(ws__,ds__12)
%plot(ws__,ds__11)

legend('80%','75%','70%','60%','50%')
grid on
xlabel('w_s')
ylabel('D_s')

% end

clear ws__ ws_ Ds_
clear ds__14 ds__13 ds__12 ds__11



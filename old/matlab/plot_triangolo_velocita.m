
T_array = [T_hub, T_mid, T_tip];
%T_array = T_mid%REMOVE
legend_text = ['Hub'; 'Mid'; 'Tip'];
%legend_text = ['Mid']

row = length(T_array);
col = 1;
sub_index = 0;
figure(ff); ff = ff+1;

for T = T_array
    
    A = [0,T.v1M];
    
    B1 = [T.v1T,0];
    
    C1 = [T.w1T,0];
    
    B2 = [T.v2T,0];
    
    C2 = [T.w2T,0];
    
    %figure(ff); ff = ff+1;
    
    %% entrata rotore
    sub_index = sub_index + 1;
    subplot(row,col,sub_index); 
    plot([-1.5*abs(C1(1)) 1.5*abs(B1(1))],[A(2) A(2)], 'k--')
    axis equal
    hold on
    plot([-1.5*abs(C1(1)) 1.5*abs(B1(1))],[0 0], 'k--')
    plot([0 0],[-0.5*A(2) 1.5*A(2)], 'k--')
    
    v1 = '\leftarrow v_1';
    w1 = '\leftarrow w_1';
    u1 = '\uparrow u';
    text([A(1)+B1(1)]/2,[A(2)+B1(2)]/2,v1)
    text([C1(1)+A(1)]/2,[C1(2)+A(2)]/2,w1)
    text([B1(1)+C1(1)]/2,[B1(2)+C1(2)]/2,u1)
    
    plot([A(1),B1(1),C1(1),A(1)], [A(2),B1(2),C1(2),A(2)], 'b','LineWidth',2.5)
    
    %% uscita rotore
    %subplot(row,col,sub_index); sub_index = sub_index + 1;
    plot([-1.5*abs(C2(1)) 1.5*abs(B2(1))],[A(2) A(2)], 'k--')
    axis equal
    hold on
    plot([-1.5*abs(C2(1)) 1.5*abs(B2(1))],[0 0], 'k--')
    plot([0 0],[-0.5*A(2) 1.5*A(2)], 'k--')
    
    v1 = '\leftarrow v_2';
    w1 = '\leftarrow w_2';
    u1 = '\uparrow u';
    text([A(1)+B2(1)]/2,[A(2)+B2(2)]/2,v1)
    text([C2(1)+A(1)]/2,[C2(2)+A(2)]/2,w1)
    text([B2(1)+C2(1)]/2,[B2(2)+C2(2)]/2,u1)
    
    plot([A(1),B2(1),C2(1),A(1)], [A(2),B2(2),C2(2),A(2)], 'r','LineWidth',1.5)
    legend(legend_text(sub_index,:))
    set(gca, 'xdir','reverse')
end

%tightfig

clear legend_text 
clear row col sub_index
clear A B1 B2 C1 C2 
clear v1 w1 u1
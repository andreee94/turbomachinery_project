function [ ST,N,n_intermedio ] = calcolo_stadi_l_costante_da_inizio_PAR(N,ST0,DATA, DATA_iter)

%global m_P beta_TOT R Cp gamma p0 T0 ws Ds eta l v_M

%N = ceil(beta_TOT^(1/beta_approx));
%numero stadio intermedio
%n_intermedio = ceil(N/2);

%ST(1).T_IN = T0;
%ST(1).p_IN = p0;

ST(1) = ST0;

for ii = 1:N
    
    if ii >1
        ST(ii).T_IN = ST(ii-1).T_OUT;
        ST(ii).p_IN = ST(ii-1).p_OUT;
        ST(ii).rho_IN = ST(ii-1).rho_OUT;
        ST(ii).eta =  ST(ii-1).eta;%balje(ST0.ws,ST0.Ds);
    else ST(ii).eta = ST0.eta;
    end
    
    ST(ii).Q = DATA.m_P / ST(ii).rho_IN;
    ST(ii).l = ST0.l;
    ST(ii).T_OUT = ST(ii).l/DATA.Cp+ST(ii).T_IN;
    
    ST(ii).w =  ST0.w;
    ST(ii).D =  ST0.D;
    
    ST(ii).eta = 0.60;
    ST(ii).deltaH_IS =  DATA.Cp*(ST(ii).T_OUT-ST(ii).T_IN)*ST(ii).eta;
    ST(ii).ws =  ST(ii).w/ST(ii).deltaH_IS^0.75*ST(ii).Q^0.5;
    ST(ii).Ds =  ST(ii).D*ST(ii).deltaH_IS^0.25/ST(ii).Q^0.5;
    
    counter = 0;
    while true
        
        eta_old = ST(ii).eta;
        ST(ii).eta = balje_PAR(ST(ii).ws,ST(ii).Ds,DATA);
        
        if eta_old == ST(ii).eta || isnan(ST(ii).eta) || counter>5
            break;
        end
        
        ST(ii).deltaH_IS =  DATA.Cp*(ST(ii).T_OUT-ST(ii).T_IN)*ST(ii).eta;
        ST(ii).ws =  ST(ii).w/ST(ii).deltaH_IS^0.75*ST(ii).Q^0.5;
        ST(ii).Ds =  ST(ii).D*ST(ii).deltaH_IS^0.25/ST(ii).Q^0.5;
        
        counter = counter+1;
    end
    
    
    ST(ii).beta = (ST(ii).l*ST(ii).eta/DATA.Cp/ST(ii).T_IN+1)^(DATA.gamma/(DATA.gamma-1));
    
    ST(ii).T_OUT_IS = ST(ii).eta*(ST(ii).T_OUT-ST(ii).T_IN)+ST(ii).T_IN;
    ST(ii).p_OUT = ST(ii).p_IN*ST(ii).beta;
    ST(ii).rho_OUT = ST(ii).p_OUT/ST(ii).T_OUT/DATA.R;
    ST(ii).rpm = ST(ii).w*60/2/pi;
    ST(ii).U = 0.5*ST(ii).D * ST(ii).w;
    ST(ii).bsuD = DATA.m_P/(ST(ii).rho_IN*DATA_iter.v_M*ST(ii).D^2*pi);
    ST(ii).Mach = ST(ii).U / sqrt(DATA.gamma*DATA.R*ST(ii).T_IN);
   
end


end
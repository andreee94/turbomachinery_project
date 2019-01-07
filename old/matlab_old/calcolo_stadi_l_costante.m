function [ ST,N,n_intermedio ] = calcolo_stadi_l_costante(N,n_intermedio,ST0)
    
    global m_P beta_TOT R Cp gamma p0 T0 ws Ds eta l
    
    %N = ceil(beta_TOT^(1/beta_approx));
    %numero stadio intermedio
    %n_intermedio = ceil(N/2);
    
    %ST(1).T_IN = T0;
    %ST(1).p_IN = p0;
    
    ST(n_intermedio) = ST0;
    
    for ii = n_intermedio+1:N
        
        if ii >n_intermedio
            ST(ii).T_IN = ST(ii-1).T_OUT;
            ST(ii).p_IN = ST(ii-1).p_OUT;
            ST(ii).rho_IN = ST(ii-1).rho_OUT;
        end
        
        ST(ii).Q = m_P / ST(ii).rho_IN;
        ST(ii).l = ST0.l;
        ST(ii).T_OUT = ST(ii).l/Cp+ST(ii).T_IN;
        ST(ii).deltaH_IS =  Cp*(ST(ii).T_OUT-ST(ii).T_IN);
        
        ST(ii).w =  ST0.w;
        ST(ii).D =  ST0.D;
        ST(ii).ws =  ST(ii).w/ST(ii).deltaH_IS^0.75*ST(ii).Q^0.5;
        ST(ii).Ds =  ST(ii).D*ST(ii).deltaH_IS^0.25/ST(ii).Q^0.5;
        
        ST(ii).eta = 0.8;%balje(ST(ii).ws,ST(ii).Ds);
        
        ST(ii).beta = (ST(ii).l*ST(ii).eta/Cp/ST(ii).T_IN+1)^(gamma/(gamma-1));
        
        ST(ii).p_OUT = ST(ii).p_IN*ST(ii).beta;
        ST(ii).rho_OUT = ST(ii).p_OUT/ST(ii).T_OUT/R;
        ST(ii).rpm = ST(ii).w*60/2/pi;
        ST(ii).U = 0.5*ST(ii).D * ST(ii).w;
    end
    
    
end
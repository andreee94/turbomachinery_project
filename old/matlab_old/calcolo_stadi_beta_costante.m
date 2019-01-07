function [ ST,N,n_intermedio ] = calcolo_stadi_beta_costante(beta_approx)
    
    global m_P beta_TOT R Cp gamma p0 T0 ws Ds eta l 
    
    N = ceil(beta_TOT^(1/beta_approx));
    %numero stadio intermedio
    n_intermedio = ceil(N/2);
    
    ST(1).T_IN = T0;
    ST(1).p_IN = p0;
    
    for ii = 1:N
        
        if ii>1
            ST(ii).T_IN = ST(ii-1).T_OUT;
            ST(ii).p_IN = ST(ii-1).p_OUT;
        end
        ST(ii).beta = beta_approx;
        ST(ii).eta = eta;
        ST(ii).p_OUT = ST(ii).p_IN*ST(ii).beta;
        ST(ii).l = l(ST(ii).eta, ST(ii).T_IN, ST(ii).beta);
        ST(ii).T_OUT = ST(ii).l/Cp+ST(ii).T_IN;
        ST(ii).rho_IN = ST(ii).p_IN/ST(ii).T_IN/R;
        ST(ii).rho_OUT = ST(ii).p_OUT/ST(ii).T_OUT/R;
        % Q con rho media
        %ST(ii).Q = m_P / (ST(ii).rho_IN + ST(ii).rho_OUT)*2;
        ST(ii).Q = m_P / ST(ii).rho_IN;
        ST(ii).deltaH_IS =  ST(ii).l * ST(ii).eta;
        ST(ii).ws =  ws;
        ST(ii).Ds =  Ds;
        ST(ii).w =  ST(ii).ws*ST(ii).deltaH_IS^0.75/ST(ii).Q^0.5;
        ST(ii).D =  ST(ii).Ds/ST(ii).deltaH_IS^0.25*ST(ii).Q^0.5;
        ST(ii).rpm = ST(ii).w*60/2/pi;
        ST(ii).U = 0.5*ST(ii).D * ST(ii).w;
    end
    
    
end


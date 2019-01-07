function [ ST,N,n_intermedio ] = calcolo_stadi_beta_costante(beta_approx)
    
    global array_of_struct
    global m_P beta_TOT R Cp gamma p0 T0 ws Ds eta l v_M
    
    N = ceil(log(beta_TOT)/log(beta_approx));%ceil(beta_TOT^(1/beta_approx));
    %numero stadio intermedio
    n_intermedio = ceil(N/2);
    
    ST(1).T_IN = T0;
    ST(1).p_IN = p0;
    
    for ii = 1:n_intermedio
        
        if ii>1
            %init struct array
            
            ST(ii).T_IN = ST(ii-1).T_OUT;
            ST(ii).p_IN = ST(ii-1).p_OUT;
        end
        
        ST_ii = ST(ii);
        
        ST_ii.beta = beta_approx;
        ST_ii.eta = eta;
        ST_ii.p_OUT = ST_ii.p_IN*ST_ii.beta;
        ST_ii.l = l(ST_ii.eta, ST_ii.T_IN, ST_ii.beta);
        ST_ii.T_OUT = ST_ii.l/Cp+ST_ii.T_IN;
        ST_ii.rho_IN = ST_ii.p_IN/ST_ii.T_IN/R;
        ST_ii.rho_OUT = ST_ii.p_OUT/ST_ii.T_OUT/R;
        ST_ii.T_OUT_IS = ST_ii.eta*(ST_ii.T_OUT-ST_ii.T_IN)+ST_ii.T_IN;
        % Q con rho media
        %ST_ii.Q = m_P / (ST_ii.rho_IN + ST_ii.rho_OUT)*2;
        ST_ii.Q = m_P / ST_ii.rho_IN;
        ST_ii.deltaH_IS =  ST_ii.l * ST_ii.eta;
        %ST_ii.deltaH_IS2 = Cp*(ST_ii.T_OUT_IS-ST_ii.T_IN);
        ST_ii.ws =  ws;
        ST_ii.Ds =  Ds;
        ST_ii.w =  ST_ii.ws*ST_ii.deltaH_IS^0.75/ST_ii.Q^0.5;
        ST_ii.D =  ST_ii.Ds/ST_ii.deltaH_IS^0.25*ST_ii.Q^0.5;
        ST_ii.rpm = ST_ii.w*60/2/pi;
        ST_ii.U = 0.5*ST_ii.D * ST_ii.w;
        %ST_ii.bsuD = m_P/(ST_ii.rho_IN*v_M*ST_ii.D^2*pi);
        
        if ii==1
            clear ST;
        end
        
        ST(ii) = ST_ii;
        
        if ii==1
            if isempty(array_of_struct)
                array_of_struct = repmat(ST(1), 20, 1 );
            end
%             ST = repmat(ST(1), n_intermedio, 1 );
            ST = array_of_struct;
            ST(ii) = ST_ii;
        end
        
    end
    
    
end


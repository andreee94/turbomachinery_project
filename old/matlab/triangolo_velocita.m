function [ T ] = triangolo_velocita( beta1,beta2,v_M, U)
    
    T.vM = v_M;
    T.U = U;
    %
    %     rad2deg(beta1)
    %     rad2deg(beta2)
    
    T.beta1 = beta1;
    T.beta2 = beta2;
    
    T.alpha1 = -beta2;
    T.alpha2 = -beta1;
    T.deltaBeta = beta2-beta1;
    T.deltaBeta_deg = rad2deg(T.deltaBeta);
    
    T.beta1_deg = rad2deg(T.beta1);
    
    T.w1M = v_M;
    T.w1T = v_M * tan(T.beta1);
    T.w1 = sqrt(T.w1M^2 + T.w1T^2);
    
    T.w2M = v_M;
    T.w2T = v_M * tan(T.beta2);
    T.w2 = sqrt(T.w2M^2 + T.w2T^2);
    
    T.v1M = v_M;
    T.v1T = v_M * tan(T.alpha1);
    T.v1 = sqrt(T.v1M^2 + T.v1T^2);
    
    T.v2M = v_M;
    T.v2T = v_M * tan(T.alpha2);
    T.v2 = sqrt(T.v2M^2 + T.v2T^2);
    
end


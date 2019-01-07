function [ deltaBeta, beta1 ] = howell( howell_POLY, beta2 )
  
    deltaBeta = deg2rad(polyval(howell_POLY,abs(rad2deg(beta2))));
    
    beta1 = beta2 - deltaBeta;
end


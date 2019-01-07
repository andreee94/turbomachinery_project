function [ camber, camber_deg ] = camber_da_profilo(m, p)

% m = freccia massima percentuale di corda 
% p = posizione freccia massima

camber = atan(2*m./p)+atan(2*m./(1+p));

camber_deg = rad2deg(camber);

end


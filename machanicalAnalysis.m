function [ res ] = machanicalAnalysis( res )
    
    res.power = res.mdot * (res.points(1).hT);% - res.points(end).hT);
    res.omega = res.n * 2 * pi / 60;
    res.torque = res.power / res.omega;
    
    res.tau = res.torque * 16 / pi / res.Dm^3;
    
end


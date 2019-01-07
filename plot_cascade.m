function [ res ] = plot_cascade( res, num_profile )
    
    data = [];
    
    rotX = @(x,y,alpha) x*cos(alpha)-y*sin(alpha);
    rotY = @(x,y,alpha) x*sin(alpha)+y*cos(alpha);
    
    figure
    hold on
    axis equal
    grid on
    
    x_TR = 0;
    
    for ii = 1:length(res.stages)
        
        stage = res.stages(ii);
        
        alpha0 = inlineif(ii==1, 0, res.MID.velocity_triangle.alpha2);
        
        [P_stator, data] = generateBlade(data, alpha0, res.MID.velocity_triangle.alpha1, res.points(2 * ii - 1).c , stage.tmax_su_c, stage.statorSolidity);
        [P_rotor, data] = generateBlade(data, res.MID.velocity_triangle.beta1, res.MID.velocity_triangle.beta2, res.points(2 * ii).c , stage.tmax_su_c, stage.rotorSolidity);
        
        res.stages(ii).statorAlpha1Geom = P_stator.alpha2geom;
        res.stages(ii).rotorBeta2Geom = P_rotor.alpha2geom;
        
        res.stages(ii).statorAlpha1GeomDEG = rad2deg(P_stator.alpha2geom);
        res.stages(ii).rotorBeta2GeomDEG = rad2deg(P_rotor.alpha2geom);
        
        %[P_stator, data] = generateBlade(data, deg2rad(60), deg2rad(-30), res.points(2 * ii - 1).c , stage.tmax_su_c);
        %[P_rotor, data] = generateBlade(data, deg2rad(-30), deg2rad(60), res.points(2 * ii).c , stage.tmax_su_c);
        
        pitch_stator = res.points(2 * ii - 1).c / stage.statorSolidity * 100;%centimeter
        pitch_rotor = res.points(2 * ii).c / stage.statorSolidity * 100;%centimeter
        
        r_stator = -alpha0 - atan(P_stator.mC(1));
        r_rotor = -res.MID.velocity_triangle.beta1 - atan(P_rotor.mC(1));
        
        %r_stator =  deg2rad(60) - atan(P_stator.mC(1));
        %r_rotor =  deg2rad(-30) - atan(P_rotor.mC(1));
        
        y_TR_stator = 0;
        y_TR_rotor = 0;
        
        for kk = 1:num_profile
            plot(x_TR + rotX(P_stator.xUP, P_stator.yUP, r_stator), rotY(P_stator.xUP, P_stator.yUP, r_stator) + y_TR_stator, 'k');
            plot(x_TR + rotX(P_stator.x, P_stator.yC, r_stator), rotY(P_stator.x, P_stator.yC, r_stator) + y_TR_stator, 'k--');
            plot(x_TR + rotX(P_stator.xDOWN, P_stator.yDOWN, r_stator), rotY(P_stator.xDOWN, P_stator.yDOWN, r_stator) + y_TR_stator, 'k');
            y_TR_stator = y_TR_stator + pitch_stator;
        end
        
        x_TR = x_TR + res.points(2 * ii - 1).c * 100;%centimeter
        
        for kk = 1:num_profile
            plot(x_TR + rotX(P_rotor.xUP, P_rotor.yUP, r_rotor), rotY(P_rotor.xUP, P_rotor.yUP, r_rotor) + y_TR_rotor, 'r');
            plot(x_TR + rotX(P_rotor.x, P_rotor.yC, r_rotor), rotY(P_rotor.x, P_rotor.yC, r_rotor) + y_TR_rotor, 'r--');
            plot(x_TR + rotX(P_rotor.xDOWN, P_rotor.yDOWN, r_rotor), rotY(P_rotor.xDOWN, P_rotor.yDOWN, r_rotor) + y_TR_rotor, 'r');
            y_TR_rotor = y_TR_rotor + pitch_rotor;
        end
        
        x_TR = x_TR + res.points(2 * ii).c * 100;%centimeter
    end
    
    
end


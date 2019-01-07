function [ p_opt, fit_opt, res ] = MAIN_GA()
    
    
    % initial data
    gl.mdot = 200;
    gl.betaTT = 2.5;
    gl.TT0 = 700+273.15;
    gl.pT0 = 150;
    
    % perfect gas data
    gl.MM = 2 * 1 + 16;
    gl.R = 8314 / gl.MM;
    
    % mathematical data
    gl.toll = 1e-6;
    gl.maxiter = 50;
    gl.secant_entropy_delta0 = 0.01;% distance between x0 and x1 in secants algorithm
    gl.NGauss = 3;
    [gl.XGauss, gl.WGauss] = lgwt(gl.NGauss, 0, 1);
    
    % velocity triangle span evolution
    gl.velocity_evolution = @velocity_evolution_ConstantAngle;
    
    % design parameters
    % gl.chiMID = 0.4652;
    % gl.b_su_DmMIN = 0.05;
    % gl.n = 6000;
    % gl.usePhi = true;
    % gl.phi = 0.4465; % vA / u
    % gl.v1t_by_kp_su_u = 1; % kp = u / v1 = coef * u / v1t % This is an alternative to define phi
    % gl.autoNStages = false;
    % gl.nStages = 4;
    % gl.etaGuess = 0.93;
    % gl.epsilon = 1;
    % gl.chord_su_b = .5;
    % gl.solidity = 1.219;
    % gl.clearance = 0.02*0.025; % 0.02 * approx b
    % gl.seals_num = 2;
    % gl.tmax_su_c = 0.2;
    
    constraints = [0.2, 0.8;...
        0.03, 0.2;...
        3000, 18000;...
        0.3, 1.5;...
        1, 10;...
        0.5, 2;...
        0.5, 2;...
        0.0003, 0.009;
        1, 4;
        0.1,0.3];
    var_type = [myga_const.VAR_TYPE_DOUBLE;
        myga_const.VAR_TYPE_DOUBLE;
        myga_const.VAR_TYPE_INT;
        myga_const.VAR_TYPE_DOUBLE;
        myga_const.VAR_TYPE_INT;
        myga_const.VAR_TYPE_DOUBLE;
        myga_const.VAR_TYPE_DOUBLE;
        myga_const.VAR_TYPE_DOUBLE;
        myga_const.VAR_TYPE_INT
        myga_const.VAR_TYPE_DOUBLE];
    
    options.population_number = 50;
    options.selection_quantity = 0.350; % from 0 to 1
    options.mutation_min = 0.2;
    options.mutation_max = 0.8;
    options.mutation_threshold = 0.2;
    options.max_iterations = 10;
    options.crossover_strategy = myga_const.CROSSOVER_STATEGY_ALWAYS_BEST;
    options.crossover_first_threshold = 0.65;
    
    [ p_opt, fit_opt, res ] = myga(@func, constraints, var_type, options) ;
    
end

function [eta] = func(x)
    
    test.beta_from_n = false;
    test.b_from_Dm = false;
    test.optimal_velocity_triangle = false;
    test.hainley_mathieson_plot = false;
    test.print = true;
    gl.mdot = 200;
    gl.betaTT = 2.5;
    gl.TT0 = 700+273.15;
    gl.pT0 = 150;
    gl.toll = 1e-6;
    gl.maxiter = 50;
    gl.secant_entropy_delta0 = 0.01;% distance between x0 and x1 in secants algorithm
    gl.NGauss = 3;
    [gl.XGauss, gl.WGauss] = lgwt(gl.NGauss, 0, 1);
    
    % velocity triangle span evolution
    gl.velocity_evolution = @velocity_evolution_ConstantAngle;
 
    ff = 1;
    gl.chiMID = x(ff); ff = ff + 1;
    gl.b_su_DmMIN =  x(ff); ff = ff + 1;
    gl.n = x(ff); ff = ff + 1;
    gl.usePhi = true;
    gl.phi = x(ff); ff = ff + 1;
    gl.v1t_by_kp_su_u = 1; % kp = u / v1 = coef * u / v1t % This is an alternative to define phi
    gl.autoNStages = false;
    gl.nStages = x(ff); ff = ff + 1;
    gl.etaGuess = 0.9;
    gl.epsilon = 1;
    gl.chord_su_b = x(ff); ff = ff + 1;
    gl.solidity = x(ff); ff = ff + 1;
    gl.clearance = x(ff); ff = ff + 1;
    gl.seals_num = x(ff); ff = ff + 1;
    gl.tmax_su_c = x(ff); ff = ff + 1;
    res = solve_turbine(gl, test);
    eta = res.etaTT;
end


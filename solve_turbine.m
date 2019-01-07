function [ res ] = solve_turbine( gl, test )
    
    if gl.NGauss ~= length(gl.XGauss)
        [gl.XGauss, gl.WGauss] = lgwt(gl.NGauss, 0, 1);
    end
    
    res = [];
    
    func_b = @(mdot, rho, Dm, n, vA, epsilon) mdot / (rho * Dm * pi * vA * epsilon);
    
    iter.eta = 0;
    err.eta = inf;
    % first guess of some iterating variables
    temp.etaTT = gl.etaGuess;
    temp.rho0 = XSteam('rho_pt', gl.pT0, gl.TT0 - 273.15);
    
    % calculation of work per stage
    hT0 = XSteam('h_pt', gl.pT0, gl.TT0 - 273.15) * 1000;
    sT0 = XSteam('s_pt', gl.pT0, gl.TT0 - 273.15); % = s0
    hToutIS = XSteam('h_ps', gl.pT0 / gl.betaTT, sT0) * 1000;
    
    while abs(err.eta) > gl.toll && iter.eta <= gl.maxiter
        
        hTout = hT0 - temp.etaTT * (hT0 - hToutIS);
        if ~gl.autoNStages % TODO implement else case
            l = (hTout - hT0) / gl.nStages;
        end
        nStages = gl.nStages;
        
        % 1 - (v2t + v1t) / 2 / u = chi
        % l = u * (v1t - v2t)
        % mdot = rho * va * pi * b/dm * (60 / 2/pi/n)^2 * u^2
        % va / u = phi
        
        u = (gl.mdot / (gl.epsilon * temp.rho0 * gl.phi * pi * (gl.b_su_DmMIN) * (60 / 2 / pi / gl.n) ^ 2)) ^ (1/3);
        vA = gl.phi * u;
        deltaVT = -l / u; % v1t - v2t
        sumVT = 2 * u * (1 - gl.chiMID); % v1t + v2t
        lambda = -l / u^2;
        v1t = 0.5 * (deltaVT + sumVT);
        Dm = u / (2 * pi * gl.n / 60);
        
        %print(test, u, vA, deltaVT, sumVT, lambda, Dm)
        
        res.MID.velocity_triangle = velocity_triangle(u, vA, deltaVT, atan(v1t / vA));
        
        % v0 is vA since we have an axial inlet flow
        P0 = static_from_totalSTEAM(gl, gl.pT0, gl.TT0, vA, false);
        temp.rho0 = P0.rho;
        
        P0.b = func_b(gl.mdot, temp.rho0, Dm, gl.n, vA, gl.epsilon);
        P0.c = P0.b * gl.chord_su_b;
        
        %print(test, P0.b)
        
        res.points = [P0];
        res.stages = [];
        
        %         res.YTotStator = zeros(1, nStages);
        %         res.YTotRotor = zeros(1, nStages);
        %         chi_h = zeros(1, nStages);
        %         etas = zeros(1, nStages);
        %         betas = zeros(1, nStages);
        %         betasTT = zeros(1, nStages);
        
        iterating = true;
        
        for ii = 1 : nStages
            
            stage.chord_su_b = gl.chord_su_b;
            stage.clearance = gl.clearance;
            stage.seals_num = gl.seals_num;
            stage.tmax_su_c = gl.tmax_su_c;
            stage.epsilon = gl.epsilon;
            
            iter.stator = 0;
            err.b_stator = inf;
            temp.b = P0.b * 1.05; % first guess
            
            % stator iteration
            while err.b_stator > gl.toll && iter.stator <= gl.maxiter
                if ii == 1
                    is_first = true;
                else
                    P0 = P2;
                    is_first = false;
                end
                
                stage.statorNBlades = round(pi * Dm / (gl.chord_su_b * temp.b) * gl.solidity);
                stage.statorSolidity = stage.chord_su_b * temp.b / pi / Dm * stage.statorNBlades;
                
                is_mid = true;
                % compute losses in the mid span
                [stage.statorYTot, stage.statorYClearance, stage.statorYSec, stage.statorYProfile] = AinleyMathiesonLosses(test, iterating, 'stator', is_mid, is_first, res.MID.velocity_triangle, P0, stage.statorSolidity, stage.chord_su_b * temp.b, stage.clearance, stage.seals_num, stage.tmax_su_c, temp.b);
                % compute thermodynamical properties of point 1 at midspan
                P1 = solve_statorREAL(gl, iterating, P0, res.MID.velocity_triangle, stage.statorYTot);%, beta);
                
                is_mid = false;
                temp.r = gl.XGauss * temp.b + Dm / 2  - temp.b / 2; % min + (max - min) * percentage
                [temp.VA1, temp.VT1, temp.alpha1deg, temp.VA2, temp.VT2, temp.alpha2deg] = gl.velocity_evolution(res.MID.velocity_triangle, Dm / 2, 'stator', true, true, temp.r);% only evaluate point 1
                
                temp.u = res.MID.velocity_triangle.u / Dm * 2 * temp.r;
                temp.rho1 = zeros(gl.NGauss, 1);
                temp.s1 = zeros(gl.NGauss, 1);
                temp.solidity = stage.chord_su_b * temp.b / 2 / pi ./ temp.r * stage.statorNBlades;
                temp.Yprofile =zeros(gl.NGauss, 1);
                for kk = 1:gl.NGauss
                    % create triangle from span evolution
                    temp.triangle = temp_triangle_stator_func(temp.VA1(kk), temp.VT1(kk), temp.alpha1deg(kk), temp.VA2(kk), temp.VT2(kk), temp.alpha2deg(kk));
                    % compute the new profile losses
                    temp.Yprofile(kk) = AinleyMathiesonLosses(test, iterating, 'stator', is_mid, is_first, temp.triangle, P0, temp.solidity(kk), stage.chord_su_b * temp.b, stage.clearance, stage.seals_num, stage.tmax_su_c, temp.b);
                    % compute the new density along the span
                    [temp.rho1(kk), temp.s1(kk)] = solve_stator_streamline(gl, P0, temp.triangle.v1, temp.Yprofile(kk) + stage.statorYClearance + stage.statorYSec);
                end
                
                % compute the blade height with a weighted average along the span
                % sum[mdot / (pi * dm * rho(r_i) * va(r_i)) * w_i]
                b = (gl.mdot / pi / Dm ./ temp.rho1 ./ temp.VA1 / stage.epsilon)' * gl.WGauss;
                
                % evaluate error with prevoius iteration and update temp.b
                err.b_stator = abs(b - temp.b);
                temp.b = b;
                iter.stator = iter.stator + 1;
            end % while err.b_stator > gl.toll
            
            P1.b = temp.b;
            P1.c = P1.b * stage.chord_su_b;
            P1.span.r = temp.r;
            P1.span.u = temp.u;
            P1.span.vA1 = temp.VA1;
            P1.span.vT1 = temp.VT1;
            P1.span.alpha1 = temp.alpha1deg;
            
            P1.span.vA2 = temp.VA2;
            P1.span.vT2 = temp.VT2;
            P1.span.alpha2 = temp.alpha2deg;
            
            P1.span.rho = temp.rho1;
            P1.span.s = temp.s1;
            P1.span.s_mean = P1.span.s' * gl.WGauss;
            P1.span.solidity = temp.solidity;
            stage.span.statorYprofile = temp.Yprofile;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            iter.rotor = 0;
            err.b_rotor = inf;
            temp.b = P1.b * 1.05; % first guess
            P2 = P1; %first guess
            
            % rotor iteration
            while err.b_rotor > gl.toll && iter.rotor <= gl.maxiter
                is_first = false;
                
                stage.rotorNBlades = round(pi * Dm / (stage.chord_su_b * temp.b) * gl.solidity);
                stage.rotorSolidity = stage.chord_su_b * temp.b / pi / Dm * stage.rotorNBlades;
                
                is_mid = true;
                % compute losses in the mid span
                [stage.rotorYTot, stage.rotorYClearance, stage.rotorYSec, stage.rotorYProfile] = AinleyMathiesonLosses(test, iterating, 'rotor', is_mid, is_first, res.MID.velocity_triangle, P1, stage.rotorSolidity, stage.chord_su_b * temp.b, stage.clearance, stage.seals_num, stage.tmax_su_c, temp.b);
                % compute thermodynamical properties of point 1 at midspan
                P2 = solve_rotorREAL(gl, iterating, P1, res.MID.velocity_triangle, stage.rotorYTot);%, beta);
                
                is_mid = false;
                temp.r = gl.XGauss * temp.b + Dm / 2  - temp.b / 2; % min + (max - min) * percentage
                [temp.WA1, temp.WT1, temp.beta1deg, temp.WA2, temp.WT2, temp.beta2deg] = gl.velocity_evolution(res.MID.velocity_triangle, Dm / 2, 'rotor', true, true, temp.r);% only evaluate point 1
                
                temp.u = res.MID.velocity_triangle.u / Dm * 2 * temp.r;
                temp.rho1 = zeros(gl.NGauss, 1);
                temp.s1 = zeros(gl.NGauss, 1);
                temp.solidity = stage.chord_su_b * temp.b / 2 / pi ./ temp.r * stage.rotorNBlades;
                temp.Yprofile =zeros(gl.NGauss, 1);
                for kk = 1:gl.NGauss
                    % create triangle from span evolution
                    temp.triangle = temp_triangle_rotor_func(temp.WA1(kk), temp.WT1(kk), temp.beta1deg(kk), temp.WA2(kk), temp.WT2(kk), temp.beta2deg(kk));
                    % compute the new profile losses
                    temp.Yprofile(kk) = AinleyMathiesonLosses(test, iterating, 'rotor', is_mid, is_first, temp.triangle, P1, temp.solidity(kk), stage.chord_su_b * temp.b, stage.clearance, stage.seals_num, stage.tmax_su_c, temp.b);
                    % compute the new density along the span
                    [temp.rho1(kk), temp.s1(kk)] = solve_rotor_streamline(gl, P1, temp.triangle.w1, temp.Yprofile(kk) + stage.rotorYClearance + stage.rotorYSec);
                end
                
                % compute the blade height with a weighted average along the span
                % sum[mdot / (pi * dm * rho(r_i) * va(r_i)) * w_i]
                b = (gl.mdot / pi / Dm ./ temp.rho1 ./ temp.WA1 / stage.epsilon)' * gl.WGauss;
                
                % evaluate error with prevoius iteration and update temp.b
                err.b_rotor = abs(b - temp.b);
                temp.b = b;
                iter.rotor = iter.rotor + 1;
            end % while err.b_rotor > gl.toll
            
            P2.b = temp.b;
            P2.c = P2.b * stage.chord_su_b;
            P2.span.r = temp.r;
            P2.span.u = temp.u;
            P2.span.wA1 = temp.WA1;
            P2.span.wT1 = temp.WT1;
            P2.span.beta1 = temp.beta1deg;
            
            P2.span.wA2 = temp.WA2;
            P2.span.wT2 = temp.WT2;
            P2.span.beta2 = temp.beta2deg;
            
            P2.span.rho = temp.rho1;
            P2.span.s = temp.s1;
            P2.span.s_mean = P2.span.s' * gl.WGauss;
            P2.span.solidity = temp.solidity;
            stage.span.rotorYprofile = temp.Yprofile;
            
            % stage performances
            stage.h2IS = XSteam('h_ps', P2.p, P0.s) * 1000;
            
            stage.chi_h = (P1.h - P2.h) / (P0.h - P2.h);
            stage.eta = (P0.h - P2.h) / (P0.h - stage.h2IS);
            stage.beta = P0.p / P2.p;
            stage.betaTT = P0.pT / P2.pT;
            stage.betaStator = P0.p / P1.p;
            stage.betaTTStator = P0.pT / P1.pT;
            stage.betaRotor = P1.p / P2.p;
            stage.betaTTRotor = P1.pT / P2.pT;
            stage.ws = 2 * pi * gl.n / 60 * gl.mdot / P0.rho / (P0.hT - P2.hT)^0.75;
            stage.ds = Dm / gl.mdot * P0.rho * (P0.hT - P2.hT)^0.25;
            stage.workspan = (P2.span.u .* (P2.span.wT2 - P2.span.wT1) .* temp.rho1 .* temp.WA1 * pi * Dm * temp.b)' * gl.WGauss / gl.mdot;
            stage.l = l;
            
            % bad way to concat
            res.stages = [res.stages, stage];
            
            res.points = [res.points, [P1, P2]];
        end % for ii = 1 : nStages
        
        hTIS = XSteam('h_ps', res.points(end).pT, res.points(1).s) * 1000;
        res.etaTT =  (res.points(1).hT - res.points(end).hT) / (res.points(1).hT - hTIS);
        
        res.beta = res.points(1).p / res.points(end).p;
        res.betaTT = res.points(1).pT / res.points(end).pT;
        
        err.eta = abs(res.etaTT - temp.etaTT);
        temp.etaTT = res.etaTT;
        
        iter.eta = iter.eta + 1;
        
        %etaTT = temp.etaTT;
        %error = err.eta;
        %print(test, etaTT, error);
    end
    
    res.sumb = sum([res.points.b]);
    res.workerrors = ([res.stages.workspan] - l) / l;
    res.workerror = mean(res.workerrors);
    res.workerrorperc = abs(res.workerror * 100);
    res.hTend_real = hT0 - abs(sum([res.stages.workspan]));
    res.beta_real = gl.pT0 / XSteam('p_hs', res.hTend_real / 1000, res.points(end).s);
    
    res.u = u;
    res.n = gl.n;
    res.vA = vA;
    res.deltaVT = deltaVT;
    res.sumVT = deltaVT;
    res.lambda = lambda;
    res.Dm = Dm;
    res.mdot = gl.mdot;
    res.stagesPower = res.mdot * (res.points(1).hT - res.points(end).hT);
    
    if iter.eta > gl.maxiter && err.eta > 2.5 * gl.toll
        res.etaTT = 0;
        myprint(test, 'Main eta cycle does not converge')
    end
    
end


function res = temp_triangle_stator_func(v1a, v1t, alpha1, v2a, v2t, alpha2)
    res.v1A = v1a;
    res.v1T = v1t;
    res.v1 = norm([v1a, v1t]);
    res.alpha1_deg = alpha1;
    
    res.v2A = v2a;
    res.v2T = v2t;
    res.v2 = norm([v2a, v2t]);
    res.alpha2_deg = alpha2;
end

function res = temp_triangle_rotor_func(v1a, v1t, alpha1, v2a, v2t, alpha2)
    res.w1A = v1a;
    res.w1T = v1t;
    res.w1 = norm([v1a, v1t]);
    res.beta1_deg = alpha1;
    
    res.w2A = v2a;
    res.w2T = v2t;
    res.w2 = norm([v2a, v2t]);
    res.beta2_deg = alpha2;
end


function [ YTot, YClearance, YSec, YpReval, YpRefval, Yp1val, Yp2val ] = AinleyMathiesonLosses(test, iterating, is_stator, is_mid, is_first, triangle, P1, solidity, chord, clearance, seals_num, tmax_su_c, blade_height)
    
    
    %AinleyMathiesonLosses(test, iterating, is_mid, alpha0, alpha1, solidity, tmax_su_c, rho, viscosity, chord, v1, clearance, blade_height)
    
    %AinleyMathiesonLosses(test, iterating, is_mid, is_first, triangle, P1, solidity, chord, clearance, tmax_su_c, blade_height)
    %test, iterating, is_mid, alpha0, alpha1, solidity, tmax_su_c, rho, viscosity, chord, v1, clearance, blade_height)
    if test.hainley_mathieson_plot == false
        if strcmp('stator', is_stator)
            alpha0 = inlineif(is_first, 0, triangle.alpha2_deg);
            alpha1 = triangle.alpha1_deg;
            viscosity = P1.mu; %mean([P1.mu, P2.mu]);
            rho = P1.rho;
            v1 = mean([inlineif(is_first, triangle.v1A, triangle.v2), triangle.v1]);
            clearance = 0;
        else
            alpha0 = inlineif(is_first, 0, triangle.beta1_deg);
            alpha1 = triangle.beta2_deg;
            viscosity = P1.mu; %mean([P1.mu, P2.mu]);
            rho = P1.rho;
            v1 = mean([triangle.w2, triangle.w1]);
        end
        
        
        %     if nargin == 1
        %         iterating = tes.iterating;
        %         is_mid = test.is_mid;
        %         alpha0 = test.alpha0;
        %         alpha1 = test.alpha1;
        %         solidity = test.solidity;
        %         tmax_su_c = test.tmax_su_c;
        %         rho = test.rho;
        %         viscosity = test.viscosity;
        %         chord = test.chord;
        %         v1 = test.v1;
        %         clearance = test.clearance;
        %         blade_height = test.blade_height;
        %         test = test.test;
        %     end
        
        % interpolations for alpha0 = 0
        %     A1 = @(a1p) inlineif(a1p < 27, 0.025 + (27 - a1p) / 530, 0.025 + (27 - a1p) / 3085);
        %
        %     B1 = @(a1p) inlineif(a1p <= 30, 0.1583 - a1p / 1640, 0.1583 - a1p / 1640);
        %
        %     C1 = @(a1p) 0.08 * ( (a1p/30)^2 -1);
        %
        %     n1 = @(a1p) 1 + a1p / 30;
        %
        %     min_inv_solidity1 = @(a1p) inlineif(a1p <= 30, 0.46 + a1p / 77, 0.614 + a1p / 130);
        %
        %     X1 = @(a1p, solidity) 1 ./ solidity - min_inv_solidity1(a1p);
        %
        %     Yp1 = @(a1p, solidity) inlineif(a1p <= 30, ...
        %         A1(a1p) + B1(a1p) .* X1(a1p, solidity).^2 + C1(a1p) .* X1(a1p, solidity).^3, ...
        %         A1(a1p) + B1(a1p) .* abs(X1(a1p, solidity).^n1(a1p)));
        %
        %
        %     % interpolations for alpha0 = - alpha1
        %
        %     A2 = @(a1p) 0.242 - a1p / 151 + (a1p / 127)^2;
        %
        %     B2 = @(a1p) inlineif(a1p <= 30, 0.3 + (30 - a1p) / 50, 0.3 + (30 - a1p) / 275);
        %
        %     C2 = @(a1p) 0.88 - a1p / 42.4 + (a1p / 72.8)^2;
        %
        %     min_inv_solidity2 = @(a1p) 0.224 + 1.575 * (a1p / 90) - (a1p / 90)^2;
        %
        %     X2 = @(a1p, solidity) 1 ./ solidity - min_inv_solidity2(a1p);
        %
        %     Yp2 = @(a1p, solidity) A2(a1p) + B2(a1p) .* X2(a1p, solidity).^2 - C2(a1p) .* X2(a1p, solidity).^3;
        
        % compute the losses for a specific value of alpha0 and alpha1 and
        % solidity
        
        Yp1val = Yp1(90 - abs(alpha1), solidity);
        Yp2val = Yp2(90 - abs(alpha1), solidity);
        
        ma = -(alpha0 / alpha1);
        if ma >= 1.01
            print(test, 'ma > 1')
            ma = 1 / ma;
        end
        
        % reference value
        YpRefval = (Yp1val + ma^2 * (Yp2val - Yp1val)) * (tmax_su_c / 0.2)^ma;
        
        %reynolds correction
        Re = rho * v1 * chord / viscosity;
        
        YpReval = YpRefval * (2e5 ./ Re).^0.2;
        
        if YpReval > 1e4
            print(test, 'Y > 1e4')
        end
        
        if ~is_mid
            YTot = YpReval;
        end
        
        % in case of tip or hub it is not required
        if is_mid
            
            %secondary losses and clearance
            alpham = atand( 0.5 * (tand(alpha0) + tand(alpha1)));
            
            %cL = 2 ./ solidity * abs(tand(alpha0) - tand(alpha1)) * cosd(alpham);
            
            %coeff =  chord ./ blade_height .* (cL .* solidity).^2 * cosd(alpha1)^2 / cosd(alpham)^3;
            
            coeff = 4 * chord ./ blade_height .* (tand(alpha0) - tand(alpha1)).^2 .* cosd(alpha1).^2 ./ cosd(alpham);
            
            YSec = 0.0334 * cosd(alpha1) / cosd(alpha0) * coeff;
            
            shrouded = true;
            
            B = inlineif(shrouded, 0.37, 0.47);
            
            clearance = clearance / seals_num ^ 0.42;
            
            YClearance = B * (clearance ./ chord).^0.78 .* coeff;
            % total losses
            
            YTot = YpReval + YClearance + YSec;
            
            
            
        else
            % YTot = [];
            YClearance = [];
            YSec = [];
            YpReval = [];
            YpRefval = [];
            Yp1val = [];
            Yp2val = [];
        end
    end
    
    %% ---------------------------------------
    % plot to compare with real graphs
    
    if test.hainley_mathieson_plot && ~iterating
        
        %         myprint(Yp1val)
        %         myprint(Yp2val)
        %         myprint(YpRefval)
        %         myprint(YSec)
        %         myprint(YClearance)
        
        f1 = figure;
        hold on
        
        f2 = figure;
        hold on
        
        n1 = get(f1,'Number');
        n2 = get(f2,'Number');
        
        inv_solidity = linspace(0.2, 1.2, 100);
        Yp1_res = zeros(size(inv_solidity));
        Yp2_res = zeros(size(inv_solidity));
        
        for a1p = [10 15 20 25 30 40 50]
            kk = 1;
            for ii = inv_solidity
                Yp1_res(kk) = Yp1(a1p, ii^-1);
                Yp2_res(kk) = Yp2(a1p, ii^-1);
                kk = kk + 1;
            end
            
            figure(n1)
            plot(inv_solidity, Yp1_res, 'LineWidth', 1.5)
            
            figure(n2)
            plot(inv_solidity, Yp2_res, 'LineWidth', 1.5)
            
        end
        figure(n1)
        title('Profile losses (axial inlet)')
        xlabel('s / c = 1 / solidity')
        ylabel('Profile losses Y_p')
        legend('10', '15', '20', '25', '30', '40', '50')
        
        figure(n2)
        title('Profile losses (symmetric blades)')
        xlabel('s / c = 1 / solidity')
        ylabel('Profile losses Y_p')
        legend('10', '15', '20', '25', '30', '40', '50')
    end
end

function res = A1(a1p)
    res = inlineif(a1p < 27, 0.025 + (27 - a1p) / 530, 0.025 + (27 - a1p) / 3085);
end

function res = B1(a1p)
    res = inlineif(a1p <= 30, 0.1583 - a1p / 1640, 0.1583 - a1p / 1640);
end

function res = C1(a1p)
    res = 0.08 * ( (a1p/30)^2 -1);
end

function res = n1(a1p)
    res =  1 + a1p / 30;
end

function res = min_inv_solidity1(a1p)
    res = inlineif(a1p <= 30, 0.46 + a1p / 77, 0.614 + a1p / 130);
end

function res = X1(a1p, solidity)
    res = 1 ./ solidity - min_inv_solidity1(a1p);
end

function res = Yp1(a1p, solidity)
    res = inlineif(a1p <= 30, ...
        A1(a1p) + B1(a1p) .* X1(a1p, solidity).^2 + C1(a1p) .* X1(a1p, solidity).^3, ...
        A1(a1p) + B1(a1p) .* abs(X1(a1p, solidity).^n1(a1p)));
end

% interpolations for alpha0 = - alpha1

function res = A2(a1p)
    res = 0.242 - a1p / 151 + (a1p / 127)^2;
end

function res = B2(a1p)
    res = inlineif(a1p <= 30, 0.3 + (30 - a1p) / 50, 0.3 + (30 - a1p) / 275);
end

function res = C2(a1p)
    res = 0.88 - a1p / 42.4 + (a1p / 72.8)^2;
end

function res = min_inv_solidity2(a1p)
    res = 0.224 + 1.575 * (a1p / 90) - (a1p / 90)^2;
end

function res = X2(a1p, solidity)
    res = 1 ./ solidity - min_inv_solidity2(a1p);
end

function res = Yp2(a1p, solidity)
    res = A2(a1p) + B2(a1p) .* X2(a1p, solidity).^2 - C2(a1p) .* X2(a1p, solidity).^3;
end

function [ P, data ] = generateBlade( data, alpha1, alpha2, c, tmax_su_c, solidity )
    
    if isempty(data)
        data = importdata('graph/A3K7.csv', ';', 1);
    end
    
    
     f = @(x) blade_func(x, alpha2, alpha1, tmax_su_c, c, solidity);
     
%     xxx = deg2rad(linspace(10,100,1000));
%     fxxx = zeros(size(xxx));
%     
%     for ii = 1:length(xxx)
%         fxxx(ii) = f(xxx(ii));
%     end
%     
%     figure(2)
%     grid on
%     plot(xxx, fxxx);
    
    alpha2_new = secants(f, sign(alpha2)*(abs(alpha2)), sign(alpha2)*(abs(alpha2)+0.1));
    
    alpha2_new_deg = rad2deg(alpha2_new);
    
    %     temp.err = inf;
    %     temp.old = inf;
    %     temp.alpha2 = alpha2;
    
    
    profile = create_profile(alpha2_new - alpha1, tmax_su_c * c);
    profile.alpha2flow = alpha2;
    profile.alpha2geom = alpha2_new;
    
    alpha2 = alpha2_new;
    
    c = c * 100;% centimeter
    
    profile.x =  profile.x * c;
    profile.yC =  profile.yC * c;
    profile.rC =  profile.rC * c;
    
    x = data.data(:, 1) / 100 * c; % normalized to 1
    yup = data.data(:, 3); % normalized to 1
    
    %x = x; % centimeter
    yup = yup / max(yup) * tmax_su_c / 2  * c;
    %ydown = - data.data(:, 3) / 100;
    
    %     P = polyfit(x, yup, 5);
    %
    %     test.print = true;
    %
    %     for ii = 1:length(x)
    %         myprint(test, [num2str(x(ii)), ' ', num2str(yup(ii))]);
    %     end
    %
    %     t = polyval(P, profile.x );
    
    % perfetto
    t = spline(x, yup, profile.x) ;
    profile.xUP = profile.x - t .* sin(atan(profile.mC));
    profile.yUP = profile.yC + t .* cos(atan(profile.mC));
    profile.xDOWN = profile.x + t .* sin(atan(profile.mC));
    profile.yDOWN = profile.yC - t .* cos(atan(profile.mC));
    % end perfetto
    
%     abc = [profile.xUP(95) profile.yUP(95);
%         profile.xUP(97) profile.yUP(97);
%         profile.xUP(99) profile.yUP(99)];
%     
    %profile.e = fit_circle_through_3_points(abc);
%     
    %profile.e = profile.rC(100);
%     
    %alpha2 = sign(alpha2) *deg2rad(interp1([35 70], [30 70], abs(rad2deg(alpha2)), 'linear', 'extrap') - 4 * c / solidity / profile.e);
    
    profile.alpha1 = alpha1;
    profile.alpha2 = alpha2;
    profile.T = t;
    
    P = profile;
    
    
    %     figure
    %     grid on
    %     hold on
    % %     plot(profile.x, profile.yC + t, 'k')
    % %     plot(profile.x, profile.yC -t, 'k')
    %
    %     % perfetto
    %     plot(profile.xUP, profile.yUP, 'r')
    %     plot(profile.xDOWN, profile.yDOWN, 'r')
    %     % end perfetto
    %
    %      plot(profile.x, profile.yC , '.k')
    % %
    % %     plot(x, yup, '--k')
    % %     plot(x, ydown, '--k')
    %     axis equal
    
    
end

function res = blade_func(alpha2,alpha2desired, alpha1, tmax_su_c, c, solidity)
    
    profile = create_profile(alpha2 - alpha1, tmax_su_c * c);
    %profile.rC =  profile.rC * c;
    %e = profile.rC(end);
    e = profile.e;
    %alpha2deg = rad2deg(alpha2);
    %alpha2degdesired = rad2deg(alpha2desired);
    %m = profile.m;
    
    res = alpha2desired - sign(alpha2) * deg2rad(interp1([35 70], [30 70], abs(rad2deg(alpha2)), 'linear', 'extrap') - 4 * c / solidity / e);
    
end






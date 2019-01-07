function [X,Y,Z ] = plot_profile_single( P,beta1,z, doNotPlot)
    
    y_TR = 0;
    r = 0;
    Z = [];
    plot2D = false;
    
    
    if nargin >1
        r = -beta1 - atan(2*P.m/P.p);
        %y_TR = P.x*tan(-beta1);
    end
    if nargin < 3
        plot2D = true;
    end
    
    if nargin < 4
        doNotPlot = false;
    end
    
    rotX = @(x,y,alpha) x*cos(alpha)-y*sin(alpha);
    rotY = @(x,y,alpha) x*sin(alpha)+y*cos(alpha);
    
    if ~doNotPlot
        if plot2D
            plot(rotX(P.xU, P.yU, r), rotY(P.xU, P.yU, r) + y_TR, 'k');
            hold on
            axis equal
            %grid on
            %mid line
            %plot(rotX(P.x, P.yC, r), rotY(P.x, P.yC, r) + y_TR, 'k--');
            plot(rotX(P.xL, P.yL, r), rotY(P.xL, P.yL, r) + y_TR, 'k');
        else
            plot3(rotX(P.xU, P.yU, r), rotY(P.xU, P.yU, r) + y_TR, ones(size(P.xU))*z);%, 'k');
            hold on
            axis equal
            %grid on
            %mid line
            %plot(rotX(P.x, P.yC, r), rotY(P.x, P.yC, r) + y_TR, 'k--');
            plot3(rotX(P.xL, P.yL, r), rotY(P.xL, P.yL, r) + y_TR, ones(size(P.xU)));%*z, 'k');
        end
    end
    X = [rotX(P.xU, P.yU, r), flip(rotX(P.xL, P.yL, r))];
    Y = [rotY(P.xU, P.yU, r), flip(rotY(P.xL, P.yL, r))];
    if ~plot2D
        Z = ones(size(X))*z;
    end
end
function [ ] = plot_profile( P,beta1,num_profile, statore)
    
    y_TR = 0;
    r = 0;
    
    if nargin >1
        r = -beta1 - atan(2*P.m/P.p);
        %y_TR = P.x*tan(-beta1);
    end
    if nargin <4
        statore = true;
    end
    if nargin <3
        num_profile = 1;
    end
    
    rotX = @(x,y,alpha) x*cos(alpha)-y*sin(alpha);
    rotY = @(x,y,alpha) x*sin(alpha)+y*cos(alpha);
    
    for ii = 1:num_profile
        %rotore
        plot(rotX(P.xU, P.yU, r), rotY(P.xU, P.yU, r) + y_TR, 'k');
        hold on
        axis equal
        grid on
        plot(rotX(P.x, P.yC, r), rotY(P.x, P.yC, r) + y_TR, 'k--');
        plot(rotX(P.xL, P.yL, r), rotY(P.xL, P.yL, r) + y_TR, 'k');
        
        %dx = rotX(P.x, P.yC, r);
        %dy = rotY(P.x, P.yC, r)+y_TR;
        
        %plot([-0.3, dx(1)], [dy(1)-0.3*tan(-beta1), dy(1)], 'b')
        %non so come disegnarlo
        %plot(dx(end) + P.x, dy(end) + tan(-beta1+atan(2*P.m/(1+P.p)))*P.x)
        
        %statore
        if statore
            plot(1 + rotX(P.xU, P.yU, r),  y_TR - rotY(P.xU, P.yU, r), 'b');
            plot(1 + rotX(P.x, P.yC, r),  y_TR - rotY(P.x, P.yC, r), 'b--');
            plot(1 + rotX(P.xL, P.yL, r), y_TR - rotY(P.xL, P.yL, r), 'b');
        end
        
        y_TR = y_TR + tan(-beta1);
    end
end


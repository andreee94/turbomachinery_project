function [  ] = plot_machine( res, points )
    
    space = 1e-3; % 1mm
    
    n = length(points);
    
    figure
    hold on
    daspect([1 1 1])
    colors = get(gca,'colororder');
    
    %% plotting stators and rotors
    
    clearence = 0.97;
    translation = 0;
    
    for ii = 1:n-1
        
        p = points(ii);
        pnext = points(ii + 1);
        
        if mod(ii, 2) == 1 % stator
            tempx = [0, 0, p.c * cos(res.triangle.alpha1), p.c * cos(res.triangle.alpha1)];
            tempy = [p.b/2, -p.b/2 * clearence, -pnext.b/2 * clearence, pnext.b/2];
            color = colors(1,:);
        end
        
        if mod(ii, 2) == 0 % stator
            tempx = [0, 0, p.c * cos(res.triangle.alpha1), p.c * cos(res.triangle.alpha1)];
            tempy = [-p.b/2, p.b/2 * clearence,  pnext.b/2 * clearence, -pnext.b/2];
            color = colors(2,:);
        end
        
        fill(tempx + translation, res.Dm + tempy, color+0.1, 'LineWidth', 1);
        
        translation = translation + space + tempx(end);
        
    end
    
    
    %% plotting of the casing
    
    xx = zeros(n + 2, 1);
    yy = zeros(n + 2, 1);
    
    kk = 1;
    xx(kk) = 0;
    yy(kk) = 0;
    kk = kk + 1;
    
    
    translation = 0;
    for p = points(1:end)
        
        if kk > 2
            xx(kk) = translation;
            yy(kk) = p.b / 2;
            kk = kk + 1;
            translation = translation + space;
        end
        
        xx(kk) = translation;
        yy(kk) = p.b / 2;
        kk = kk + 1;
        translation = translation + p.c * cos(res.triangle.alpha1);
        
    end
    
    xx(end) = xx(end - 1);
    yy(end) = 0;
    
    plot(xx, res.Dm + yy, 'k', 'LineWidth', 3);
    plot(xx, res.Dm - yy, 'k', 'LineWidth', 3);
    
    plot(xx(2:end-1), res.Dm + yy(2:end-1), '.k', 'MarkerSize', 30);
    plot(xx(2:end-1), res.Dm - yy(2:end-1), '.k', 'MarkerSize', 30);
    
    %plot(xx, -res.Dm + yy);
    %plot(xx, -res.Dm - yy);
    
    
end




























function [ ] = plot_pTrho_machine( res )
    
    space = 1e-3; % 1mm
    
    n = length(res.points);
    
    ii_array = 1:3;
    
    for ii = ii_array
        
        figure
        hold on
        daspect([1 1 1])
        %colors = get(gca,'colororder');
        
        %% plotting of the casing
        xx = zeros(n + 2, 1);
        yy = zeros(n + 2, 1);
        
        kk = 1;
        xx(kk) = 0;
        yy(kk) = 0;
        kk = kk + 1;
        
        
        translation = 0;
        for p = res.points
            
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% plotting of the quantity
        
        m = 10; % points per section
        
        xx = zeros(m, 2 * n - 2);
        yy = zeros(size(xx));
        value = zeros(size(xx));
        
        kk = 1;
        
        translation = 0;
        for p = res.points
            
            if kk > 1
                xx(:, kk) = translation;
                yy(:, kk) = res.Dm + linspace(-p.b / 2, p.b / 2, m);
                value(:, kk) = value(:, kk-1);
                kk = kk + 1;
                translation = translation + space;
            end
            
            if kk < 2 * n  - 2
                xx(:, kk) = translation;
                yy(:, kk) = res.Dm + linspace(-p.b / 2, p.b / 2, m);
                value(:, kk) = inlineswitch(ii, ii_array, [p.p, p.T, p.rho]);
                kk = kk + 1;
                translation = translation + p.c * cos(res.triangle.alpha1);
            end
            
        end
        
        surf(xx, yy, value);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
end


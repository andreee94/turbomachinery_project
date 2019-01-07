function [xmin, xmax] = plotMinMax(x, y, plot_min, plot_max)
    
    plot(x, y)
    hold on
    
    if plot_min
        min_mask = y==min(y);
        xmin = x(min_mask);
        plot(xmin, y(min_mask), 'ro', 'LineWidth',2);
    end
    
    if plot_max
        max_mask = y==max(y);
        xmax = x(max_mask);
        plot(xmax, y(max_mask), 'ro', 'LineWidth',2);
    end
end


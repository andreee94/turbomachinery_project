function [xmin, xmax] = plotMinMax(x, y, plot_min, plot_max)
    
    h = plot(x, y, 'LineWidth', 1.5);
    color = get(h, 'Color');
    hold on
    grid on
    
    if plot_min
        min_mask = y==min(y);
        xmin = x(min_mask);
        plot(xmin, y(min_mask), 'o', 'LineWidth', 2, 'Color', color);
    end
    
    if plot_max
        max_mask = y==max(y);
        xmax = x(max_mask);
        plot(xmax, y(max_mask), 'o', 'LineWidth', 2, 'Color', color);
    end
end


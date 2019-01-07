function [ ] = plot_span_velocity( span )
    figure
    grid on
    hold on
    yyaxis left
    ylabel('Velocity (m/s)')
    plot(span.r, span.vA1, '-', 'LineWidth', 2)
    plot(span.r, span.vT1,':', 'LineWidth', 2)
    plot(span.r, span.u, '--', 'LineWidth', 2)
    
    yyaxis right
    ylabel('Angle (degrees)')
    plot(span.r, (span.alpha1),'-', 'LineWidth', 2)
    plot(span.r, (span.alpha2),'--', 'LineWidth', 2)
    
    legend('v_{A1}', 'v_{T1}', 'U', 'alpha0', 'alpha_1')
    xlabel('Radius (m)')
   
%     [ax, h1, h2 ] = plotyy(span.r, [span.vA1, span.vT1, span.u] , span.r, [span.alpha1, span.alpha2]);
%     [ax, h1, h2 ] = plotyy(span.r, [span.vA1] , span.r, [span.alpha1]);
%     [ax, h1, h2 ] = plotyy(span.r, [span.vA2] , span.r, [span.alpha2]);
%     
%     h1(1).LineWidth = 2;
%     h1(2).LineWidth = 2;
%     h1(3).LineWidth = 2;
%     h1(1).LineStyle = '-';
%     h1(2).LineStyle = ':';
%     h1(3).LineStyle = '--';
%     
%     ax(1)
%     
%     xlabel(ax(1), 'Radius (m)')
%     xlabel(ax(2), 'Radius (m)')
%     
%     ylabel(ax(1), 'Velocity (m/s)');
%     ylabel(ax(2), 'Angle (degrees)')
%     
%     legend('v_{A1}', 'v_{T1}', 'U', 'alpha_1', 'alpha_2')
end


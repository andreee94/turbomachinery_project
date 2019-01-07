function [ T0, h0 ] = static_temperature( TT0, v0, p0, cp0 )
    
    %h = @(T) XSteam('h_pT', p0, T - 273.15) * 1000;
    
    
    cp = @(T) XSteam('Cp_pT', p0, T - 273.15) * 1000;
    
    if nargin < 4
    %    cp = @(T) XSteam('Cp_pT', p0, T - 273.15) * 1000;
        cp0 = cp(TT0);
    end
    
    %f = @(T) cp(T) * T + 0.5 * v0^2 - cp(T) * TT0;
    f = @(T) T + 0.5 * v0^2 / cp(T) - TT0;
    
    %     tt = [900:1:1000];
    %     cps = zeros(size(tt));
    %     j = 1;
    %     for ttt = tt
    %         cps(j) = f(ttt);
    %         j = j + 1;
    %     end
    %     plot(tt, cps)
    
    
    % as first guesses we take
    % TT0
    % TT0 - v0^2 / 2 / cp0;
    % usually we are very close to the guess
    x0 = TT0 - 0.5 * v0^2 / cp0;
    T0 = secants(f, x0, x0 + 0.1 );
    h0 = XSteam('h_pT', p0, T0 - 273.15) * 1000;%h(T0);
    %cp0 = cp(T0);
    
end


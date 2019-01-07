
function [] = plotTS(res)
    
    grid on;
    ss = [res.points.s];
    TT = [res.points.T];
    plot(ss, TT, 'LineWidth', 1.5);
    hold on
    grid on
    curveS_liq = zeros(1,100);
    curveS_vap = zeros(1,100);
    curveT_liq = linspace(0,373.5,100);
    curveT_vap = linspace(373.5,0,100);
    kk = 1;
    for T = curveT_liq
        curveS_liq(kk) = XSteam('sL_T',T);
        kk = kk+1;
    end
    kk = 1;
    for T = curveT_vap
        curveS_vap(kk) = XSteam('sV_T',T);
        kk = kk+1;
    end
    plot([curveS_liq,curveS_vap(1),curveS_vap],[curveT_liq,curveT_vap(1),curveT_vap],  'LineWidth', 1.5);
    xlabel('Entropy $s$')
    ylabel('Temperature $T(K)$')
end
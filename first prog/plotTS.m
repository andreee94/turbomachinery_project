ss = [result.points.s];
TT = [result.points.T];
plot(flip(ss), TT)
hold on
curveS_liq = zeros(1,500);
curveS_vap = zeros(1,500);
curveT_liq = linspace(0,374,500);
curveT_vap = linspace(374,0,500);
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
plot([curveS_liq,curveS_vap],[curveT_liq,curveT_vap]);

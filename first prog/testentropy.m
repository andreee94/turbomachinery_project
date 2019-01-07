
myN = 50;

T = linspace(620, 700+273.15, myN);

p = linspace(10, 150, myN);


[arrayp, arrayT] = meshgrid(p, T);

ss = zeros(myN, myN);

for kk = 1:myN^2
    ss(kk) = XSteam('s_pT', arrayp(kk), arrayT(kk) - 273);
end

figure
surf(arrayp, arrayT, ss)

hold on

plot3([points.p], [points.T], 9* ones(size([points.p])), 'ro', 'LineWidth', 3)

view(0,90)  % XY
% 
% P0_.p = 150;
% P0_.T = 700+273.15;
% P0_.h = XSteam('h_pT', P0_.p, P0_.T - 273.15);
% P0_.s = XSteam('s_ph', P0_.p, P0_.h);
% 
% P2_.p = 42;%P0_.p / 2.5;
% P2_.T = P0_.T * (P2_.p / P0_.p)^0.20;
% %P2_.T = XSteam('T_ps', P2_.p, P0_.s);
% P2_.hIS = XSteam('h_ps', P2_.p, P0_.s);
% 
% plot3([P0_.p P2_.p], [P0_.T P2_.T], 9* ones(2,1), 'k-', 'LineWidth', 3)
% 
% P2_.T = XSteam('T_ph', P2_.p, P2_.hIS) + 273.15;

% k = log(P2_.T/P0_.T) / log(P2_.p/P0_.p);
% 
% plot3([P0_.p P2_.p], [P0_.T P2_.T], 9* ones(2,1), 'b-', 'LineWidth', 3)
% 
% P2_.s = XSteam('s_pT', P2_.p, P2_.T - 273.15);
% P2_.s = XSteam('s_ph', P2_.p, P2_.hIS);
% 
% 
% P2_.h = P0_.h - 0.92 * (P0_.h-P2_.hIS);
% 
% P2_.T = XSteam('T_ph', P2_.p, P2_.h) + 273.15;
% 
% k = log(P2_.T/P0_.T) / log(P2_.p/P0_.p);
% 
% 
% plot3([P0_.p P2_.p], [P0_.T P2_.T], 9* ones(2,1), 'r-', 'LineWidth', 3)

%% entropy calculation on a stator with  Yp = 0
% 
% gamma = 1.33;
% cp = gl.cp0;
% R = gl.R;%gl.cp0 - gl.cv0;%gl.R;
% %cp = R * gamma / (gamma-1);
% 
% 
% v1 = result.triangle.v1;
% v1q = v1^2 / 2;
% v0q = result.triangle.v1A^2 / 2;
% hT1 = cp * gl.TT0;%result.points(1).hT;
% TT1 = gl.TT0;%result.points(1).hT;
% 
% cp/(cp-R)
% 
% deltas = @(Y) cp * log((hT1 - v1q)/(hT1 - v0q)) - R * log(Y * (2*R/cp*(hT1-v0q) + 2*v0q) / (2*R/cp*(hT1-v1q) + 2*v1q) * (hT1 - v1q)/(hT1 - v0q));
% 
% figure
% plot(0.9:0.0001:1, deltas(0.9:0.0001:1));
% grid on









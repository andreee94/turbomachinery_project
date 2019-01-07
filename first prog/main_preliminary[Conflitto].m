
test.beta_from_n = false;
test.b_from_Dm = false;
test.optimal_velocity_triangle = false;

gl.mdot = 200;
gl.betaTT = (2.5)^(1/1);
gl.TT0 = 700+273.15;
gl.pT0 = 150;
gl.MM = 2 * 1 + 16;
gl.R = 8314 / gl.MM;
gl.toll = 1e-6;

% supposing n = 3000rpm
% non si può fare perché il diametro è di 5 metri
% n = 3000;
% omega = 2 * pi * n / 60;

% steam approximation
% gamma = 1.33;
gl.cv0 = XSteam('Cv_pT', gl.pT0, gl.TT0 - 273.15) * 1000; % 1000 to convert in joule
gl.cp0 = XSteam('Cp_pT', gl.pT0, gl.TT0 - 273.15) * 1000; % 1000 to convert in joule
gl.gamma0 = gl.cp0 / gl.cv0;
% k = (gamma - 1)/gamma
gl.k0 = 1 - gl.gamma0^-1;
gl.rho0T = XSteam('rho_pT', gl.pT0, gl.TT0 - 273.15);
%hT0 = cp0 * gl.TT0;

%first guess
%every quantity is at MID
%lambda = 1.8;
%phi = 0.6;
%etaTT = 0.92;

%workIS = cp0 * gl.TT0 * (gl.betaTT ^ -k0 - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test n respect to work
if test.beta_from_n
    n_array = [1500:100:15000];
    beta_res = zeros(size(n_array));
    u_res = zeros(size(n_array));
    j = 1;
    for n = n_array
        [beta_res(j), u_res(j)] = beta_from_n(n, 1.5, gl.TT0, gl.cp0, gl.k0, etaTT);
        j = j + 1;
    end
    
    figure
    grid on
    plot(n_array, beta_res)
    figure
    grid on
    plot(n_array, u_res)
    
    % I take 6000 rpm as a good value
    % u = 471
    % betaTT = 1.419
    % 3 ideal stages we geta betaTT = 2.74
end

% n = 6000;
% betaTT = gl.betaTT ^ (1 / 3);
%
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% work = workIS * etaTT;
% %lamda = work / u^2
% u = sqrt(abs(work) / lambda);
% %work = u * deltaVT
% deltaVT = work / u;
% %phi = VA / u
% VA = u * phi;
% %deltaVT assuming chi = 0.5 is
% % -deltaVT = u + 2 * VT2 so
% VT2 = (-deltaVT - u) * 0.5;
% %deltaVT = VT2 - VT1 so
% VT1 = VT2 - deltaVT;
% %Periferal mach number
% Mu = u / sqrt(gamma0 * R * gl.TT0);
% % dMID from u and omega
% %D = 2 * u / omega;

mdot_f = @(rho, b_su_Dm, n, phi, u, epsilon) pi * rho * b_su_Dm * 60^2 / pi^2 / n^2 * phi * u^3 * epsilon;

u_f = @(mdot, rho, b_su_Dm, n, phi, epsilon) (mdot / (pi * rho * b_su_Dm * 60^2 / pi^2 / n^2 * phi * epsilon)) ^ (1/3);

beta_f = @(lambda, u, eta, cp, TT0, k) (1-lambda * u^2 / (eta * cp * TT0)) ^(-1 / k);

b_f = @(mdot, rho, Dm, n, vA, epsilon) mdot / (rho * Dm * pi * vA * epsilon);

%-------------------------------------------------------
%TODO must iterate


%[smith92_X, smith92_Y, smith92_POLY] = load_csv_poly('graph\smith_92.csv');

smitharray_X = linspace(0.4, 1, 50);
smitharray_Y = linspace(1, 2, 50);

[smitharray_X, smitharray_Y] = meshgrid(smitharray_X, smitharray_Y);

res = [];
res(1).phi = 0;
res(1).lambda = 0;
res(1).uold =0;
res(1).nstages = 0;
res(1).mdot = 0;
res(1).nstagesround = 0;
res(1).beta = 0;
res(1).u = 0;
res(1).vA = 0;
res(1).alpha1 = 0;
res(1).alpha2 = 0;
res(1).deltaBeta = 0;
res(1).b_su_Dm = 0;
res(1).Dm = 0;
res(1).lu1 = 0;
res(1).lu2 =0;
res(1).P0 = [];
res(1).triangle = [];

res = repmat(res(1), size(smitharray_X));

n = 12000;
epsilon = 1;
etaTT = 0.92;

parfor kk = 1:length(smitharray_X(:))
    res(kk) = solve_triangles_and_point0(gl, test, smitharray_X(kk), smitharray_Y(kk), etaTT, epsilon, n);
end

%%

surf(smitharray_X, smitharray_Y, vec2mat([res.nstagesround], length(smitharray_X))')

[smith92_X, smith92_Y, smith92_POLY] = load_csv_poly('graph\smith_92.csv');
[smith94_X, smith94_Y, smith94_POLY] = load_csv_poly('graph\smith_94.csv');

hold on

hh(2) = plot3(smith92_X, smith92_Y, 8*ones(size(smith92_X)),'r', 'LineWidth', 3);
hh(3) = plot3(smith94_X, smith94_Y, 8*ones(size(smith94_X)),'r', 'LineWidth', 3);

uistack(hh(2),'top')
uistack(hh(3),'top')


%%

result = solve_triangles_and_point0(gl, test, 0.6, 1.3, 0.94, 1, 12000);
plot_velocity_triangle(result.triangle);





% P0.TT = gl.TT0;
% P0.pT = gl.pT0;
% P0.v = triangle.v2;
% [P0.T, P0.h] = static_temperature(P0.TT, P0.v, P0.pT, cp0);% it is an approximation since we need also p0 to get cp0
%
% P0.p = P0.pT / (P0.TT / P0.T) ^ k0;
% P0


%     %b = @() gl.mdot ./ (pi .* Dm * rho0TT * phi * pi * n / 60 .* Dm);
%     %b = gl.mdot / (pi * 1.5 * VA * rho0TT)
%     %b_su_D = b ./ 1.5;
%     figure
%     grid on
%     title('b/Dm respect to Dm')
%     Dm_array = [0.5:0.1:3];
%     plot(Dm_array, b(Dm_array) ./ Dm_array)
%
%     hold on
%     plot(Dm_array, 0.025 * ones(size(Dm_array)))

%%








% VA1 = VA;
% VA2 = VA;
% V1 = norm([VT1, VA1]);
% V2 = norm([VT2, VA2]);
% hT1= hT0;
% h1 = hT1 - 0.5 * V1^2;
% %we approximate using cp0 instead of 1 because we need p1 other than h1
% %to calculate T1
% cp1 = cp0;
% gamma1 = gamma0;
% T1 = h1 / cp1;
% M1 = V1 / sqrt(gamma1 * R * T1);
% alpha1 = rad2deg(atan(VT1 / VA1));
% alpha2 = rad2deg(atan(VT2 / VA2));
% alpha0 = alpha2;
%
% %balje chart analysis
%
% V0 = V2;
% T0 = gl.TT0 - V0^2 / 2 /cp0;
% p0 = gl.pT0 / (gl.TT0/T0) ^k0;
% rho0 = XSteam('rho_pT', p0, T0-273.15);
% Q0 = gl.mdot / rho0;
%
% ns = 1;%0.6;
% ws = 2 * pi * ns / 60;
% ds = 3;%4;
%
% w = ws * (-workIS)^0.75 / sqrt(Q0);
% wrpn = w * 60 / 2 /pi;
% dTIP = ds / (-workIS)^0.25 * sqrt(Q0);
%
% % mdot = b * Dm * pi * VA * rho
%
% P0.T = T0;
% P0.TT = gl.TT0;
% P0.V = V0;
% P0.p = p0;
% P0.pT = gl.pT0;
% P0.rho = rho0;
% P0.alpha = alpha0;



% print(work)
% print(u)
% print(Mu)
% print(deltaVT)
% print(VA)
% print(VT1)
% print(V1)
% print(VT2)
% print(V2)
% print(hT1)
% print(h1)
% print(T1)
% print(M1)
% print(alpha1)
% print(alpha2)
% print(w)
% print(dTIP)
% %print(D)
% print(wrpn)

test.beta_from_n = false;
test.b_from_Dm = false;
test.optimal_velocity_triangle = false;
test.hainley_mathieson_plot = false;

gl.mdot = 200;
gl.betaTT = (2.5)^(1/1);
gl.TT0 = 700+273.15;
gl.pT0 = 150;
gl.MM = 2 * 1 + 16;
gl.R = 8314 / gl.MM;
gl.toll = 1e-6;
gl.maxiter = 20;
gl.secant_entropy_delta0 = 0.01;% distance between x0 and x1 in secants algorithm
gl.NGauss = 3;
gl.velocity_evolution = @velocity_evolution_ConstantAngle;
[gl.XGauss, gl.WGauss] = lgwt(gl.NGauss, 0, 1);

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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%[smith92_X, smith92_Y, smith92_POLY] = load_csv_poly('graph\smith_92.csv');

smitharray_X = linspace(0.4, 1, 100);
smitharray_Y = linspace(1, 2, 100);

[smitharray_X, smitharray_Y] = meshgrid(smitharray_X, smitharray_Y);

res = [];
res(1).n = 0;
res(1).phi = 0;
res(1).lambda = 0;
res(1).uold =0;
res(1).nstages = 0;
res(1).mdot = 0;
res(1).YTotRotor = 0;
res(1).YTotStator = 0;
res(1).nstagesround = 0;
res(1).beta = 0;
res(1).u = 0;
res(1).vA = 0;
res(1).alpha1 = 0;
res(1).alpha2 = 0;
res(1).deltaBeta = 0;
res(1).b_su_Dm = 0;
res(1).Dm = 0;
res(1).b = 0;
res(1).lu1 = 0;
res(1).lu2 =0;
%res(1).P0 = [];
%res(1).P1 = [];
%res(1).P2 = [];
res(1).points = [];
res(1).eta = 0;
res(1).etaTT = 0;
res(1).etas = [];
res(1).chi_h = [];
res(1).chi_tr = [];
res(1).triangle = [];
res(1).betas =[];
res(1).betasTT = [];
res(1).betaRes = 0;
res(1).betaResTT = 0;
res(1).ii = [];
res(1).iii = [];

res = repmat(res(1), size(smitharray_X));

n = 6000;
epsilon = 1;
etaTT = 0.92;
chi = 0.47;
chord_su_b = 0.5;
sigma = 1.25;
tmax_su_c = 0.2;

%parpool('local', 4)
tic
parfor kk = 1:length(smitharray_X(:))
    res(kk) = solve_triangles_and_point0(gl, test, smitharray_X(kk), smitharray_Y(kk), etaTT, epsilon, n, chi, chord_su_b, sigma, tmax_su_c);
    print(kk);
end

t = toc;
disp(['Time for ' num2str(length(smitharray_X(:))) ' iterations = ', num2str(t)]);

%%

% pT_stages = zeros(size(smitharray_X(:)));
%
% for kk = 1:length(smitharray_X(:))
%     tempP0 = res(kk).P0;
%     tempP2 = res(kk).P2;
%     if isfield(tempP0,'p')
%         pT_stages(kk) = tempP2.p / tempP0.p;
%     end
% end

[smith92_X, smith92_Y, smith92_POLY] = load_csv_poly('graph\smith_92.csv');
[smith94_X, smith94_Y, smith94_POLY] = load_csv_poly('graph\smith_94.csv');

figure
%etas = [res(:).eta];
%surf(smitharray_X, smitharray_Y, vec2mat([etas.stageSS], length(smitharray_X)))
%surf(smitharray_X, smitharray_Y, vec2mat([res.nstages], length(smitharray_X))')
%hold on
%contourf(smitharray_X, smitharray_Y, vec2mat([res.nstagesround] .* [res.etaTT], length(smitharray_X))',100)

surf(smitharray_X, smitharray_Y, vec2mat([res.nstagesround], length(smitharray_X))')

hold on

hh(2) = plot3(smith92_X, smith92_Y, 10*ones(size(smith92_X)),'r', 'LineWidth', 3);
hh(3) = plot3(smith94_X, smith94_Y, 10*ones(size(smith94_X)),'r', 'LineWidth', 3);

uistack(hh(2), 'top')
uistack(hh(3), 'top')

view(0,90)  % XY
%caxis([8 9])

figure
surf(smitharray_X, smitharray_Y, vec2mat([res.etaTT].*[res.nstagesround], length(smitharray_X))')
hold on
hh(2) = plot3(smith92_X, smith92_Y, 10*ones(size(smith92_X)),'r', 'LineWidth', 3);
hh(3) = plot3(smith94_X, smith94_Y, 10*ones(size(smith94_X)),'r', 'LineWidth', 3);
view(0,90)  % XY
caxis([3 4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test on best chi

chord_su_b = 0.5;
sigma = 1.25;
tmax_su_c = 0.2;
phi = 0.6;
lambda = 1.3;
chi = 0.5;

best_chi = chi;
best_sigma = sigma;
best_chord_su_b = chord_su_b;
%%

chi_array = linspace(0.45, 0.55, 100);
etaTT_res = zeros(size(chi_array));
nstage_res = zeros(size(chi_array));

parfor kk = 1:length(chi_array)
    
    chi = chi_array(kk);
    
    result = solve_triangles_and_point0(gl, test, phi, lambda, 0.90, 1, 6000, chi, chord_su_b, sigma, tmax_su_c);
    etaTT_res(kk) = result.etaTT;
    nstage_res(kk) = result.nstagesround;
    
    print(kk)
end

figure
yyaxis left
[~, best_chi] = plotMinMax(chi_array, etaTT_res, true, true);
%yyaxis right
%[~, ~] = plotMinMax(chi_array, nstage_res, true, true);
title('eta_T_T function of chi')


%% test on best sigma

sigma_array = linspace(0.8, 2, 100);
etaTT_res = zeros(size(sigma_array));
kk = 1;

parfor  kk = 1:length(sigma_array)
    
    ss = sigma_array(kk);
    
    result = solve_triangles_and_point0(gl, test, phi, lambda, 0.90, 1, 6000, best_chi, chord_su_b, ss, tmax_su_c);
    etaTT_res(kk) = result.etaTT;
    
    print(kk)
end

figure
[~, best_sigma] = plotMinMax(sigma_array, etaTT_res, true, true);
title('eta_T_T function of sigma')

%% test on best chord/b

chord_su_b_array = linspace(0.25, 5, 100);
etaTT_res = zeros(size(chord_su_b_array));
kk = 1;

parfor  kk = 1:length(chord_su_b_array)
    
    cc = chord_su_b_array(kk);
    
    result = solve_triangles_and_point0(gl, test, phi, lambda, 0.90, 1, 6000, best_chi, cc, best_sigma, tmax_su_c);
    etaTT_res(kk) = result.etaTT;
    
    print(kk)
end

figure
[~, best_chord_su_b] = plotMinMax(chord_su_b_array, etaTT_res, true, true);
title('eta_T_T function of chord/b')


%% test on best n

n_array = linspace(1500, 25000, 101);
etaTT_res = zeros(size(n_array));
nstages_res = zeros(size(n_array));
nstagesround_res = zeros(size(n_array));

parfor  kk = 1:length(n_array)
    
    nn = n_array(kk);
    
    result = solve_triangles_and_point0(gl, test,  phi, lambda, 0.90, 1, nn, best_chi, chord_su_b, best_sigma, tmax_su_c);
    etaTT_res(kk) = result.etaTT;
    nstages_res(kk) = result.nstages;
    nstagesround_res(kk) = ceil(nstages_res(kk));
    
    print(kk)
end

figure
yyaxis left
[~, best_n] = plotMinMax(n_array, etaTT_res, true, true);
yyaxis right
plotMinMax(n_array, nstagesround_res, true, true);

%% test on best lambda

lambda_array = linspace(0.5, 3, 100);
etaTT_res = zeros(size(lambda_array));
nstages_res = zeros(size(lambda_array));
nstagesround_res = zeros(size(lambda_array));

parfor  kk = 1:length(lambda_array)
    
    %ss = sigma_array(kk)
    ll = lambda_array(kk);
    
    result = solve_triangles_and_point0(gl, test, phi, ll, 0.90, 1, 6000, best_chi, chord_su_b, best_sigma, tmax_su_c);
    etaTT_res(kk) = result.etaTT;
    nstages_res(kk) = result.nstages;
    nstagesround_res(kk) = ceil(nstages_res(kk));
    
    print(kk)
end
figure
yyaxis left
[~, best_lambda] = plotMinMax(lambda_array, etaTT_res, true, true);
yyaxis right
plotMinMax(lambda_array, nstagesround_res, true, true);
title('eta_T_T function of lambda')


%% test on best phi

phi_array = linspace(0.3, 2, 100);
etaTT_res = zeros(size(phi_array));
nstages_res = zeros(size(phi_array));
nstagesround_res = zeros(size(phi_array));

parfor  kk = 1:length(phi_array)
    
    %ss = sigma_array(kk)
    pp = phi_array(kk);
    
    result = solve_triangles_and_point0(gl, test, pp, lambda, 0.90, 1, 6000, best_chi, chord_su_b, best_sigma, tmax_su_c);
    etaTT_res(kk) = result.etaTT;
    nstages_res(kk) = result.nstages;
    nstagesround_res(kk) = ceil(nstages_res(kk));
    
    print(kk)
end

figure
yyaxis left
[~, best_phi] = plotMinMax(phi_array, etaTT_res, true, true);
yyaxis right
plotMinMax(phi_array, nstagesround_res, true, true);
title('eta_T_T function of phi')


%% test on N gauss

ngauss_array = 3:12;
b_res = zeros(size(ngauss_array));
kk = 1;

parfor  kk = 1:length(ngauss_array)
    
    gl_new = gl;
    gl_new.NGauss = ngauss_array(kk);
    [gl_new.XGauss, gl_new.WGauss] = lgwt(gl_new.NGauss, 0, 1);
    
    result = solve_triangles_and_point0(gl_new, test, phi, lambda, 0.90, 1, 6000, best_chi, chord_su_b, best_sigma, tmax_su_c);
    b_res(kk) = result.points(end).b;%sum([result.points.b]);
    print(kk)
end

figure
semilogy(ngauss_array, 100*b_res / b_res(1) - 100);
%[~, ~] = plotMinMax(ngauss_array, b_res / b_res(1) - 1, true, true);
title('guass points respect to sum of b of all points')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% solve_triangles_and_point0(gl, test, phi, lambda, etaTT, epsilon, n, chiDEF, chord_su_b, sigma, tmax_su_c )
tic

% confronto con documento trovato
gl.mdot = 22.9036;
gl.betaTT = 1.69;
gl.pT0 = 150;
gl.TT0 = 1114.15 -273.15;

result = solve_triangles_and_point0(gl, test, 1 , 2.8, 0.92, 1, 7500, 0.47, 1, sigma, tmax_su_c);

%result = solve_triangles_and_point0(gl, test, 0.4 , 1.6, 0.85, 1, 7500, 0.47, 1, 1.25, 0.2); 

t = toc;
disp(['Time for iteration = ', num2str(t)]);
disp([result.points.b])

points = result.points;%[result.P0 result.P1 result.P2];

figure
subplot(2, 2, 1)
plot([0:result.nstagesround*2], [points.h]);
legend('static enthalpy h')
subplot(2, 2, 2)
plot([0:result.nstagesround*2], [points.s]);
legend('static entropy s')
subplot(2, 2, 3)
plot([0:result.nstagesround*2], [points.p]);
legend('static pressure p')
subplot(2, 2, 4)
plot([0:result.nstagesround*2], [points.T]);
legend('static temperature T')

testentropy

plot_velocity_triangle(result.triangle);

plot_machine(result, result.points);

plot_pTrho_machine(result);

kk = 2;
bb = linspace(-0.5 * result.points(kk).b, 0.5 * result.points(kk).b, 50);
%triangles =  result.span.triangle(bb);
%figure
%hold on
%plot(bb, triangles.v1A)
%plot(bb, triangles.alpha1_deg)
%plot(bb, triangles.alpha2_deg)

%%
figure
hold on
for tmax_su_c = [0.2:0.2:1]
    
    [Ytot] = AinleyMathiesonLosses(test, result.triangle.alpha2_deg, result.triangle.alpha1_deg, 0.2:0.01:2, tmax_su_c, result.points(1).rho, result.P0.mu, result.b, result.triangle.v1, 0, result.b);
    
    plot([0.2:0.01:2], Ytot)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% solidity rotor and stator
figure
hold on
sigma = linspace(1.2, 1.3, 100);
for tmax_su_c = [0.2:0.2:1]
    
    [YtotStat] = AinleyMathiesonLosses(test, 0, result.triangle.alpha1_deg, sigma, tmax_su_c, result.P0.rho, result.P0.mu, result.b, result.vA, 0, result.b);
    [YtotRot] = AinleyMathiesonLosses(test, result.triangle.beta1_deg, result.triangle.beta2_deg, sigma, tmax_su_c, result.P0.rho, result.P0.mu, result.b, result.triangle.w1, result.b/100, result.b);
    
    %test.hainley_mathieson_plot = true;
    %[YtotRot] = AinleyMathiesonLosses(test, result.triangle.beta1_deg, result.triangle.beta2_deg, sigma, tmax_su_c, result.P0.rho, result.P0.mu, result.b, result.triangle.w1, 0, result.b);
    
    %test.hainley_mathieson_plot = false;
    
    %figure
    plotMinMax(sigma, YtotStat, true, true)
    hold on
    plotMinMax(sigma, YtotRot, true, true)
    hold on
    %plot(sigma, YtotRot)
end


%%
figure
hold on
c = linspace(0.5, 2, 100) * result.b;

for tmax_su_c = [0.2:0.2:1]
    
    [Ytot] = AinleyMathiesonLosses(test, result.triangle.alpha2_deg, result.triangle.alpha1_deg, 1, tmax_su_c, result.P0.rho, result.P0.mu, c, result.triangle.v1, 0, result.b);
    
    
    plot(c /  result.b, Ytot)
end

%P1 = solve_stator(gl, result.P0, result.triangle, result.YTotStator, result.beta);


%plot( linspace(0.2, 1.8, 100), Ytot)


%%
vMeanStat = 0.5 * (result.vA + result.triangle.v1);
[YTotStator] = AinleyMathiesonLosses(test, 0, result.triangle.alpha1_deg, 1, 0.2, result.P0.rho, result.P0.mu, 0.5*result.b, vMeanStat, 0, result.b);

%statorP1 = solve_stator(gl, result.P0, result.triangle, 0, result.beta)

%eta_statorSS = (result.P0.h - statorP1.h) / (result.P0.h - statorP1.hIS)
%eta_statorTT = (result.P0.hT - statorP1.hT) / (result.P0.hT - statorP1.hIS + 0.5*result.vA^2 - 0.5*result.triangle.v1^2)

%%






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
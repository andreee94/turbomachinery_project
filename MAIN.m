


test.beta_from_n = false;
test.b_from_Dm = false;
test.optimal_velocity_triangle = false;
test.hainley_mathieson_plot = false;
test.print = true;

% initial data
gl.mdot = 200;
gl.betaTT = 2.5;
gl.TT0 = 700+273.15;
gl.pT0 = 150;

% perfect gas data
gl.MM = 2 * 1 + 16;
gl.R = 8314 / gl.MM;

% mathematical data
gl.toll = 1e-6;
gl.maxiter = 50;
gl.secant_entropy_delta0 = 0.01;% distance between x0 and x1 in secants algorithm
gl.NGauss = 10;
[gl.XGauss, gl.WGauss] = lgwt(gl.NGauss, 0, 1);

% velocity triangle span evolution
gl.velocity_evolution = @velocity_evolution_ConstantAngle;

% design parameters
gl.chiMID = 0.4652;
gl.b_su_DmMIN = 0.05;
gl.n = 6000;
gl.usePhi = true;
gl.phi = 0.4465; % vA / u
gl.v1t_by_kp_su_u = 1; % kp = u / v1 = coef * u / v1t % This is an alternative to define phi
gl.autoNStages = false;
gl.nStages = 4;
gl.etaGuess = 0.93;
gl.epsilon = 1;
gl.chord_su_b = .5;
gl.solidity = 1.219;
gl.clearance = 0.02*0.025; % 0.02 * approx b
gl.seals_num = 2;
gl.tmax_su_c = 0.2;

res = solve_turbine(gl, test);
res = machanicalAnalysis(res);

% gl.chiMID = 0.2;
% res2 = solve_turbine(gl, test);
%
% gl.chiMID = 0.7;
% res3 = solve_turbine(gl, test);

%%

stagesArray = 1:1:8;
stagesArrayFINE = 1:1:14;
solidityArray = linspace(1, 1.4, 10);
solidityArrayFINE = linspace(1, 1.4, 100);
phiArray = linspace(0.3, 0.9, 7);
phiArrayFINE = linspace(0.3, 1.5, 50);
chiArray = linspace(0.35, 0.55, 5);
chiArrayFINE = linspace(0.35, 0.55, 100);
chordsubArray = linspace(0.25, 2, 8);
bsudArray = linspace(0.025, 0.5, 10);
clearanceArray = linspace(1, 10, 10) * 1e-4;
NGaussArray = 1:10;
sealsArray = 1:10;
nArray = 3000:1500:12000;%linspace(3000, 12000, 7);
nArrayFINE = 3000:250:12000;%linspace(3000, 12000, 7);

%% real case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solidity and phi
figure
etas = solveAndPlot(gl, test, phiArray, solidityArray, 'phi', 'solidity', 'etaTT', false, true);

% solidity and rpm
figure
etas = solveAndPlot(gl, test, nArray, solidityArray, 'n', 'solidity', 'etaTT', false, true);

% solidity and nstages
figure
etas = solveAndPlot(gl, test, stagesArray, solidityArrayFINE, 'nStages', 'solidity', 'etaTT', false, true);

% solidity and chi mid
figure
etas = solveAndPlot(gl, test, chiArray, solidityArray, 'chiMID', 'solidity', 'etaTT', false, true);

% solidity and chord su b
figure
etas = solveAndPlot(gl, test, chordsubArray, solidityArray, 'chord_su_b', 'solidity', 'etaTT', false, true);

% solidity and number of stages wrong interval
figure
etas = solveAndPlot(gl, test, stagesArray, linspace(0.55, 3, 100), 'nStages', 'solidity', 'etaTT', false, true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chi and phi
figure
etas = solveAndPlot(gl, test, phiArray, chiArray, 'phi', 'chiMID', 'etaTT', false, true);

% chi and rpm
figure
etas = solveAndPlot(gl, test, nArray, chiArray, 'n', 'chiMID', 'etaTT', false, true, false, 'reverse');

% chi and nstages
figure
etas = solveAndPlot(gl, test, stagesArray, chiArray, 'nStages', 'chiMID', 'etaTT', false, true, false, 'reverse');

% solidity and chimid
figure
etas = solveAndPlot(gl, test, gl.solidity, chiArrayFINE, 'solidity', 'chiMID', 'etaTT', false, true);

% chi and chord su b
figure
etas = solveAndPlot(gl, test, chordsubArray, chiArray, 'chord_su_b', 'chiMID', 'etaTT', false, true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seals_number and rpm
figure
etas = solveAndPlot(gl, test, sealsArray, nArray, 'seals_num', 'n', 'etaTT', false, true, false, 'reverse');
% seals_number and nstages
figure
etas = solveAndPlot(gl, test, sealsArray, stagesArray, 'seals_num', 'nStages', 'etaTT', false, true, false, 'reverse');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearance and rpm
figure
etas = solveAndPlot(gl, test, clearanceArray, nArray, 'clearance', 'n', 'etaTT', false, true, false, 'reverse');
% clearance and nstages
figure
etas = solveAndPlot(gl, test, clearanceArray, stagesArray, 'clearance', 'nStages', 'etaTT', false, true, false, 'reverse');
% clearance and seals_number
figure
etas = solveAndPlot(gl, test, clearanceArray, sealsArray, 'clearance', 'seals_num', 'etaTT', false, true, false, 'reverse');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% rpm and number of stages
figure
etas = solveAndPlot(gl, test, nArray, stagesArrayFINE, 'n', 'nStages', 'etaTT', false, true);
figure
etas = solveAndPlot(gl, test, stagesArray, nArrayFINE, 'nStages', 'n', 'etaTT', false, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi and rpm
figure
etas = solveAndPlot(gl, test, [4000:500:7000], phiArrayFINE, 'n', 'phi',  'etaTT', false, true);

% phi and number of stages
figure
etas = solveAndPlot(gl, test, [3:5], phiArrayFINE,  'nStages', 'phi', 'etaTT', false, true);

% solidity
figure
etas = solveAndPlot(gl, test,gl.n, linspace(1.2, 1.3, 100), 'n', 'solidity', 'etaTT', false, true);

% chi mid
figure
etas = solveAndPlot(gl, test,gl.n, linspace(0.45, 0.5, 100), 'n', 'chiMID', 'etaTT', false, true);





%%

pause(0.1)
figure
sumb = solveAndPlot(gl, test, stagesArray, NGaussArray, 'nStages', 'NGauss' ,  'sumb', false, false, true);

pause(0.1)
figure
workerrorperc = solveAndPlot(gl, test,nArray, stagesArray,'n' , 'nStages',   'workerrorperc', false, true);

pause(0.1)
figure
workerrorperc = solveAndPlot(gl, test, bsudArray, stagesArray, 'b_su_DmMIN' , 'nStages', 'beta_real', true, true);

pause(0.1)
figure
workerrorperc = solveAndPlot(gl, test, nArray, stagesArray, 'n' , 'nStages', 'beta_real', true, true);

pause(0.1)
figure
etas = solveAndPlot(gl, test, bsudArray, stagesArray, 'b_su_DmMIN' , 'nStages',  'etaTT', false, true);

pause(0.1)
figure
etas = solveAndPlot(gl, test, bsudArray, nArray, 'b_su_DmMIN' , 'n',  'etaTT', false, true);
%%

pause(0.1)
figure
etas = solveAndPlot(gl, test, chordsubArray, stagesArray, 'chord_su_b' , 'nStages',  'etaTT', false, true);
%etas = solveAndPlot(gl, test, stagesArray, chordsubArray, 'nStages', 'chord_su_b', 'etaTT', false, true);

pause(0.1)
figure
etas = solveAndPlot(gl, test, chiArray, stagesArray, 'chiMID' , 'nStages',  'etaTT', false, true);

pause(0.1)
figure
etas = solveAndPlot(gl, test, phiArray, stagesArray, 'phi' , 'nStages',  'etaTT', false, true);

pause(0.1)
figure
etas = solveAndPlot(gl, test, nArray, stagesArray, 'n' , 'nStages',  'etaTT', false, true);

%%

figure
etas = solveAndPlot(gl, test, nArray, chiArray, 'n', 'chiMID', 'etaTT', false, true);

pause(0.1)
figure
etas = solveAndPlot(gl, test, stagesArray, nArray, 'nStages', 'n', 'etaTT', false, true);

pause(0.1)
figure
etas = solveAndPlot(gl, test, sealsArray, nArray, 'seals_num', 'n', 'etaTT', false, true);

pause(0.1)
figure
etas = solveAndPlot(gl, test, stagesArray, solidityArray, 'nStages', 'solidity', 'etaTT', false, true);

%figure
%etas = solveAndPlot(gl, test, stagesArray, solidityArray, 'nStages', 'solidity', 'etaTT', false, true);

%plotMinMax(stagesArray, stagesEta, false, true)

%%

%
% function [etas] = solveAndPlot(gl, test, array1, array2, prop1, prop2, propOut, plotmin, plotmax)
%
%     if ~isempty(array2)
%         etas = zeros(length(array1), length(array2));
%         legendString = cell(length(array1), 1);
%         for ii = 1:length(array1)
%             gl = setfield(gl, prop1, array1(ii));
%             tic
%             for kk = 1:length(array2)
%                 %kk
%                 tempgl = gl;
%                 tempgl = setfield(tempgl, prop2, array2(kk));
%
%                 res = solve_turbine(tempgl, test);
%                 etas(ii, kk) = getfield(res, propOut);
%             end
%             toc
%             disp(['Progress of outer loop = ', num2str(ii / length(array1) * 100), '%'])
%             plotMinMax(array2, etas(ii, :), plotmin, plotmax)
%             legendString{ii} = [prop1, ' = ', num2str(array1(ii))];
%         end
%
%         xlabel(prop2)
%         ylabel(propOut)
%         f=get(gca,'Children');
%         step = 1+plotmin+plotmax;
%         f = f(end:-step:1);
%         %f(1+plotmin+plotmax:1+plotmin+plotmax:end);
%         legend(f,legendString)
%     end
%
% end


%% plot blade

P = generateBlade( [], res.MID.velocity_triangle.beta1, res.MID.velocity_triangle.beta2, res.points(end).c, gl.tmax_su_c, res.stages(end).rotorSolidity );

figure
grid on
hold on
axis equal
plot(P.xUP, -P.yUP, 'k', 'LineWidth', 2);
plot(P.xDOWN, -P.yDOWN, 'k', 'LineWidth', 2);
plot(P.x, -P.yC, '--k', 'LineWidth', 1);

%% 

colors = get(gca,'colororder');

% stator losses
figure
grid on
y = [ [res.stages.statorYProfile]', [res.stages.statorYSec]', [res.stages.statorYClearance]'] * 100;
b = bar(y,'stacked', 'FaceColor','flat');
b(1).Parent.Parent.Colormap =  colors(1:3,:);
%b.CData(:,:) = colors(1:3,:);
legend('Profile losses', 'Secondary losses', 'Clearance losses')

% rotor losses
figure
grid on
y = [ [res.stages.rotorYProfile]', [res.stages.rotorYSec]', [res.stages.rotorYClearance]'] * 100;
b = bar(y,'stacked', 'FaceColor','flat');
b(1).Parent.Parent.Colormap =  colors(1:3,:);
%b.CData(:,:) = colors(1:3,:);
legend('Profile losses', 'Secondary losses', 'Clearance losses')

%%
figure
grid on
hold on
spans = [res.stages.span];
y = [spans.statorYprofile]' * 100;

for ii = 1:length(res.stages)
    plot(gl.XGauss, y(ii, :), 'LineWidth', 1.5);
end
xlabel('Span coordinate')
ylabel('Profile losses')
legend('Stator 1', 'Stator 2', 'Stator 3', 'Stator 4')


figure
grid on
hold on
spans = [res.stages.span];
y = [spans.rotorYprofile]' * 100;

for ii = 1:length(res.stages)
    plot(gl.XGauss, y(ii, :), 'LineWidth', 1.5);
end
xlabel('Span coordinate')
ylabel('Profile losses')
legend('Rotor 1', 'Rotor 2', 'Rotor 3', 'Rotor 4')






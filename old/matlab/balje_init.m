
global balje_data_80_X balje_data_80_Y
global balje_data_75_X balje_data_75_Y
global balje_data_70_X balje_data_70_Y
global balje_data_60_X balje_data_60_Y
global balje_data_50_X balje_data_50_Y
global balje_data_axial_X balje_data_axial_Y
global balje_data_axial_X_false balje_data_axial_Y_false

balje_data_80 = csvread('data\balje_80.csv');
balje_data_75 = csvread('data\balje_75.csv');
balje_data_70 = csvread('data\balje_70.csv');
balje_data_60 = csvread('data\balje_60.csv');
balje_data_50 = csvread('data\balje_50.csv');
balje_data_axial = csvread('data\balje_axial.csv');
balje_data_axial_false = csvread('data\balje_axial_limit.csv');

balje_data_80_X = balje_data_80(:,1);
balje_data_80_Y = balje_data_80(:,2);

balje_data_75_X = balje_data_75(:,1);
balje_data_75_Y = balje_data_75(:,2);

balje_data_70_X = balje_data_70(:,1);
balje_data_70_Y = balje_data_70(:,2);

balje_data_60_X = balje_data_60(:,1);
balje_data_60_Y = balje_data_60(:,2);

balje_data_50_X = balje_data_50(:,1);
balje_data_50_Y = balje_data_50(:,2);

balje_data_axial_X = balje_data_axial(:,1);
balje_data_axial_Y = balje_data_axial(:,2);

balje_data_axial_X_false = balje_data_axial_false(:,1);
balje_data_axial_Y_false = balje_data_axial_false(:,2);

P = balje_data_70_X(1);
%[balje_data_70_X, balje_data_70_Y] = points2contour(balje_data_70_X, balje_data_70_Y,P,'ccw');
 
% if graph
%     figure(ff); ff = ff+1;
%     loglog(balje_data_80_X,balje_data_80_Y)%,'*');
%     hold on
%     axis equal
%     loglog(balje_data_75_X,balje_data_75_Y)%,'*');
%     loglog(balje_data_70_X,balje_data_70_Y)%,'*');
%     loglog(balje_data_60_X,balje_data_60_Y)%,'*');
%     loglog(balje_data_50_X,balje_data_50_Y)%,'*');
% end


% s_80 = find(balje_data_80_Y(1) < balje_data_axial_Y,1,'first');
% e_80 = find(balje_data_80_Y(end) > balje_data_axial_Y,1,'last');
%
% s_75 = find(balje_data_75_Y(1) < balje_data_axial_Y,1,'first');
% e_75 = find(balje_data_75_Y(end) > balje_data_axial_Y,1,'last');
%
% s_70 = find(balje_data_70_Y(1) < balje_data_axial_Y,1,'first');
% e_70 = find(balje_data_70_Y(end) > balje_data_axial_Y,1,'last');
%
% s_60 = find(balje_data_60_Y(1) < balje_data_axial_Y,1,'first');
% e_60 = find(balje_data_60_Y(end) > balje_data_axial_Y,1,'last');
%
% s_50 = find(balje_data_50_Y(1) < balje_data_axial_Y,1,'first');
% e_50 = find(balje_data_50_Y(end) > balje_data_axial_Y,1,'last');
%
% plot(balje_data_axial_X(s_80:e_80),balje_data_axial_Y(s_80:e_80),'*');
%
% r = linspace(0,20,50);
% for ii = 1:50;
%
%     if inpolygon(r(ii),3.5,balje_data_50_X,balje_data_50_Y);
%         plot(r(ii),3.5,'ko');
%     else
%         plot(r(ii),3.5,'bo');
%     end
% end

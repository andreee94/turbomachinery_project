
%  c = parcluster('local');
%  c.NumWorkers = 16;
%  parpool(c, c.NumWorkers);

global m_P beta_TOT R Cp gamma p0 T0 ws Ds eta l graph v_M

graph = true;
ff = 1;

howell_init
balje_init

DATA.howell_POLY = howell_POLY;

DATA.balje_data_80_X = balje_data_80_X;
DATA.balje_data_80_Y = balje_data_80_Y;
DATA.balje_data_75_X = balje_data_75_X;
DATA.balje_data_75_Y = balje_data_75_Y;
DATA.balje_data_70_X = balje_data_70_X;
DATA.balje_data_70_Y = balje_data_70_Y;
DATA.balje_data_60_X = balje_data_60_X;
DATA.balje_data_60_Y = balje_data_60_Y;
DATA.balje_data_50_X = balje_data_50_X;
DATA.balje_data_50_Y = balje_data_50_Y;

m_P = 65;%kg/s
beta_TOT = 13.5;
R = 287;%K/kgK
gamma = 1.4;
Cp = gamma/(gamma-1)*R;%K/kgK
p0 = 0.85*1e5;%pa
T0 = 268;%k
rho0 = p0/R/T0;

l = @(eta, T_in, beta) Cp*T_in/eta*(beta^((gamma-1)/gamma)-1);

%da balje
% ws_Ds_array_error = [];
% ws_Ds_error_index = 1;

% if exist('data\ws_ds_wrong.mat', 'file') == 2
%     load('data\ws_ds_wrong.mat')
%     ws_Ds_error_index = size(ws_Ds_array_error,1);
% end

% ws_array = linspace(2.3,2.5,100);
% Ds_array = linspace(1.9,2.2,100);

ws_array = linspace(1,5,150);
Ds_array = linspace(1,4,150);

ii = 0;
% p0_ = p0;
% m_P_ = m_P;


item_temp.ws = [];
item_temp.Ds = [];
item_temp.beta_TOT = 0;
item_temp.eta_TOT = 0;
item_temp.eta_TOT_in_out = 0;
item_temp.N_TOT = 0;
item_temp.beta_MEDIO = 0;
item_temp.D = 0;
item_temp.w = 0;
item_temp.rpm = 0;
item_temp.T_OUT = 0;
item_temp.errorU = [];
item_temp.errorBalje = [];
item_temp.bsuD_MIN = 0;
item_temp.bsuD_MAX = 0;

ws_array_length = length(ws_array);
Ds_array_length = length(Ds_array);

iter_res = repmat(item_temp,ws_array_length* Ds_array_length,1);

DATA.m_P = m_P;
DATA.l = l;
DATA.p0 = p0;
DATA.T0 = T0;
DATA.R = R;
DATA.Cp = Cp;
DATA.beta_TOT = beta_TOT;
DATA.gamma = gamma;


eta = 0.8; [ST,~,n_intermedio] = calcolo_stadi_beta_costante(1.2);
DATA.array_of_struct = repmat(ST(1), 20, 1 );

beta_approx_array = linspace(1.1,1.4,300);

tic
parfor kk = 1:ws_array_length*Ds_array_length
    
    %kk
    ii = mod(kk, Ds_array_length);
    if ii == 0
        ii = ws_array_length;
    end
    jj = ceil(kk/Ds_array_length);
    
    %         if ~isempty(ws_Ds_array_error) && ismember([ws, Ds], ws_Ds_array_error, 'rows')
    %             continue;
    %         end
    DATA_iter = struct();
    DATA_iter.ws = ws_array(ii);
    DATA_iter.Ds = Ds_array(jj);
    
    
    item = struct();
    item.ws = DATA_iter.ws;
    item.Ds = DATA_iter.ws;
    
    DATA_iter.eta = balje_PAR(DATA_iter.ws,DATA_iter.ws,DATA);%da balje
    
    if isnan(DATA_iter.eta)
        item.errorBalje = true;
        item.errorU = false;
        continue;
    end
    
    
    %% OTTIMIZZAZIONE BETA CON U MASSIMO POSSIBILE
    
    uu_opt = zeros(300,1);
    jj_index = 1;
    %beta_approx_array = linspace(1.19,1.21,500);
    
    for beta_approx = beta_approx_array;
        
        [ST,~,n_intermedio] = calcolo_stadi_beta_costante_PAR(beta_approx,DATA,DATA_iter);
        
        uu_opt(jj_index) = ST(n_intermedio).U;
        jj_index = jj_index + 1;
        
    end
    %plot(beta_approx_array,uu_opt);
    
    beta_approx = beta_approx_array(find(uu_opt<500,1,'last'));
    
    if isempty(beta_approx) || max(uu_opt) < 500
        item.errorU = true;
        continue
    end
    
    item.errorU = false;
    
    %% CALCOLO STADI CON BETA OTTIMIZZATO
    
    [ST,N,n_intermedio] = calcolo_stadi_beta_costante_PAR(beta_approx,DATA,DATA_iter);
    
    ST_stima_beta_cost = ST;
    
    jj_index = 1;
    v_m = zeros(500,1);
    alpha1_array = linspace(0.1,deg2rad(50),500);
    
    for alpha1 = alpha1_array
        
        beta2 = -alpha1;
        
        deltaBeta = howell(DATA.howell_POLY,beta2);
        %deltaBeta = deg2rad(deltaBeta);
        
        v_m(jj_index) = ST(n_intermedio).l/ST(n_intermedio).U/(abs(tan(abs(beta2)+deltaBeta))-tan(abs(beta2)));
        %v_m(jj) = ST(n_intermedio).l/ST(n_intermedio).U/(tan(beta2+deltaBeta)-tan(beta2));
        jj_index = jj_index +1;
    end
    
    DATA_iter.v_M = max(v_m);%(end)
    
    beta2_deg = rad2deg(-alpha1_array(DATA_iter.v_M == v_m));
    
    [deltabeta, beta1] = howell(DATA.howell_POLY, deg2rad(beta2_deg));
    
    deltabeta_deg = rad2deg(deltabeta);
    beta1_deg = rad2deg(beta1);
    
    ST0 = ST(1);
    ST0.l = ST(n_intermedio).l;
    ST0.U = ST(n_intermedio).U;
    ST0.D = ST(n_intermedio).D;
    ST0.w = ST(n_intermedio).w;
    
    ST = calcolo_stadi_l_costante_da_inizio_PAR(N+2,ST0,DATA,DATA_iter);
    
    if isnan(ST(end).p_OUT)% || max([ST.beta]) > 1.4
        item.errorBalje = true;
        continue
    end
    item.errorBalje = false;
    
    L_punto = sum([ST.deltaH_IS])*DATA.m_P;
    L_punto_entrante = length(ST)*ST(1).l*DATA.m_P;
    
    item.beta_TOT = ST(end).p_OUT/DATA.p0;
    item.eta_TOT = L_punto/L_punto_entrante;
    item.eta_TOT_in_out = DATA.l(1,DATA.T0,ST(end).p_OUT/DATA.p0)*DATA.m_P/L_punto_entrante;
    item.N_TOT = length(ST);
    item.beta_MEDIO = item.beta_TOT^(1/length(ST));
    item.T_OUT = ST(end).T_OUT;
    item.D = ST(end).D;
    item.w = ST(end).w;
    item.rpm = ST(end).rpm;
    item.bsuD_MIN = min([ST.bsuD]);
    item.bsuD_MAX = max([ST.bsuD]);
    
    iter_res(kk) = item;%%mod(kk, Ds_array_length),ceil(kk/Ds_array_length)) = item;
    
    %disp([num2str(ii/length(ws_array)/Ds_array_length*100), '% di completamento'])
end


disp('Terminato')
toc

%% plot figure

% for iter_res_item = iter_res
%     if iter_res_item.errorU || iter_res_item.errorBalje
%         if isempty(ws_Ds_array_error) || ~ismember([ws, Ds], ws_Ds_array_error, 'rows')
%             ws_Ds_array_error(ws_Ds_error_index,:) = [iter_res_item.ws, iter_res_item.Ds];
%             ws_Ds_error_index = ws_Ds_error_index+1;
%         end
%     end
% end

%save('data\ws_ds_wrong.mat','ws_Ds_array_error')

iter_res_beta_medio = reshape([iter_res.beta_MEDIO],[ws_array_length,Ds_array_length]);
iter_res_beta_medio(iter_res_beta_medio==0) = 1;
imagesc(ws_array,Ds_array,iter_res_beta_medio)
set(gca, 'YDir', 'normal');

hold on
axis equal
loglog(balje_data_80_X,balje_data_80_Y)%,'*');
loglog(balje_data_75_X,balje_data_75_Y)%,'*');
loglog(balje_data_70_X,balje_data_70_Y)%,'*');
loglog(balje_data_60_X,balje_data_60_Y)%,'*');
loglog(balje_data_50_X,balje_data_50_Y)%,'*');


% iter_res_noempty = remove_empty_iter_res(iter_res);
% 
% 
% Z_iter_res = zeros(length(ws_array),length(Ds_array));
% 
% for iter_res_item = iter_res
%     
%     %Z_iter_res(find(ws_array==iter_res_item.ws),find(Ds_array==iter_res_item.Ds)) = iter_res_item.beta_TOT;
%     if ~isnan( iter_res_item.beta_TOT)
%         Z_iter_res(find(ws_array==iter_res_item.ws),find(Ds_array==iter_res_item.Ds)) = iter_res_item.beta_TOT;
%     else Z_iter_res(find(ws_array==iter_res_item.ws),find(Ds_array==iter_res_item.Ds)) = 1;
%     end
%     
% end
% 
% %imagesc(([iter_res.ws]),([iter_res.Ds]),Z_iter_res);
% 
% 
% imagesc(ws_array,Ds_array,Z_iter_res);

%% save data
% name = ['data\data',num2str(ws_array(1)),'-',num2str(ws_array(end)),'_',num2str(Ds_array(1)),'-',num2str(Ds_array(end)),' ',num2str(length(ws_array)),'X',num2str(length(Ds_array))];
% 
% savefig([name,'.fig']);
% 
% save([name,'.mat'],'Ds_array','ws_array','iter_res','iter_res_noempty');

%% ottimo di beta medio
mask_opt = [iter_res.beta_MEDIO]==max([iter_res.beta_MEDIO]);
item_opt = iter_res(mask_opt);
item_opt = item_opt(1);
item_beta_medio_OPT  = iter_res(mask_opt);

disp('----------------------------------------------------')
disp(['Ottimo di beta medio:'])
disp(['ws = ', num2str(item_opt.ws)])
disp(['Ds = ', num2str(item_opt.Ds)])
disp(['w = ', num2str(item_opt.w)])
disp(['D = ', num2str(item_opt.D)])
disp(['rpm = ', num2str(item_opt.rpm)])
disp(['Beta medio = ', num2str(item_opt.beta_MEDIO)])
disp(['N con beta medio = ', num2str(log(beta_TOT)/log(item_opt.beta_MEDIO))])
disp(['Eta totale = ', num2str(item_opt.eta_TOT)])
disp(['Eta totale in out = ', num2str(item_opt.eta_TOT_in_out)])
disp(['T out = ', num2str(item_opt.T_OUT)])

%% ottimo di rendimento
mask_opt = [iter_res.eta_TOT]==max([iter_res.eta_TOT]);
item_opt = iter_res(mask_opt);
item_opt = item_opt(1);
item_eta_tot_OPT  = iter_res(mask_opt);

disp('----------------------------------------------------')
disp(['Ottimo di rendimento:'])
disp(['ws = ', num2str(item_opt.ws)])
disp(['Ds = ', num2str(item_opt.Ds)])
disp(['w = ', num2str(item_opt.w)])
disp(['D = ', num2str(item_opt.D)])
disp(['rpm = ', num2str(item_opt.rpm)])
disp(['Beta medio = ', num2str(item_opt.beta_MEDIO)])
disp(['N con beta medio = ', num2str(log(beta_TOT)/log(item_opt.beta_MEDIO))])
disp(['Eta totale = ', num2str(item_opt.eta_TOT)])
disp(['Eta totale in out = ', num2str(item_opt.eta_TOT_in_out)])
disp(['T out = ', num2str(item_opt.T_OUT)])

%% ottimo di rendimento in out
mask_opt = [iter_res.eta_TOT_in_out]==max([iter_res.eta_TOT_in_out]);
item_opt = iter_res(mask_opt);
item_opt = item_opt(1);
item_eta_tot_in_out_OPT  = iter_res(mask_opt);

disp('----------------------------------------------------')
disp(['Ottimo di rendimento in out:'])
disp(['ws = ', num2str(item_opt.ws)])
disp(['Ds = ', num2str(item_opt.Ds)])
disp(['w = ', num2str(item_opt.w)])
disp(['D = ', num2str(item_opt.D)])
disp(['rpm = ', num2str(item_opt.rpm)])
disp(['Beta medio = ', num2str(item_opt.beta_MEDIO)])
disp(['N con beta medio = ', num2str(log(beta_TOT)/log(item_opt.beta_MEDIO))])
disp(['Eta totale = ', num2str(item_opt.eta_TOT)])
disp(['Eta totale in out = ', num2str(item_opt.eta_TOT_in_out)])
disp(['T out = ', num2str(item_opt.T_OUT)])



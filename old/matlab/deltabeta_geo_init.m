
[m_data_X, m_data_Y, m_POLY] = load_csv_poly('data\m.csv');
[n_data_X, n_data_Y, n_POLY] = load_csv_poly('data\n_10.csv');
[b_data_X, b_data_Y, b_POLY] = load_csv_poly('data\b.csv');

[k_i_t_data_X, k_i_t_data_Y, k_i_t_POLY] = load_csv_poly('data\k_i_t.csv');
[k_delta_t_data_X, k_delta_t_data_Y, k_delta_t_POLY] = load_csv_poly('data\k_delta_t.csv');
[i0_10_data_X, i0_10_data_Y, i0_10_POLY] = load_csv_poly('data\i0_10_10.csv');
[delta0_10_data_X, delta0_10_data_Y, delta0_10_POLY] = load_csv_poly('data\delta0_10_10.csv');

ii = 0;
%T_array_temp = [];
for T = T_array
    
    ii = ii+1;
    
    T.GEO.m = polyval(m_POLY, abs(T.beta1));
    T.GEO.n = polyval(n_POLY, abs(T.beta1));
    
    T.GEO.sigma = 1;%%%%
    T.GEO.b = polyval(b_POLY, abs(T.beta1));
    T.GEO.t_su_c = 0.1;
    
    T.GEO.k_delta_sh = 1;
    T.GEO.k_delta_t = polyval(k_delta_t_POLY, T.GEO.t_su_c);
    T.GEO.delta0_10 = polyval(delta0_10_POLY, abs(T.beta1));
    T.GEO.delta0 = T.GEO.k_delta_sh*T.GEO.k_delta_t*T.GEO.delta0_10;
    
    T.GEO.k_i_sh = 1;
    T.GEO.k_i_t = polyval(k_i_t_POLY, T.GEO.t_su_c);
    T.GEO.i0_10 =  polyval(i0_10_POLY, abs(T.beta1));
    T.GEO.i0 = T.GEO.k_i_sh*T.GEO.k_i_t*T.GEO.i0_10;
    
    T.deltaBeta_geo = (T.deltaBeta - T.GEO.i0-T.GEO.delta0)/(1-T.GEO.m/T.GEO.sigma^T.GEO.b);
    T.deltaBeta_geo_deg = rad2deg(T.deltaBeta_geo);
    
    T.delta = T.GEO.delta0 + T.GEO.m *T.deltaBeta_geo/T.GEO.sigma^T.GEO.b;
    T.delta_deg = rad2deg(T.delta);
    T.i = T.GEO.i0 + T.GEO.n*T.deltaBeta_geo;
    T.i_deg = rad2deg(T.i);
    T.profile = create_profile(T.deltaBeta_geo, T.GEO.t_su_c);
    
    T_array_temp(ii) = T;
end

T_array = T_array_temp;
%T_array = [T_hub, T_mid, T_tip];
T_hub = T_array(1);
T_mid = T_array(2);
T_tip = T_array(3);


%DF = (1-abs(T.w2/T.w1)) + (T.v2T-T.v1T)/2/T.w1/sigma

%DF = (1-cos(T.beta1)/cos(T.beta2)) + abs((tan(-T.beta1)-tan(T.beta2))*cos(T.beta1)/2/sigma)

clear T_array_temp ii
clear m_data_X m_data_Y %m_POLY
clear n_data_X n_data_Y %n_POLY
clear b_data_X b_data_Y %b_POLY
clear k_i_t_data_X k_i_t_data_Y %k_i_t_POLY
clear k_delta_t_data_X k_delta_t_data_Y %k_delta_t_POLY
clear i0_10_data_X i0_10_data_Y %i0_10_POLY
clear delta0_10_data_X delta0_10_data_Y %delta0_10_POLY
clear delta delta_deg i i_deg
clear m n b t_su_c
clear k_delta_sh k_delta_t = 1 delta0_10 delta0
clear k_i_sh k_i_t i0_10 i0




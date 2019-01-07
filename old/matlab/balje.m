function [ eta ] = balje( ws,ds )
    
    global my_warning
    global balje_data_80_X balje_data_80_Y
    global balje_data_75_X balje_data_75_Y
    global balje_data_70_X balje_data_70_Y
    global balje_data_60_X balje_data_60_Y
    global balje_data_50_X balje_data_50_Y
    %global balje_data_axial_X balje_data_axial_Y
    
    if inpolygon(ws,ds,balje_data_80_X,balje_data_80_Y);
        eta = 0.8;
    elseif inpolygon(ws,ds,balje_data_75_X,balje_data_75_Y);
        eta = 0.75;
    elseif inpolygon(ws,ds,balje_data_70_X,balje_data_70_Y);
        eta = 0.7;
    elseif inpolygon(ws,ds,balje_data_60_X,balje_data_60_Y);
        eta = 0.6;
    elseif inpolygon(ws,ds,balje_data_50_X,balje_data_50_Y);
        eta = 0.5;
    else eta = nan;
        if my_warning
            disp(['attenzione ws = ', num2str(ws),' e ds = ',num2str(ds),' causano problemi']);
        end
    end
    
end


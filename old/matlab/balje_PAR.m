function [ eta ] = balje_PAR( ws,ds,DATA )

    if inpolygon(ws,ds,DATA.balje_data_80_X,DATA.balje_data_80_Y);
        eta = 0.8;
    elseif inpolygon(ws,ds,DATA.balje_data_75_X,DATA.balje_data_75_Y);
        eta = 0.75;
    elseif inpolygon(ws,ds,DATA.balje_data_70_X,DATA.balje_data_70_Y);
        eta = 0.7;
    elseif inpolygon(ws,ds,DATA.balje_data_60_X,DATA.balje_data_60_Y);
        eta = 0.6;
    elseif inpolygon(ws,ds,DATA.balje_data_50_X,DATA.balje_data_50_Y);
        eta = 0.5;
    else eta = nan;
    end
    
end

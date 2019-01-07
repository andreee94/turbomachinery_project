function [ X,Y,POLY ] = load_csv_poly( name )

    data = csvread(name);
    X = data(:,1);
    Y = data(:,2);
    POLY = polyfit(X,Y,6);

end


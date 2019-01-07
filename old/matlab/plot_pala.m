%richiede evoluzione triangoli

l = length(deltabeta_evoluzione);

figure(ff); ff = ff+1;



for ii = 1:l
    P = create_profile(deltabeta_evoluzione(ii), 0.1);
    [x,y,z] = plot_profile_single(P, beta1_evoluzione(ii), d_array(ii),true);
    if ii == 1
        X = zeros(l,length(x));
        Y = zeros(l,length(y));
        Z = zeros(l,length(z));
    end
    X(ii,: )=x;
    Y(ii,: )=y;
    Z(ii,: )=z;
end

surf(X*0.2,Y*0.2,Z);%,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
%camlight left
colormap winter
axis equal



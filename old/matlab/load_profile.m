function [ P ] = load_profile(name )
    
    
    %     fid = fopen(name);
    %
    %     tline = fgetl(fid);
    %     while ischar(tline)
    %         %disp(tline);
    %         tline = fgetl(fid);
    %     end
    %
    %     fclose(fid);
    
    % potrebbe essere necessario aggiungere .dat
    [x, y] = textread(name,'%n %n','headerlines',1);
    
    l = length(y);
    if mod(l,2) == 1
        midline = 0.5*(y(1:(l+1)/2) + flip(y((l+1)/2:end)));
    else 
        midline = 0.5*(y(1:l/2) + flip(y(l/2+1:end)));
    end
    
    if all(midline==0)
        m = inf;
        p = 0;
    else
        m = max(midline);
        p = x(midline==m);
    end
    
    % x, y, m, p
    %using struct
    P.name = name;
    P.x = x;
    P.y = y;
    P.midline = midline;
    P.m = m;
    P.p = p;
    [P.camber, P.camber_deg] = camber_da_profilo(m,p);
    
end


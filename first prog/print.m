function [] = print(var)
    % implementare allineamento dell'uguale a numero
    % fissato di caratteri
    if isstring(var)
        disp(string);
    else
        disp([inputname(1), ' = ', num2str(round(var, 10))]);
    end
end


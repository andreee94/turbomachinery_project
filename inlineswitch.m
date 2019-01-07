function c = inlineswitch(condition, params , values)
    
    c = values(1);
    
    for p = params
        if condition == p
            c = values(params == p);
        end
    end
    
end


function [ array ] = remove_empty_iter_res( array_in)
    
    empty_elems = arrayfun(@(s) ~isempty(s.beta_TOT),array_in);
    
    array=array_in(empty_elems);
    
end


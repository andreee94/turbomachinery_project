function [] = myprint(test, varargin)
    % implementare allineamento dell'uguale a numero
    % fissato di caratteri
    
    
    if ~isfield(test, 'print')
       print =true;
       varargin = test;
    else 
        print = test.print;
    end
    
    if print
        for k = 1:nargin-1
            var = varargin{k};
            if ischar(var)
                disp(var);
            else
                %name = evalin('caller',['inputname(' num2str(k+1) ');']);
                name = inputname(k+1);
                disp([name, ' = ', num2str(round(var, 10))]);
            end
            
        end
    end
end


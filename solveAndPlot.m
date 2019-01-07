function [etas] = solveAndPlot(gl, test, array1, array2, prop1, prop2, propOut, plotmin, plotmax, normalize, reverse)
    
    if nargin < 10
        normalize = false;
    end
    
    if nargin < 11
        reverse = 'no';
    end
    
    if ~isempty(array2)
        etas = zeros(length(array1), length(array2));
        legendString = cell(length(array1), 1);
        for ii = 1:length(array1)
            gl = setfield(gl, prop1, array1(ii));
            tic
            parfor kk = 1:length(array2)
                %kk
                tempgl = gl;
                tempgl = setfield(tempgl, prop2, array2(kk));
                
                res = solve_turbine(tempgl, test);
                etas(ii, kk) = getfield(res, propOut);
            end
            toc
            disp(['Progress of outer loop = ', num2str(ii / length(array1) * 100), '%'])
            if normalize
                etas(ii, :) = etas(ii, :) / etas(ii, 1);
            end
            plotMinMax(array2(etas(ii, :) ~= 0), etas(ii, etas(ii, :) ~= 0), plotmin, plotmax);
            legendString{ii} = [prop1, ' = ', num2str(array1(ii))];
        end
        
        xlabel(prop2)
        ylabel(propOut)
        
        if strcmp(propOut, 'etaTT')
            yticks(0.7:0.01:1);
            if min(min(etas)) > 0.9 && max(max(etas) < 1)
                ylim([0.9 1])
            end
        end
        
        f=get(gca,'Children');
        step = 1+plotmin+plotmax;
        f = f(end:-step:1);
        %f(1+plotmin+plotmax:1+plotmin+plotmax:end);
        legend(f,legendString)
    end
    
    if strcmp(reverse, 'reverse')
        figure
        solveAndPlot(gl, test, array2, array1, prop2, prop1, propOut, plotmin, plotmax, normalize, 'no')
    end
    
end
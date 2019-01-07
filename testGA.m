[ p_opt, fit_opt, res ] = MAIN_GA();

%%
%fits = zeros(size(res()));

for ii = 1:50
    plot(1:20, res(ii).fitness)
    hold on
    %fits(ii) = res(ii).fits 
end
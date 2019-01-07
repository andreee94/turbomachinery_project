% for now constraint only between min and max
% myga find the MAXIMUM
function [ p_opt, fit_opt, res ] = myga(f, constraints, var_type, options)
    
    % options:
    % options.population_number
    % options.selection_quantity % from 0 to 1
    % options.mutation_min
    % options.mutation_max
    % options.mutation_threshold
    % options.max_iterations
    % options.crossover_first_threshold
    
    pop = init_population(options.population_number, constraints, var_type);
    fitness = fitness_population(options.population_number, pop, f);
    [fitness, pop] = sort_population(pop, fitness);
    disp('random population generated')
    kk = 1;
    res(kk).pop = pop;
    res(kk).fitness = fitness;
    
    convergence = false;
    while ~convergence
        disp('apply selection')
        [fitness, pop] = apply_selection(options.population_number, pop, fitness, options);
        
        disp('complete population')
        [fitness, pop] = complete_population(f, options.population_number, pop, fitness, options);
        
        disp('sort population')
        [fitness, pop] = sort_population(pop, fitness);
    
        disp(['iteration ' , num2str(kk), 'done.'])
        kk = kk + 1;
        res(kk).pop = pop;
        res(kk).fitness = fitness;
        if kk >= options.max_iterations
           convergence = true; 
        end
    end
    
    p_opt = pop(1, :);
    fit_opt = fitness(1);
    
end

function [fit, pop] = complete_population(f, n, pop, fit, options)
    %myga_const.CROSSOVER_STATEGY_ALWAYS_BEST = 11;
    first_nan_index = find(isinf(fit) == 1, 1); % only the first nan elements
    
    if options.crossover_strategy == myga_const.CROSSOVER_STATEGY_ALWAYS_BEST
        for ii = first_nan_index:n
            index = mod(ii - first_nan_index + 1, n - first_nan_index) + 1;
            pop(ii, :) = crossover(pop(1, :), pop(ii, :), options.crossover_first_threshold);
            fit(ii) = f(pop(ii, :));
        end
    end
end

function p = crossover(p1, p2, threshold)
    rands = rand(size(p1));
    
    mask = rands < threshold;
    %antimask = 1 - mask;
    
    p = p2;
    p(mask) = p1(mask);
    %p(antimask) = p2(antimask);
end

function [fit, pop] = apply_selection(n, pop, fit, options)
    rands = rand(n, 1);
    comparison = linspace(1, options.selection_quantity / 2, n)';
    
    fit(rands >= comparison) = -inf; % part of the population
    [fit, pop] = sort_population(pop, fit);
end

function pop = apply_mutation(n, pop, options)
    % TODO now is linear but can be more complex or also user defined
    mutation_probability = linspace(options.mutation_min, options.mutation_max, n); 
    rands = rand(n, 1);
end

function fit = fitness_population(n, pop, f)
    fit = zeros(n, 1);
    for ii = 1:n
        fit(ii) = f(pop(ii, :));
    end
end

function [fit, pop] = sort_population(pop, fit)
   [fit, index] = sort(fit, 1, 'descend');
   pop = pop(index, :);
end

function pop = init_population(n, constraints, var_type)
    %VAR_TYPE_INT = 11;
    %VAR_TYPE_DOUBLE = 9;
    % var_type = 11 INT
    % var_type = 9 DOUBLE
    var_num = size(constraints, 1);
    pop = zeros(n, var_num);
    
    rands = rand(size(pop));
    
    for ii = 1:var_num
       pop(:, ii) = constraints(ii, 1) + rands(:, 1) .* (constraints(ii, 2) - constraints(ii, 1));
       if var_type(ii) == myga_const.VAR_TYPE_INT
            pop(:, ii) = round(pop(:, ii)); 
       end
    end
end




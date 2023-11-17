function is_invaded = sim__invasibility_map(params, outfile, p1s, ...
    log10deltas, rho_pre)
% Simulates invasibility tests in a 2-species co-culture for the given
% nutrient profiles and adaptor's sensing tolerances values in 'p1s' and
% 'log10deltas', respectively.
% At each parameter combination, a batch starts with a relative population
% 'rho_pre' (~< 1) of a 'pre-existing' species, while introducing a
% miniscule relative population 1 - rho_pre of the other species into the
% culture. Species switch roles, and invasibility is characterized using
% batch-end populations at both cases. The pre-existing species is invaded
% if it does not grow in relative population throughout the batch.
% Returns a matrix 'is_invaded' of which the ij-th element corresponds to
% invasibility character for log10deltas(i), p1s(j). The values are:
% (1) Mutual invasibility.
% (2) Non-adaptor invaded only.
% (3) Adaptor invaded only.
% It will optionally plot the resulting invasibility map (if params.plt =
% 1), and will save the data in a file in path 'outfile' (if not empty).

%%% NOTE: Here, a co-culture of a non-adaptor (Sp. 1) and an adaptor
%%% (Sp. 2) is considered. However, this can be easily generalized to any
%%% 2-species combination.

%% Check for existing data
if ~ isempty(outfile)
    if exist(outfile, 'file') == 2
        disp('Invasibility data exists, returning');
        is_invaded = [];
        return
    end
end

%% Warning 
if sum(is_adaptor == 1) && p ~= 2
    warning('For adaptation, supply 2 nutrients (params.p = 2)!')
    is_invaded = [];
    return
end

%% Simulate invasibility tests
plt = params.plt; % Keep for after tests
params.plt = 0;   % Avoid plotting individual batches

is_invaded = nan(length(log10deltas), length(p1s));
for i = 1:length(log10deltas)
    for j = 1:length(p1s)
        % Plug in moving-parameter values
        params.log10delta(2) = log10deltas(i);
        params.P = [p1s(j); 1 - p1s(j)];
        
        % 1st species pre-existing
        params.b0 = [rho_pre; 1 - rho_pre];
        rho = sim__invasibility(params);
        if 1 - rho_pre <= rho(2)
            is_invaded1 = 1;
        else
            is_invaded1 = 0;
        end
        
        % 2st species pre-existing
        params.b0 = [1 - rho_pre; rho_pre];
        rho = sim__invasibility(params);
        if 1 - rho_pre <= rho(1)
            is_invaded2 = 1;
        else
            is_invaded2 = 0;
        end
        
        % Determine invasibility
        if is_invaded1 + is_invaded2 == 2
            is_invaded(i, j) = 1;
        elseif is_invaded1 == 1
            is_invaded(i, j) = 2;
        else
            is_invaded(i, j) = 3;
        end
    end
end

%% Plot
if plt
    plot__invasibility_map(p1s, log10deltas, is_invaded, params, rho_pre)
end

%% Save results
if ~ isempty(outfile)
    save(outfile, 'is_invaded', 'params', 'p1s', 'log10deltas', 'rho_pre')
    disp(['Data has been saved to ', outfile])
end

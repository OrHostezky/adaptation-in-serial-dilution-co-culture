function output = sim__serial__interbatch(params, outfile)
% Simulates a serial-dilution course for the given parameter values in
% 'params', using MATLAB's built-in ODE solver. With an appropriate value
% for 'max_batches', the process will continue until reaching a steady
% state of inter-batch dynamics, or other distinct long-time behavior. It
% will also stop at any point if only 1 species survives.
% Returns the inter-batch data (values at dilutions) as a structure
% 'output'. It will optionally: (1) plot the inter-batch dynamics and the
% following intra-batch dynamics (if params.plt = 1), (2) compute growth
% integrals (if params.ss = 1), (3) plot total enzyme concetrations (if
% params.pltE = 1), (4) save the data in a file in path 'outfile' (if not
% empty) every 10^4 batches and at a steady state (or other distinct long-
% time behavior).

%% Set simulation parameters
% Nutrients
c0 = 10^params.log10c0;
K = params.K;
p = params.p;
P = params.P;

% Species
rho0 = params.rho0;
m = params.m;
b0 = params.b0;

% Enzymatics
E = params.E;

if ~ isfield(params, 'alpha0')
    if isrow(E)
        alpha0 = zeros(m, p) + E' / p;
    else
        alpha0 = zeros(m, p) + E / p;
    end
else
    alpha0 = params.alpha0;
end

if isfield(params, 'log10delta')
    delta = 10.^params.log10delta;
else
    delta = c0 * ones(1, m);
end

if isfield(params, 'ctrl0')
    ctrl0 = params.ctrl0;                % Given
else
    ctrl0 = (P(1) >= P(2)) * ones(m, 1); % Abundant nutrient
end

is_adaptor = params.is_adaptor;
model = params.model;

% Control
if isfield(params, 'max_batches')
    max_batches = params.max_batches;
else
    max_batches = round(1e6 * (c0 / rho0)^(- 1 / 6));
end

%% Warning 
if sum(is_adaptor == 1) && p ~= 2
    warning('For adaptation, supply 2 nutrients (params.p = 2)!')
    output = [];
    return
end

%% Initialize
% Initialize storage vectors
tb = 1; % Batch no. variable
rho = nan(m, max_batches + 1);
rho(:, tb) = b0;
alpha = nan(m, p, max_batches + 1);
alpha(:, :, tb) = alpha0;
ctrl = nan(m, max_batches + 1);
ctrl(:, tb) = ctrl0;
switch_store = nan(m, max_batches);

% Plug in existing data (if such exists)
if ~ isempty(outfile)
    if exist(outfile, 'file') == 2
        disp('Attempts to read existing outfile')
        try
            load(outfile,'output')
            if isfield(output, 'ss')
                disp('Steady-state data already exists, returning')
                output = [];
                return
            end
            
            % Collect data
            tb = size(output.rho, 2);
            rho(:, 1:tb) = output.rho;
            alpha(:, :, 1:tb) = output.alpha;
            ctrl(:, 1:tb) = output.ctrl;
            switch_store(:, 1 : tb - 1) = output.switch_store;
            if isfield(output, 'pop_frac')
                pop_frac = output.pop_frac;
            end
            
            disp(['Continues existing simulation from batch ', int2str(tb)])
        
        catch
            warning('Failed loading, returning')
            output = [];
            return
        end
    end
end

%% Simulate serial-dilution
ss1 = 0; ss2 = 0;                          % Steady-state indicators
cond = round(5e2 * (c0 / rho0)^(- 1 / 3)); % Steady-state condition
xprev = [rho(:, tb); ctrl(:, tb)];
pop_frac_mean = rho(:, tb);

while tb <= max_batches && ss1 < cond && ss2 < cond && sum(xprev(1:m) > 0) > 1 
    % Choose initial control states for coming batch
    if model == 1
        ctrl_b = xprev(end - m + 1 : end);
    elseif model == 2 || model == 3
        ctrl_b = ctrl0;
    end

    % Simulate batch
    [xfinal, ~, ~, switch_num] = sim__batch(c0, K, P, rho0, rho(:, tb), E, ...
        delta, alpha(:, :, tb), ctrl_b, is_adaptor, model, 0, tb);
    
    % Discard nutrient concentrations, neglect extinct populations
    xfinal(1:p) = [];
    xfinal(1:m) = xfinal(1:m) / sum(xfinal(1:m));
    xfinal(xfinal(1:m) < 1e-20) = 0;
    xfinal((1 - xfinal(1:m)) < 1e-20) = 1;
    
    % Plug final values into dynamic variables
    rho(:, tb + 1) = xfinal(1:m);
    for i = 1:m
        alpha(i, :, tb + 1) = xfinal(p * (i - 1) + m + 1 : p * i + m)';
    end
    ctrl(:, tb + 1) = xfinal(end - m + 1 : end);
    switch_store(:, tb) = switch_num;
    
    % Detect "exact" populations steady state
    if max(abs(xfinal(1:m) - xprev(1:m))) < 1e-10
        ss1 = ss1 + 1;
    else
        ss1 = 0;
    end
    
    % Detect convergence of running logarithmic averages (weighing down extinctions)
    pop_frac_mean_temp = pop_frac_mean;
    pop_frac_mean = (tb * pop_frac_mean - log10(xfinal(1:m) + 1e-20)) / (tb + 1);
    if max(abs(pop_frac_mean - pop_frac_mean_temp)) < 1e-7
        ss2 = ss2 + 1; 
    else
        ss2 = 0;
    end
    
    % If a 2nd-type steady state has been detected, calculate the mean 
    % populations over the shortest of the following:
    % One inter-batch cycle / 1e4 [*] / max_batches - tb [*] ([*] - batches)
    if ss2 == cond
        x0 = xfinal(1:m);
        x = nan(2 * m, 1);
        pop_frac = xfinal(1:m);
        tbb = 1; % Sampling-batch no.
        tbb_max = min(1e4, max_batches - tb);

        while tbb <= tbb_max && ~ isequal(x, x0)
            % Choose initial control states for the coming batch
            if model == 1
                ctrl_b = xfinal(end - m + 1 : end);
            elseif model == 2 || model == 3
                ctrl_b = ctrl0;
            end
            
            % Simulate batch
            [xfinal, ~, ~, switch_num] = sim__batch(c0, K, P, rho0, ...
                rho(:, tb + tbb), E, delta, alpha(:, :, tb + tbb), ctrl_b, ...
                is_adaptor, model, 0, tb + tbb);
            
            % Discard nutrient concentrations, neglect extincted populations
            xfinal(1:p) = [];
            xfinal(1:m) = xfinal(1:m) / sum(xfinal(1:m));
            xfinal(xfinal(1:m) < 1e-20) = 0;
            xfinal((1 - xfinal(1:m)) < 1e-20) = 1;
            
            % Plug final values into dynamic variables
            rho(:, tb + tbb + 1) = xfinal(1:m);
            for i = 1:m
                alpha(i, :, tb + tbb + 1) = xfinal(p * (i - 1) + m + 1 : p * i + m)';
            end
            ctrl(:, tb + tbb + 1) = xfinal(end - m + 1 : end);
            switch_store(:, tb + tbb) = switch_num;
            
            % Running averages
            pop_frac = (tbb * pop_frac + xfinal(1:m)) / (tbb + 1);
            x = xfinal(1:m);
            
            tbb = tbb + 1;
        end

        % Check variability in population fractions
        rho_end = rho(:, tb + (1:tbb));                         % Long time populations
        rho_end_im = max(rho_end, [], 2) - min(rho_end, [], 2); % Fluctuation ranges
        
        tb = tb + tbb - 1;
    end
    
    % Update loop variables
    tb = tb + 1;                              % Batch no.
    xprev = xfinal([1:m, end - m + 1 : end]); % [rho; ctrl]
    
    % Save every 10^4 * n batches (n = 1, 2, 3, ...)
    if rem(tb - 1, 1e4) == 0
        % Output
        output.rho = rho(:, 1:tb);
        output.alpha = alpha(:, :, 1:tb);
        output.ctrl = ctrl(:, 1:tb);
        output.switch_store = switch_store(:, 1 : tb - 1);
                
        if ~ exist('pop_frac', 'var') || sum(xprev(1:m) == 0) == m - 1
            output.pop_frac = xprev(1:m);
        else
            output.pop_frac = pop_frac;
        end
        
        % Save data
        if ~ isempty(outfile)
            save(outfile, 'output', 'params')
            disp(['Partial data (until batch #', int2str(tb - 1), ...
                ') has been saved to ', outfile])
        end
    end
end

%% Output
output.rho = rho(:, 1:tb);
output.alpha = alpha(:, :, 1:tb);
output.ctrl = ctrl(:, 1:tb);
output.switch_store = switch_store(:, 1 : tb - 1);

% Relative populations
if ~ exist('pop_frac', 'var') || sum(output.rho(:, end) == 0) == m - 1
    output.pop_frac = output.rho(:, end);
else
    output.pop_frac = pop_frac;
end

%% Plot inter-batch and following intra-batch dynamics
if params.plt   
    if tb > 1
        plot__interbatch(output.rho, output.alpha, c0, P, E, delta, ctrl0, ...
            model, tb - 1)
    end

    if model == 1
        ctrl_b = ctrl(:, tb);
    elseif model == 2 || model == 3
        ctrl_b = ctrl0;
    end
    [~, t, x, ~] = sim__batch(c0, K, P, rho0, rho(:, tb), E, delta, ...
        alpha(:, :, tb), ctrl_b, is_adaptor, model, 1, tb); 
end

%% Additional checks
% Compute growth integrals
if params.ss
    [~, output.growth_int] = find_growth(t, x, K, p, m);
end

% Derive enzyme budgets from dynamic concentrations ('sanity check')
if params.pltE
    figure
    for i = 1:m
        plot(t, sum(x(p * i + m + 1 : p * (i + 1) + m, :)), '-', ...
            'DisplayName', ['$E_', int2str(i), '$'], 'LineWidth', 1)
        hold on
    end
    
    xlim([0, t(end)])
    legend('Interpreter', 'latex')
    title('Total adaptors strategy budgets', 'Interpreter', 'latex')
end

%% Determine long-time behavior ('steady-state')
if ss1 == cond || sum(output.pop_frac == 0) == m - 1
    ss = 1;     % "Exact"
elseif ss2 == cond    
    if max(rho_end_im) < 0.02
        ss = 2; % Small fluctuations
    elseif max(rho_end_im) < 0.05
        ss = 3; % Moderate fluctuations
    else
        ss = 4; % Large fluctuations
    end

    for i = 1:m
        uniq = unique(rho_end(i, :));
        % Unwanted convergence of running averages (monotonic dynamics)
        if isequal(rho_end(i, :), uniq) || isequal(rho_end(i, :), flip(uniq))
            ss = - 1;
        end 
    end
else
    ss = 0; % No steady state / unclear (simulate more batches)
end
output.ss = ss;

%% Save data
if ~ isempty(outfile)
    save(outfile, 'output', 'params')
    disp(['Complete data has been saved to ', outfile])
end

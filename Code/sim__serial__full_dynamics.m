function output = sim__serial__full_dynamics(params, outfile)
% Simulates a serial-dilution course for the given parameter values in
% 'params', using MATLAB's built-in ODE solver. The process will continue
% for the given number of batches / until only 1 species survives.
% Returns the full time-course data as a structure 'output'. It will
% optionally plot the full dynamics of the system (if params.plt = 1), and
% save the data in a file in path 'outfile' (if not empty) every 10^4
% batches and at the end of the simulation.

%%% NOTE: Uses v7.3 for saving data.

%% Check for existing data
if ~ isempty(outfile)
    if exist(outfile, 'file') == 2
        disp('Data exists, returning');
        output = [];
        return
    end
end

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

%% Simulate serial-dilution
% Initialize
tb = 1;         % Batch no. variable
t_tot = 0;      % Total time (accumulating)
t = []; x = []; % Storage vectors

rho = b0;
alpha = alpha0;
xb = [rho; ctrl0]; % Solution at dilution

while tb <= max_batches && sum(xb(1:m) > 0) > 1 
    % Choose initial control states for coming batch
    if model == 1
        ctrl = xb(end - m + 1 : end);
    elseif model == 2 || model == 3
        ctrl = ctrl0;
    end

    % Simulate batch
    [xb, tadd, xadd, ~] = sim__batch(c0, K, P, rho0, rho, E, delta, alpha, ...
        ctrl, is_adaptor, model, 0, tb);

    % Discard final nutrient concentrations, neglect extinct end populations
    xb(1:p) = [];
    xb(1:m) = xb(1:m) / sum(xb(1:m));
    xb(xb(1:m) < 1e-20) = 0;
    xb((1 - xb(1:m)) < 1e-20) = 1;
    
    % Plug final values into dynamic variables
    rho = xb(1:m);
    for i = 1:m
        alpha(i, :) = xb(p * (i - 1) + m + 1 : p * i + m)';
    end
    
    % Isolate (relative) populations and strategies
    xadd = xadd(p + 1 : end - m, :);
    rho_tot = sum(xadd(1:m, :));
    for i = 1:m
        xadd(i, :) = xadd(i, :) ./ rho_tot;
    end

    % Tailor solution
    t = [t, tadd + t_tot];
    x = [x, xadd];
    t_tot = t(end);
    
    % Update loop variables
    tb = tb + 1;                       % Batch no.
    xb = xb([1:m, end - m + 1 : end]); % [rho; ctrl]
    
    % Save every 10^4 * n batches (n = 1, 2, 3, ...)
    if rem(tb - 1, 1e4) == 0
        % Output
        output.t = t;
        output.rho = x(1:m, :);
        output.alpha = nan(m, p, length(t));
        for i = 1:m
            output.alpha(i, :, :) = permute(x(p * (i - 1) + m + 1 : p * i + m, :), ...
                [3, 1, 2]);
        end
        
        % Save data
        if ~ isempty(outfile) && exist('output', 'var')
            save(outfile, 'output', 'params', '-v7.3')
            disp(['Partial data (until batch #', int2str(tb - 1), ...
                ') has been saved to ', outfile])
        end
    end
end

%% Output
output.t = t;
output.rho = x(1:m, :);
output.alpha = nan(m, p, length(t));
for i = 1:m
    output.alpha(i, :, :) = permute(x(p * (i - 1) + m + 1 : p * i + m, :), ...
        [3, 1, 2]);
end

%% Plot
if params.plt
    plot__full_dynamics(t, x, c0, P, E, delta, model, tb - 1)
end

%% Save data
if ~ isempty(outfile) && exist('output', 'var')
    save(outfile, 'output', 'params', '-v7.3')
    disp(['Complete data has been saved to ', outfile])
end

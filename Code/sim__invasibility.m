function rho = sim__invasibility(params)
% Simulates a single, high-precision batch of growth for the given
% parameter values in 'params', using MATLAB's built-in ODE solver, with
% the purpose of charactarizing invasibility.
% Returns relative species populations at the end of the batch 'rho'. It
% will optionally plot the intra-batch dynamics (if params.plt = 1).

%% Set simulations parameters
% Nutrients
c0 = 10^params.log10c0; 
p = params.p;
P = params.P;

% Species
m = params.m;

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

%% Warning 
if sum(is_adaptor == 1) && p ~= 2
    warning('For adaptation, supply 2 nutrients (params.p = 2)!')
    rho = [];
    return
end

%% Simulate and output
% Simulate batch
[xfinal, ~, ~, ~] = sim__batch(c0, params.K, P, params.rho0, params.b0, ...
    E, delta, alpha0, ctrl0, params.is_adaptor, params.model, params.plt, ...
    1, 1);

% Isolate populations, neglect extinct populations
rho = xfinal(p + (1:m));
rho = rho / sum(rho);
rho(rho < 1e-20) = 0;
rho(1 - rho < 1e-20) = 1;

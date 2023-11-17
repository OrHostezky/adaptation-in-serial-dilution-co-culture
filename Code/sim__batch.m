function [xfinal, t, x, switch_num] = sim__batch(c0, K, P, rho0, b0, E, ...
    delta, alpha0, ctrl0, is_adaptor, model, plt, batch, high_precision)
% Simulates a single batch of growth for the given parameter values, using
% MATLAB's built-in ODE solver.
% Returns a full solution x(t), the solution at batch end, and the number
% of enzyme production switches throughout the time-course. It will
% optionally increase solver precision (if high_precision = 1), and plot
% the intra-batch dynamics (if plt = 1).

%% Initialize
p = length(P);    % Nutrient no.
m = length(b0);   % Species no.
tspan = [0, inf]; % Time span

% Determine integration tolerances
reltol = 1e-11;
abstol = 1e-10;
if exist("high_precision", "var")
    if high_precision
        reltol = 2.5e-14;
        abstol = 1e-13;
    end
end

% Set initial conditions to state vector by order:
% [c(1), ... c(p), rho(1), ... rho(m), alpha(1, 1), ... alpha(1, p), ...
%     alpha(2, 1), ... alpha(m, p)]
x0 = [c0 * P', rho0 * b0'];
for i = 1:m
    x0 = [x0, alpha0(i, :)];
end

%% Simulate
if model == 1
    % Choose enzyme production (adaptation)
    delta_c = abs(P(1) - P(2)) ./ (1 + delta / c0);   % Sensory input
    ctrl0(is_adaptor .* delta_c > 0.5) = P(1) > P(2); % Adaptation
    
    % Integrate until 1st event
    options = odeset('NonNegative', 1, 'RelTol', reltol, 'AbsTol', abstol, ...
        'Events', @(t, x) eventfun1(t, x, K, p, E, delta, ctrl0, is_adaptor));
    [t, x, te, xe, ie] = ode89(@(t, x) odefun(t, x, K, p, E, ctrl0, is_adaptor),  ...
        tspan, x0, options);
    ctrl = repmat(ctrl0', [length(t), 1]);
    
    % Repeat until non-switch event
    while max(ie) <= m
        % Initialize next integration
        tspan(1) = te(1); x0 = xe(1, :); % Start from previous event time & solution
        ctrl0(ie) = 1 - ctrl0(ie);       % Switch production of species ie
        
        % Integrate
        options = odeset('NonNegative', 1, 'RelTol', reltol, 'AbsTol', abstol, ...
            'Events', @(t, x) eventfun1(t, x, K, p, E, delta, ctrl0, is_adaptor));
        [tadd, xadd, te, xe, ie] = ode89(@(t, x) odefun(t, x, K, p, E, ctrl0, ...
            is_adaptor), tspan, x0, options);
        
        % Tailor solution
        t = [t; tadd(2:end)];
        x = [x; xadd(2:end, :)];
        ctrl = [ctrl(1 : end - 1, :); repmat(ctrl0', [length(tadd), 1])];
    end

    % Forced halt
    if max(ie) == m + 2
        disp(['Batch #', int2str(batch), ' was terminated at t = 1000'])
    end


elseif model == 2
    % Choose enzyme production (adaptation)
    delta_c = abs(P(1) - P(2)) ./ (1 + delta / c0);    % Sensory input
    ctrl_t = ctrl0;
    ctrl_t(is_adaptor .* delta_c > 0.5) = P(1) > P(2); % Adaptation
    
    % Integrate until 1st event
    options = odeset('NonNegative', 1, 'RelTol', reltol, 'AbsTol', abstol, ...
        'Events', @(t, x) eventfun2(t, x, K, p, E, delta, ctrl0, ctrl_t, is_adaptor));
    [t, x, te, xe, ie] = ode89(@(t, x) odefun(t, x, K, p, E, ctrl_t, is_adaptor), ...
            tspan, x0, options);
    ctrl = repmat(ctrl_t', [length(t), 1]);
    
    % Repeat until non-switch event
    while max(ie) <= m
        % Initialize next integration
        tspan(1) = te(1); x0 = xe(1, :); % Start from previous event time & solution
        ctrl_t(ie) = 1 - ctrl_t(ie);     % Switch production of species ie
        
        % Integrate
        options = odeset('NonNegative', 1, 'RelTol', reltol, 'AbsTol', abstol, ...
            'Events', @(t, x) eventfun2(t, x, K, p, E, delta, ctrl0, ctrl_t, is_adaptor));
        [tadd, xadd, te, xe, ie] = ode89(@(t, x) odefun(t, x, K, p, E, ctrl_t, ...
            is_adaptor), tspan, x0, options);
        
        % Tailor solution
        t = [t; tadd(2:end)];
        x = [x; xadd(2:end, :)];
        ctrl = [ctrl(1 : end - 1, :); repmat(ctrl_t', [length(tadd), 1])];
    end

    % Forced halt
    if max(ie) == m + 2
        disp(['Batch #', int2str(batch), ' was terminated at t = 1000'])
    end

    
elseif model == 3
    % Choose enzyme production (adaptation)
    delta_c = P(p - ctrl0)' ./ (P(p - ctrl0)' + delta / c0); % Sensory input
    ctrl_t = ctrl0;
    idx = is_adaptor .* delta_c < 0.5;
    ctrl_t(idx) = 1 - ctrl0(idx);                            % Adaptation

    % Integrate until1 1st event
    options = odeset('NonNegative', 1, 'RelTol', reltol, 'AbsTol', abstol, ...
        'Events', @(t, x) eventfun3(t, x, K, p, E, delta, ctrl0, ctrl_t, is_adaptor));
    [t, x, te, xe, ie] = ode89(@(t, x) odefun(t, x, K, p, E, ctrl_t, is_adaptor),  ...
        tspan, x0, options);
    ctrl = repmat(ctrl_t', [length(t), 1]);
    
    % Repeat until non-switch event
    while max(ie) <= m
        % Initialize next integration
        tspan(1) = te(1); x0 = xe(1, :); % Start from previous event time & solution
        ctrl_t(ie) = 1 - ctrl_t(ie);     % Switch production of species ie
        
        % Integrate
        options = odeset('NonNegative', 1, 'RelTol', reltol, 'AbsTol', abstol, ...
            'Events', @(t, x) eventfun3(t, x, K, p, E, delta, ctrl0, ctrl_t, is_adaptor));
        [tadd, xadd, te, xe, ie] = ode89(@(t, x) odefun(t, x, K, p, E, ctrl_t, ...
            is_adaptor), tspan, x0, options);
        
        % Tailor solution
        t = [t; tadd(2:end)];
        x = [x; xadd(2:end, :)];
        ctrl = [ctrl(1 : end - 1, :); repmat(ctrl_t', [length(tadd), 1])]; 
    end

    % Forced halt
    if max(ie) == m + 2
        disp(['Batch #', int2str(batch), ' was terminated at t = 1000'])
    end

end

%% Output
x = [x, ctrl];             % Plug in control states to solution
t = t'; x = x';            % Propagate time through columns
xfinal = x(:, end);        % Final values
state = x(1, :) > x(2, :); % Nutrient "state"

% Count enzyme production switches
switch_num = nan(m, 1);
for i = 1:m
    switch_num(i) = sum(abs(diff(x(end - m + i, :))));
end

%% Plot
if plt
    plot__intrabatch(t, x, state, c0, K, P, E, delta, ctrl0, model, batch)
end

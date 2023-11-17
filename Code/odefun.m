function dxdt = odefun(~, x, K, p, E, ctrl, is_adaptor)
% An ODE system in the form of dx/dt = f(x, params) ('~' - unused time
% variable). Solves for nutrients 'c', species populations 'rho', and their
% metabolic strategies 'alpha'.
%%% NOTE: alpha is constant when p is not 2.

m = length(E); % Species no.

%% Extract variables from state vector
c = x(1:p);
rho = x(p + 1 : p + m);
alpha = nan(m, p);
for i = 1:m
    alpha(i, :) = x(p * i + m + 1 : p * (i + 1) + m)';
end

%% Compute derivative
dxdt = zeros(size(x));

% Populations growth
growth_rate = alpha * (c ./ (K + c));
dxdt(p + 1 : p + m) = rho .* growth_rate;

% Nutrient depletion rates
dxdt(1:p) = - (c ./ (K + c)) .* transpose(sum(alpha .* repmat(rho, [1, p])));

% Adaptors strategy changes
if p == 2
    for i = 1:m
        dxdt(p * i + m + 1 : p * (i + 1) + m) = is_adaptor(i) * growth_rate(i) ...
            * (E(i) * [ctrl(i); 1 - ctrl(i)] - alpha(i, :)');
    end
end

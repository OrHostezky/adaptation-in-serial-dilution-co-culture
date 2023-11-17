function [iszero, isterm, dir] = eventfun1(t, x, K, p, E, delta, ctrl, is_adaptor)
% Stops integration when: (1) one of the species is about to switch enzyme
% production (Model-1 adaptation), (2) the norm of the changes in the
% nutrients and populations is smaller than 10^-8, (3) t = 1000.

m = length(E); % Species no.

if p == 2
    c = x(1:p);
    delta_c = diff(c) ./ (delta + sum(c)); % Sensory input
else
    delta_c = 0;
end

dxdt = odefun(t, x, K, p, E, ctrl, is_adaptor); % Derivative

%% Output
iszero = [0.5 - sign(ctrl' - 0.5) .* is_adaptor .* delta_c, ...
    norm(dxdt(1 : p + m)) - 1e-8, 1000 - t];
isterm = ones(1, m + 2);
dir = - ones(1, m + 2); % Last can be 0, doesn't matter

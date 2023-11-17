function [growth_rate, growth_int] = find_growth(t, x, K, p, m)
% Computes species' growth rates and their time integrals (total growth
% over time), based on growth-batch solution given by x(t), and used
% parameter values.

c = x(1:p, :); % Nutrient contentrations
x(1 : p + m, :) = [];

alpha_i = nan(m, p);
growth_rate = nan(m, length(t));
for ti = 1:length(t)
    for i = 1:m
        alpha_i(i, :) = x(p * (i - 1) + 1 : p * i, ti)';
    end
    c_i = c(:, ti);
    growth_rate(:, ti) = alpha_i * (c_i ./ (K + c_i));
end
growth_int = trapz(t, growth_rate, 2);

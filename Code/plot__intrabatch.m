function plot__intrabatch(t, x, state, c0, K, P, E, delta, ctrl0, model, batch)
% Plots the intra-batch dynamics given by the solution x(t), and used
% parameter values, icluding: (1) nutrient concentrations and species
% populations, (2) produced enzymes, metabolic strategies, and nutrient
% state, (3) species' growth rates and growth integrals (total growth over
% batch). Saves the resulting figure in '../Plots/Raw/'.

%% Initialize
p = length(P); % Nutrient no.
if p == 2
    nutrient_cmap = ['k'; '#77AC30'];
else
    nutrient_cmap = colormap(pink(p));
end

m = length(E); % Species no.
if m == 2
    cmap = ['r'; 'b'];
else
    cmap = colormap(parula(m));
end

figure
sgtitle(['Intra-batch dynamics for batch #', int2str(batch)], ...
    'Interpreter', 'latex')

%% Plot nutrients and species
subplot(1, 3, 1)
labels = strings(p + m, 1);
for i = 1:p
    plot(t, x(i, :), 'Color', nutrient_cmap(i, :))
    labels(i) = ['$c_{', int2str(i), '}(t)$'];
    hold on
end
for i = 1:m
    plot(t, x(i + p, :), 'Color', cmap(i, :), 'LineWidth', 1)
    labels(i + p) = ['$\rho_{', int2str(i), '}(t)$'];
    hold on
end

xlim([0, t(end)])
ylim([0, inf])
xlabel('Time', 'Interpreter', 'latex')
ylabel('Amount', 'Interpreter', 'latex')
legend(labels, 'Interpreter', 'latex', 'Location', 'best')

%% Plot control states and strategies
ax = subplot(1, 3, 2);
labels = strings(2 * m + 1, 1);
for i = 1:m
    plot(t, x(end - m + i, :), '--', 'Color', cmap(i, :), 'LineWidth', 2)
    labels(i) = ['$Ctrl_{', num2str(i), '}(t)$'];
    hold on
end
for i = 1:m
    plot(t, x(p * i + m + 1, :) / E(i), 'Color', cmap(i, :), 'LineWidth', 1)
    labels(i + m) = ['$\alpha_{', num2str(i), ',1}(t) / E_{', int2str(i), '}$'];
    hold on
end
plot(t, state, 'Color', 'k')
labels(end) = '$c_1(t) > c_2(t)$';

xlim([0, t(end)])
ylim([0, 1])
xlabel('Time', 'Interpreter', 'latex')
ylabel('Enzyme-1 fraction', 'Interpreter', 'latex')
legend(labels, 'Interpreter', 'latex', 'Location', 'southeast')

%% Plot growth rates
% Compute growth rates
[growth_rate, growth_int] = find_growth(t, x, K, p, m);

% Display growth integrals
str = strings(m + 1, 1);
str(1) = 'Growth integrals:';
for i = 1:m
    str(i + 1) = ['Sp. ', int2str(i), ': $', num2str(growth_int(i)), '$'];
end
annotation('texbatchox', [.72, 0, .3, .25], 'String', str, 'Interpreter', ...
    'latex', 'FibatchoxToText', 'on', 'FontSize', 8)

% Plot
subplot(1, 3, 3)
labels = strings(m, 1);
for i = 1:m
    plot(t, growth_rate(i, :), 'Color', cmap(i, :), 'LineWidth', 1)
    labels(i) = ['$g_{', int2str(i), '}(t)$'];
    hold on
end

%set(gca, 'YScale', 'log'); % Use when rates decay exponentially
xlim([0, t(end)])
xlabel('Time', 'Interpreter', 'latex')
ylabel('Growth rate', 'Interpreter', 'latex')
legend(labels, 'Interpreter', 'latex')

%% Subtitle
sbtitle_prefix = ['$\textrm{Model} \;', int2str(model), '\quad c_0 = \; 10^{', ...
    num2str(log10(c0)), '} \quad c_1(0)/c_0 = \;', num2str(P(1)), ...
    '\quad E_\sigma = \;', mat2str(E), '\quad'];
sbtitle_suffix = ['\Delta_\sigma = \; 10^{', mat2str(log10(delta)), '}$'];

if model == 1
    sbtitle = [sbtitle_prefix, sbtitle_suffix];
elseif model == 2 || model == 3
    sbtitle = [sbtitle_prefix, 'c_{preferred}^{\sigma} = \;', ...
        mat2str(2 - ctrl0'), '\quad', sbtitle_suffix];
end

subtitle(ax, sbtitle, 'Interpreter', 'latex', 'FontSize', 10)

%% Save figure
filename_prefix = ['..', filesep, 'Plots', filesep, 'Raw', filesep, int2str(m), ...
    'sp__model_', int2str(model), '__log10c0_', num2str(log10(c0)), '__p1_', ...
    num2str(P(1), 3), '__E_', mat2str(E, 2)];
filename_suffix = ['__log10delta_', mat2str(log10(delta)), '__batch_#', ...
    int2str(batch)];

if model == 1
    filename = [filename_prefix, filename_suffix];
elseif model == 2 || model == 3
    filename = [filename_prefix, '__ctrl0_', mat2str(ctrl0'), ...
        filename_suffix];
end

saveas(gcf, [filename, '.fig'])
saveas(gcf, [filename, '.png'])

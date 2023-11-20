function plot__interbatch(rho, alpha, c0, P, E, delta, ctrl0, model, batches)
% Plots the inter-batch dynamics (values at dilutions) given by the
% solution for rho, alpha, and used parameter values, icluding relative
% species populations and their metabolic strategies. Saves the resulting
% figure in '../Plots/Raw/'.

%% Initialize
m = length(E); % Species no.
if m == 2
    cmap = ['r'; 'b'];
else
    cmap = colormap(parula(m));
end

%% Plot
figure
labels = strings(1, 2 * m);

for i = 1:m
    plot(0:batches, rho(i, :), 'Color', cmap(i, :), 'LineWidth', 2)
    labels(i) = ['$\rho_{', int2str(i), '} / \sum_\sigma \rho_\sigma$'];
    hold on
end
for i = 1:m
    alpha_i1 = reshape(alpha(i, 1, :), [1, length(alpha(i, 1, :))]);
    plot(0:batches, alpha_i1 / E(i), '--', 'Color', cmap(i, :))
    labels(i + m) = ['$\alpha_{', int2str(i), ',1} / E_{', int2str(i), '}$'];
    hold on
end

xlim([0, batches])
ylim([0, 1])
xlabel('Dilution number', 'Interpreter', 'latex')
ylabel('Population / Enzyme-1 fraction at dilution', 'Interpreter', 'latex')
legend(labels, 'Interpreter', 'latex')
title('Inter-batch dynamics', 'Interpreter', 'latex')
set(gca, 'Fontsize', 14)

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

subtitle(sbtitle, 'Interpreter', 'latex', 'FontSize', 11)

%% Save figure
filename_prefix = ['..', filesep, 'Plots', filesep, 'Raw', filesep, int2str(m), ...
    'sp__model_', int2str(model), '__log10c0_', num2str(log10(c0)), '__p1_', ...
    num2str(P(1), 3), '__E_', mat2str(E, 2)];
filename_suffix = ['__log10delta_', mat2str(log10(delta)), '__', ...
    int2str(batches), '_batches'];

if model == 1
    filename = [filename_prefix, filename_suffix];
elseif model == 2 || model == 3
    filename = [filename_prefix, '__ctrl0_', mat2str(ctrl0'), ...
        filename_suffix];
end

saveas(gcf, [filename, '.fig'])
saveas(gcf, [filename, '.png'])

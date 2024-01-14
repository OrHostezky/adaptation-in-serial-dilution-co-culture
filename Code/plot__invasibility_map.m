function plot__invasibility_map(p1s, log10deltas, is_invaded, params, rho_pre)
% Plots a heatmap displaying the invasibility character in a 2-species co-
% culture (index 1-3, see sim__invasibility_map() for details) vs nutrient
% profile bias (x) and adaptor's sensing tolerance (y). Saves the resulting
% figure in '../Plots/Raw/'.

%%% NOTE: Here, a co-culture of a non-adaptor (Sp. 1) and an adaptor
%%% (Sp. 2) is considered.

log10deltas_to_c0 = log10deltas - params.log10c0;
model = params.model;

%% Plot
figure

% Axes
ax = axes('Position', [.2, .17, .52, .67]);
if model == 1 || model == 2
    imagesc(2 * p1s - 1, log10deltas_to_c0, is_invaded)    
elseif model == 3
    imagesc(2 * p1s - 1, - log10deltas_to_c0, is_invaded)
end
set(gca, 'YDir', 'normal')

% Colormap
colormap(lines(3))
colorbar(ax, 'Position', [.73, .17, .03, .67], 'Ticks', 1:3, ...
    'TickLabels', {'Mutual invasibility','Adaptor invades', ...
    'Non-adaptor invades'}, 'TickLabelInterpreter', 'latex')
clim([0.5, 3.5])

%% Label
xlabel('$\frac{c_1(0) - c_2(0)}{c_0}$', 'Interpreter', 'latex', ...
    'FontSize', 20)
if model == 1 || model == 2
    ylabel('$\log_{10}{\left(\frac{\Delta}{c_0}\right)}$', ...
        'Interpreter', 'latex', 'FontSize', 16) 
elseif model == 3
    ylabel('$\log_{10}{\left(\frac{c_0}{\Delta}\right)}$', ...
        'Interpreter', 'latex', 'FontSize', 16)
end

sbtitle = ['$\textrm{Model} \;', int2str(model), '\quad \log_{10}(c_0) = \;', ...
    num2str(params.log10c0), '\quad E_\sigma = \;', mat2str(params.E, 2), ...
    '\quad \alpha_0 = \;', mat2str(params.alpha0)];
if model == 1
    sbtitle = [sbtitle, '\quad Ctrl_0 = \;', ...
    mat2str(params.ctrl0(params.is_adaptor == 1)), '$'];
elseif model == 2 || model == 3
    sbtitle = [sbtitle, '\quad c_{preferred} = \;', ...
    mat2str(2 - params.ctrl0(params.is_adaptor == 1)), '$'];
end

subtitle(sbtitle, 'Interpreter', 'latex', 'FontSize', 12)        
title('Adaptor VS Non-adaptor - Invasibility', ...
    'Interpreter', 'latex', 'FontSize', 18)

%% Save
filename = ['..', filesep, 'Plots', filesep, 'Raw', filesep, ...
    '1na_1a__heatmap__invasibility_vs_p1_log10delta__model_', ...
    int2str(model), '__log10c0_', num2str(params.log10c0), '__p1_', ...
    num2str(p1s(1)), 'to', num2str(p1s(end)), '__E_', mat2str(params.E, 2), ...
    '__alpha0_', mat2str(params.alpha0), '__ctrl0_', ...
    int2str(params.ctrl0(params.is_adaptor == 1)), '__log10(delta_to_c0)_', ...
    num2str(log10deltas_to_c0(1)), 'to', num2str(log10deltas_to_c0(end)), ...
    '__rho_pre_', num2str(rho_pre)];

saveas(gcf, [filename, '.fig'])
saveas(gcf, [filename, '.png'])

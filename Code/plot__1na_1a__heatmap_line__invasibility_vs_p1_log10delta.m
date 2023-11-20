function plot__1na_1a__heatmap_line__invasibility_vs_p1_log10delta(dir_name, ...
    ctrl0, log10c0s)
% Reads data from a set of 2-species co-culture, invasiblity-map simulations
% in '../Data/Raw/{dir_name}', and plots a 4-heatmap subplot-line, such that
% all maps correspond to adaptor's enzyme-production value in 'ctrl0' and
% total nutrient amounts in 'log10c0s' (MUST contain 4 values!). Each heatmap
% displays the invasibility character (index 1-3, see 'sim__invasibility_map'
% for details) vs nutrient profile bias (x) and adaptor's sensing tolerance
% (y). Saves the resulting figure in '../Plots/'.

%%% NOTE: Here, a co-culture of a non-adaptor (Sp. 1) and an adaptor (Sp. 2)
%%% is considered.

%% Initialize
figure('Renderer', 'painters', 'Position', 100 * [1, 1, 8, 2.2])
if ctrl0 == 1, bottom = .21; else, bottom = 0.37; end
[axs, ~] = tight_subplot(1, 4, [0, .04], [bottom, 12 / 22 - bottom], ...
    [.1, .28]); % File exchange

d = dir(['..', filesep, 'Data', filesep, 'Raw', filesep, dir_name, ...
    filesep, 'out*.mat']);

for dd = 1:length(d)
    %% Load
    fullname = [dir_name, filesep, d(dd).name];
    fullpath = ['..', filesep, 'Data', filesep, 'Raw', filesep, fullname];
    disp(['Reading ', fullname]);
    try
        load(fullpath, 'is_invaded', 'params', 'p1s', 'log10deltas', 'rho_pre');
    catch
        warning(['Failed loading ', fullname])
        continue
    end
    
    % Check if relevant data
    if params.ctrl0(params.isAdaptor == 1) ~= ctrl0 || ...
            ~ sum(params.log10c0 == log10c0s)
        continue
    end
    log10deltas_to_c0 = log10deltas - params.log10c0;
    model = params.model;

    %% Plot
    % Axes
    axes(axs(params.log10c0 == unique(log10c0s)))
    if model == 1 || model == 2
        imagesc(2 * p1s - 1, log10deltas_to_c0, is_invaded)    
    elseif model == 3
        imagesc(2 * p1s - 1, - log10deltas_to_c0, is_invaded)
    end
    
    % Colormap
    colormap(lines(3))
    if params.log10c0 == max(log10c0s)
        colorbar('Position', [.76, bottom, .02, 10 / 22], 'Ticks', 1:3, ...
            'TickLabels', {'Mutual invasibility', 'Adaptor invades', ...
            'Non-adaptor invades'}, 'TickLabelInterpreter', 'latex')
    end
    clim([0.5, 3.5])
    set(gca, 'YDir', 'normal', 'FontSize', 14)
    
    %% Label
    if params.log10c0 == min(log10c0s)
        if model == 1 || model == 2
            ylabel('$\log_{10}{\left(\frac{\Delta}{c_0}\right)}$', ...
                'Interpreter', 'latex', 'FontSize', 18)
        elseif model == 3
            ylabel('$\log_{10}{\left(\frac{c_0}{\Delta}\right)}$', ...
                'Interpreter', 'latex', 'FontSize', 18)
        end
    end

    if ctrl0 == 1
        subtitle(['$c_0 = 10^{', int2str(params.log10c0), '}$'], ...
            'Interpreter', 'latex', 'FontSize', 17)
    else
        xlabel('$\frac{c_1(0) - c_2(0)}{c_0}$', 'Interpreter', 'latex', ...
            'FontSize', 22)
    end
end

if model == 1
    stitle = ['$\textrm{Initially producing enzyme} \;', int2str(2 - ctrl0), '$'];
elseif model == 2 || model == 3
    stitle = ['$\textrm{Preferred nutrient} \;', int2str(2 - ctrl0), '$'];
end
sgtitle(stitle, 'Interpreter', 'latex', 'FontSize', 20)

%% Save
filename = ['..', filesep, 'Plots', filesep, ...
    '1na_1a__heatmap_line__invsibility_vs_p1_log10delta__model_', ...
    int2str(model), '__p1_', num2str(p1s(1)), 'to', num2str(p1s(end)), ...
    '__E_', mat2str(params.E,2), '__alpha0_', mat2str(params.alpha0), ...
    '__ctrl0_', mat2str(ctrl0), '__log10(delta_to_c0)_', ...
    num2str(log10deltas_to_c0(1)), 'to', num2str(log10deltas_to_c0(end)), ...
    '__rho_pre_', num2str(rho_pre)];

saveas(gcf, [filename, '.fig'])
saveas(gcf, [filename, '.png'])

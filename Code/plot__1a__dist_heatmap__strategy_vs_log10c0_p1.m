function plot__1a__dist_heatmap__strategy_vs_log10c0_p1(dir_name)
% Reads data from a set of fixed sensing tolerance, 1-adaptor full-dynamics
% simulations in '../Data/Raw/{dir_name}', and plots a heatmap, displaying
% the mean metabolic strategy (excluding transient) vs total nutrient
% amount (x), and nutrient profile bias (y), where in each x-y combination,
% the distribution of the strategy is plotted. Saves the resulting figure
% in '../Plots/'.

%%% NOTE: Uses v7.3 for saving data (can be a very large file!).

%% Initialize
% Load table
path = ['..', filesep, 'Data', filesep, 'Raw', filesep, dir_name];
try
    split_tab = readtable([path, filesep, 'to_run.csv']);
catch
    warning(["Failed loading 'to_run.csv' in ", path])
    return
end

% Set x-y variables
log10c0 = uniquetol(split_tab.log10c0, 1e-8);
log10c0_range = max(log10c0) - min(log10c0);
p1 = uniquetol(split_tab.p1, 1e-8);
p1_range = max(p1) - min(p1);

% Set general tile measures
strtpt1 = 0.11; endpt1 = 0.87;
strtpt2 = 0.1; endpt2 = 0.98;
wdth = (endpt1 - strtpt1) / length(log10c0);
hght = (endpt2 - strtpt2) / length(p1);

%% Plot
% Outer axes
figure
ax = axes('Position',[strtpt1, strtpt2, endpt1 - strtpt1, endpt2 - strtpt2], ...
    'Color', 'none');

xlim([min(log10c0) - 0.5 * (log10c0(2) - log10c0(1)), ...
    max(log10c0) + 0.5 * (log10c0(2) - log10c0(1))])
ylim(2 * [min(p1) - 0.5 * (p1(2) - p1(1)), max(p1) + 0.5 * (p1(2) - p1(1))] - 1)
ax.XTick = log10c0;
ax.YTick = 2 * p1 - 1;
xlabel('$\log_{10}(c_0)$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$\frac{c_1(0) - c_2(0)}{c_0}$', 'Interpreter', 'latex', 'FontSize', 20)

cb = colorbar(ax, 'Position', [endpt1 + 0.01, strtpt2, 0.03, endpt2 - strtpt2]);
cmap = flip(copper, 1);
colormap(cmap)
set(get(cb, 'label'), 'string', 'Mean Enzyme-1 fraction', ...
    'Interpreter', 'latex', 'FontSize', 16)

% Distributions
d = dir([path, filesep, 'out*.mat']);
for dd = 1:length(d)
    % Load
    disp(['Reading ', dir_name, filesep, d(dd).name]);
    try
       load([path, filesep, d(dd).name], 'output', 'params');
       range_strt = round(0.1 * size(output.rho, 2)); % Excluding transient
       if isfield(output,'t')
           t = output.t(range_strt:end);
       end
       strategy = output.alpha(1, 1, range_strt:end) / params.E(1);
       clear output
    catch
        warning(['Failed loading ', dir_name, filesep, d(dd).name])
        continue
    end
    
    % Create sub-axes
    left = strtpt1 + (params.log10c0 - min(log10c0)) / log10c0_range * ...
        (endpt1 - strtpt1 - wdth);
    bottom = strtpt2 + (params.P(1) - min(p1)) / p1_range * ...
        (endpt2 - strtpt2 - hght);
    axes('Position', [left, bottom, wdth, hght], 'Box', 'on');
    
    % Plot distribution
    cnts = histcounts(strategy, linspace(0, 1, 80));
    if exist('t', 'var')
        mean_strategy = trapz(t, strategy) / (t(end) - t(1));
    else
        mean_strategy = mean(strategy);
    end
    c = find(1 / 256 : 1 / 256 : 1 >= mean_strategy, 1);
    
    histogram('BinEdges', linspace(0, 1, 80),'BinCounts', cnts / max(cnts), ...
        'EdgeColor', 'none', 'FaceColor', cmap(c, :))
    set(gca, 'XTick', [])
    set(gca, 'YTick', [])
    set(gca, 'XTickLabel', {})
    set(gca, 'YTickLabel', {})
    xlim([0, 1])
    ylim([0, 1])
    set(gca, 'YScale', 'log')
end

%% Label and save
log10delta_to_c0 = params.log10delta(1) - params.log10c0;
model = params.model;

sbtitle = ['$ \textrm{Model} \;', int2str(model), "\quad E_\sigma' = \;", ...
    num2str(params.E(1), 2)];

filename_prefix = ['..', filesep, 'Plots', filesep, ...
        '1a__dist_heatmap__strategy_vs_log10c0_p1__model_', int2str(model) ...
        '__log10c0_', num2str(log10c0(1)), 'to', num2str(log10c0(end)), ...
        '__p1_', num2str(p1(1)), 'to', num2str(p1(end)), '__E_', ...
        num2str(params.E(1), 2)];
filename_suffix = ['__log10(delta_to_c0)_', num2str(log10delta_to_c0)];

if model == 1
    sbtitle = [sbtitle, '$'];
    filename = [filename_prefix, filename_suffix];
elseif model == 2 || model == 3
    sbtitle = [sbtitle, '\quad \textrm{Preferred Nutrient} \;', ...
        int2str(2 - params.ctrl0(1)), '$'];
    filename = [filename_prefix, '__ctrl0_', mat2str(params.ctrl0(1)), ...
        filename_suffix];
end

subtitle(ax, sbtitle, 'Interpreter', 'latex', 'FontSize', 16)        
%title(ax, "Adaptor's Mean Strategy distribution-heatmap", ...
%    'Interpreter', 'latex', 'FontSize', 18)

saveas(gcf, [filename '.png'])
hgsave(gcf, [filename '.fig'], '-v7.3')

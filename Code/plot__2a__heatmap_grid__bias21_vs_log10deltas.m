function plot__2a__heatmap_grid__bias21_vs_log10deltas(dir_name)
% Reads a collected-data table '../Data/Collected/collected_{dir_name}.csv'
% of an inter-batch simulation set in '../Data/Raw/{dir_name}' of 2 equal-
% budget, equal enzyme-production value adaptors, and plots a 3x3 grid of
% heatmaps - one for each combination of total nutrient amount and nutrient
% profile (there MUST be 3 * 3 = 9 combinations in the set!) -  each of
% which displaying the steady state (or other distinct long-time behavior)
% bias in population fractions, (rho*(1) - rho*(2)) / (rho*(1) + rho*(2)),
% vs adaptors' sensing tolerances (Adaptor-1 (x), Adaptor-2 (y)). Saves the
% resulting figure in '../Plots/'.

%% Initialize
col_tab = readtable(['..', filesep, 'Data', filesep, 'Collected', filesep, ...
    'collected__', dir_name, '.csv']);

% Constant parameters
model = col_tab.model;
E = [col_tab.E1(1), col_tab.E2(1)];

% Moving parameters
log10c0 = unique(col_tab.log10c0);
p1 = flip(unique(col_tab.p1));

% Heatmap variables
log10delta1_to_c0 = uniquetol(col_tab.log10delta1 - col_tab.log10c0, 1e-8);
log10delta2_to_c0 = uniquetol(col_tab.log10delta2 - col_tab.log10c0, 1e-8);

% Heatmap values
bias21 = nan(length(log10delta2_to_c0), length(log10delta1_to_c0));

% Long-time behavior categorization
imAlpha = ones(size(bias21));                                         % Non-finished runs
ss_0 = nan(length(log10delta1_to_c0) * length(log10delta2_to_c0), 2); % No detection
ss_3 = nan(length(log10delta1_to_c0) * length(log10delta2_to_c0), 2); % Moderate fluctuations
ss_4 = nan(length(log10delta1_to_c0) * length(log10delta2_to_c0), 2); % Large fluctuations

% Outer axes
figure('Renderer', 'painters', 'Position', 100 * [1, 1, 10, 5])
[axs, ~] = tight_subplot(length(p1), length(log10c0), [.01, .03], ...
    [.22, .215], [.39, .22]); % File exchange

for l = 1:length(p1)
    for k = 1:length(log10c0)
        %% Collect
        temptab = col_tab(col_tab.log10c0 == log10c0(k) & col_tab.p1 == p1(l), :);
        for row = 1:height(temptab)
            % Locate
            log10delta_row = [temptab.log10delta1(row), temptab.log10delta2(row)];
            [j, ~] = find(abs(log10delta1_to_c0 + log10c0(k) - log10delta_row(1)) < 1e-8);
            [i, ~] = find(abs(log10delta2_to_c0 + log10c0(k) - log10delta_row(2)) < 1e-8);
            
            % Categorize
            if temptab.steady_state(row) <= 0
                ss_0(find(isnan(ss_0(:, 1)), 1), :) = log10delta_row;
            elseif temptab.steady_state(row) == 4
                ss_4(find(isnan(ss_4(:, 1)), 1), :) = log10delta_row;
            else
                if temptab.steady_state(row) == 3
                    ss_3(find(isnan(ss_3(:, 1)), 1), :) = log10delta_row;
                end
                bias21(i,j) = temptab.pop_frac1(row) - temptab.pop_frac2(row);
            end
        end
        imAlpha(isnan(bias21)) = 0;
        
        %% Plot
        axes(axs(sub2ind([length(p1), length(log10c0)], k, l)))
        
        if model == 1 || model == 2
            imagesc(log10delta1_to_c0, log10delta2_to_c0, bias21, ...
                'AlphaData', imAlpha)
            hold on

            if sum(isnan(ss_0)) < numel(ss_0)
                plot(ss_0(:, 1) - log10c0(k), ss_0(:, 2) - log10c0(k), 'xk', ...
                    'MarkerSize', 8, 'HandleVisibility', 'off')
                % OR: 'DisplayName', 'No state detection'
            end
            if sum(isnan(ss_3)) < numel(ss_3)
                plot(ss_3(:, 1) - log10c0(k), ss_3(:, 2) - log10c0(k), '+k', ...
                    'MarkerSize', 8, 'DisplayName', 'Moderate fluctuations')
            end
            if sum(isnan(ss_4)) < numel(ss_4)
                plot(ss_4(:, 1) - log10c0(k), ss_4(:, 2) - log10c0(k), '*k', ...
                    'MarkerSize', 8, 'DisplayName', 'Large fluctuations')
            end
        
        elseif model == 3
            imagesc(- log10delta1_to_c0, - log10delta2_to_c0, bias21, ...
                'AlphaData', imAlpha)
            hold on

            if sum(isnan(ss_0)) < numel(ss_0)
                plot(log10c0(k) - ss_0(:, 1), log10c0(k) - ss_0(:, 2), 'xk', ...
                    'MarkerSize', 8, 'HandleVisibility', 'off')
            end
            if sum(isnan(ss_3)) < numel(ss_3)
                plot(log10c0(k) - ss_3(:, 1), log10c0(k) - ss_3(:, 2), '+k', ...
                    'MarkerSize', 8, 'DisplayName', 'Moderate fluctuations')
            end
            if sum(isnan(ss_4)) < numel(ss_4)
                plot(log10c0(k) - ss_4(:, 1), log10c0(k) - ss_4(:, 2), '*k', ...
                    'MarkerSize', 8, 'DisplayName', 'Large fluctuations')
            end
        end
        
        colormap(turbo)
        clim([-1, 1])
        set(gca, 'YDir', 'normal', 'Color', 0.9 * [1, 1, 1], 'FontSize', 20)

        %% Label
        if k == 1
            txt1 = ['$\frac{c_1(0) - c_2(0)}{c_0} = ', num2str(2 * p1(l) - 1), '$'];
        end

        if l == 1 
            if k == 1,     start = - 4;
            elseif k == 2, start = - 3.2;
            else,          start = - 2.7;
            end
            txt2 = ['$c_0 = 10^{', int2str(log10c0(k)), '}$'];
        end

        if model == 1 || model == 2
            ax.XTick = log10delta1_to_c0;
            ax.YTick = log10delta2_to_c0;

            if k == 1
                if l == 2
                    ylabel('$\log_{10}{\left(\frac{\Delta_2}{c_0}\right)}$', ...
                        'Interpreter', 'latex', 'FontSize', 28)
                end
                text(- 12.8, - 1.25, txt1, 'Interpreter', 'latex', 'FontSize', 30)
            end
            
            if k == 2 && l == length(p1)
                xlabel('$\log_{10}{\left(\frac{\Delta_1}{c_0}\right)}$', ...
                    'Interpreter', 'latex', 'FontSize', 28)
            end

            if l == 1
                text(start, 1, txt2, 'Interpreter', 'latex', 'FontSize', 30)
            end

        elseif model == 3
            ax.XTick = - log10delta1_to_c0;
            ax.YTick = - log10delta2_to_c0;

            if k == 1
                if l == 2
                    ylabel('$\log_{10}{\left(\frac{c_0}{\Delta_2}\right)}$', ...
                        'Interpreter', 'latex', 'FontSize', 28)
                end
                text(- 10.3, 1.25, txt1, 'Interpreter', 'latex', 'FontSize', 30)                
            end
            
            if k == 2 && l == length(p1)
                xlabel('$\log_{10}{\left(\frac{c_0}{\Delta_1}\right)}$', ...
                    'Interpreter', 'latex', 'FontSize', 28)
            end
            
            if l == 1
                text(start + 2.5, 3.5, txt2, 'Interpreter', 'latex', 'FontSize', 30)
            end
        end

        if l < length(p1), set(gca,'XTickLabel',{}); end
        if k > 1, set(gca,'YTickLabel',{}); end 
        
        %% Reset
        hold off
        bias21(:) = nan;
        imAlpha(:) = 1;
        ss_0(:) = nan;
        ss_3(:) = nan;
        ss_4(:) = nan;
    end
end

%% Label grid and save
lg = legend('FontSize', 21, 'Interpreter', 'latex');
lg.Position([1, 2]) = 5e-3 * [1, 1];

cb = colorbar(axs(1), 'Position', [.81, .22, .025, .565], 'Limits', [-1, 1]);
set(get(cb, 'label'), 'string', ['$\frac{\rho_1^* - \rho_2^*}', ...
    '{\rho_1^* + \rho_2^*}$'], 'Interpreter', 'latex', 'FontSize', 48)

filename_prefix = ['..', filesep, 'Plots', filesep, ...
    '2a__heatmap_grid__bias21_vs_log10deltas__model_', int2str(model), ...
    '__log10c0_', mat2str(log10c0), '__p1_', mat2str(flip(p1), 2), '__E_', ...
    mat2str(E)];
filename_suffix = ['__log10(delta_to_c0)_', num2str(log10delta1_to_c0(1)), ...
    'to', num2str(log10delta1_to_c0(end))];

if model == 1
    filename = [filename_prefix, filename_suffix];
elseif model == 2 || model == 3
    c_preferred = [col_tab.c_preferred1(1), col_tab.c_preferred2(1)];
    stitle = ['$\textrm{Model} \;', int2str(model), ...
        '\quad \textrm{Preferred nutrients:} \;', mat2str(c_preferred), '$'];
    sgtitle(stitle, 'Interpreter', 'latex', 'FontSize', 36)
    filename = [filename_prefix, '__ctrl0_', mat2str(2 - c_preferred), ...
        filename_suffix];
end

saveas(gcf, [filename, '.png'])
saveas(gcf, [filename, '.fig'])

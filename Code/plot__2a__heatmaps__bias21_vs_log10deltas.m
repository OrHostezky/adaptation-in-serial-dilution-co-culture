function plot__2a__heatmaps__bias21_vs_log10deltas(dir_name)
% Reads a collected-data table of a set of fixed nutrient amount, 2-adaptor
% co-culture, inter-batch simulations in dir_name, and plots a heatmap for
% each combination of nutrient profile and Adaptor-2's enzyme budget,
% displaying the steady state (or other distinct long-time behavior) bias
% in population fractions, (rho*(1) - rho*(2)) / sum(rho*), vs adaptors'
% sensing tolerances (Adaptor-1 (x), Adaptor-2 (y)). Saves the resulting
% figure in '../Plots/'.

%% Initialize
col_tab = readtable(['..', filesep, 'Data', filesep, 'Collected', filesep, ...
    'collected__', dir_name, '.csv']);

% Constant parameters
model = col_tab.model;
log10c0 = col_tab.log10c0(1);
E1 = col_tab.E1(1);

% Moving parameters
p1 = unique(col_tab.p1);
E2 = unique(col_tab.E2);

% Heatmap variables
log10delta1 = uniquetol(col_tab.log10delta1, 1e-8);
log10delta2 = uniquetol(col_tab.log10delta2, 1e-8);

% Heatmap values
bias21 = nan(length(log10delta2), length(log10delta1));

% Long-time behavior categorization
imAlpha = ones(size(bias21));                             % Non-finished runs
ss_0 = nan(length(log10delta1) * length(log10delta2), 2); % No detection
ss_3 = nan(length(log10delta1) * length(log10delta2), 2); % Moderate fluctuations
ss_4 = nan(length(log10delta1) * length(log10delta2), 2); % Large fluctuations

for k = 1:length(p1)
    for l = 1:length(E2)
        %% Collect
        temptab = col_tab(col_tab.p1 == p1(k) & col_tab.E2 == E2(l), :);
        for row = 1:height(temptab)
            % Locate
            log10delta_row = [temptab.log10delta1(row), ...
                temptab.log10delta2(row)];
            [j, ~] = find(abs(log10delta1 - log10delta_row(1)) < 1e-8);
            [i, ~] = find(abs(log10delta2 - log10delta_row(2)) < 1e-8);
            
            % Categorize
            if temptab.steady_state(row) <= 0
                ss_0(find(isnan(ss_0(:, 1)), 1), :) = log10delta_row;
            elseif temptab.steady_state(row) == 4
                ss_4(find(isnan(ss_4(:, 1)), 1), :) = log10delta_row;
            else
                if temptab.steady_state(row)==3
                    ss_3(find(isnan(ss_3(:,1)),1),:) = log10delta_row;
                end
                bias21(i, j) = temptab.pop_frac1(row) - temptab.pop_frac2(row);
            end
        end
        imAlpha(isnan(bias21)) = 0;
        
        %% Plot
        figure

        if model == 1 || model == 2
            imagesc(log10delta1 - log10c0, log10delta2 - log10c0, bias21, ...
                'AlphaData', imAlpha)
            hold on
    
            if sum(isnan(ss_0)) < numel(ss_0)
                plot(ss_0(:, 1) - log10c0, ss_0(:, 2) - log10c0, 'xk', ...
                    'MarkerSize', 50, 'DisplayName', 'No detection')
            end
            if sum(isnan(ss_3)) < numel(ss_3)
                plot(ss_3(:, 1) - log10c0, ss_3(:, 2) - log10c0, '+k', ...
                    'MarkerSize', 50, 'DisplayName', 'Moderate fluctuations')
            end
            if sum(isnan(ss_4)) < numel(ss_4)
                plot(ss_4(:, 1)- log10c0, ss_4(:, 2) - log10c0, '*k', ...
                    'MarkerSize', 50, 'DisplayName', 'Large fluctuations')
            end

        elseif model == 3
            imagesc(log10c0 - log10delta1, log10c0 - log10delta2, bias21, ...
                'AlphaData', imAlpha)
            hold on
    
            if sum(isnan(ss_0)) < numel(ss_0)
                plot(log10c0 - ss_0(:, 1), log10c0 - ss_0(:, 2), 'xk', ...
                    'MarkerSize', 50, 'DisplayName', 'No detection')
            end
            if sum(isnan(ss_3)) < numel(ss_3)
                plot(log10c0 - ss_3(:, 1), log10c0 - ss_3(:, 2), '+k', ...
                    'MarkerSize', 50, 'DisplayName', 'Moderate fluctuations')
            end
            if sum(isnan(ss_4)) < numel(ss_4)
                plot(log10c0 - ss_4(:, 1), log10c0 - ss_4(:, 2), '*k', ...
                    'MarkerSize', 50, 'DisplayName', 'Large fluctuations')
            end
        end
               
        %% Label
        colormap(turbo)
        h = colorbar;
        clim([-1, 1])
        set(gca, 'YDir', 'normal', 'Color', 0.9*[1, 1, 1])
        
        if model == 1 || model == 2
            xlabel('$\log_{10}{\left(\frac{\Delta_1}{c_0}\right)}$', ...
                'Interpreter', 'latex', 'FontSize', 14)
            ylabel('$\log_{10}{\left(\frac{\Delta_2}{c_0}\right)}$', ...
                'Interpreter', 'latex', 'FontSize', 14)
        elseif model == 3
            xlabel('$\log_{10}{\left(\frac{c_0}{\Delta_1}\right)}$', ...
                'Interpreter', 'latex', 'FontSize', 14)
            ylabel('$\log_{10}{\left(\frac{c_0}{\Delta_2}\right)}$', ...
                'Interpreter', 'latex', 'FontSize', 14)
        end

        set(get(h, 'label'), 'string', ['$\frac{\rho_1^* - \rho_2^*}', ...
            '{\rho_1^* + \rho_2^*}$'], 'Interpreter', 'latex', 'FontSize', 18)
        legend()
        
        %% Title and save
        sbtitle = ['$\textrm{Model} \;', int2str(model), ...
            '\quad \log_{10}(c_0) = \;', num2str(log10c0), ...
            '\quad c_1(0)/c_0 = \;', num2str(p1(k)), ...
            '\quad E_\sigma = \;', mat2str([E1, E2(l)])];
        
        filename_prefix = ['..', filesep, 'Plots', filesep, ...
            '2a__heatmap__bias21_vs_log10deltas__model_', int2str(model), ...
            '__log10c0_', num2str(log10c0), '__p1_', num2str(p1(k), 3), ...
            '__E_', mat2str([E1, E2(l)], 2)];
        filename_suffix = ['__log10deltas_', num2str(log10delta1(1)), 'to', ...
            num2str(log10delta1(end))];
        
        if model == 1
            sbtitle = [sbtitle, '$'];
            filename = [filename_prefix, filename_suffix];
        elseif model == 2 || model == 3
            c_preferred = [col_tab.c_preferred1(1), col_tab.c_preferred2(1)];
            sbtitle = [sbtitle, '\quad c_{preferred}^{\sigma} = \;', ...
                mat2str(c_preferred), '$'];
            filename = [filename_prefix, '__ctrl0', mat2str(2 - c_preferred), ...
                filename_suffix];
        end

        subtitle(sbtitle, 'Interpreter', 'latex', 'FontSize', 10)        
        title('Steady-state Populaion Bias VS Sensing Tolerances', ...
            'Interpreter', 'latex')

        saveas(gcf, [filename, '.fig'])
        saveas(gcf, [filename, '.png'])
        
        %% Reset
        hold off
        bias21(:) = nan;
        imAlpha(:) = 1;
        ss_0(:) = nan;
        ss_3(:) = nan;
        ss_4(:) = nan;
    end
end

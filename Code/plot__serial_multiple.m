function plot__serial_multiple(batches_to_plt)
% Plots the inter-batch and following intra-batch dynamics of all
% simulations in '../Data/Raw/out*.mat' data files (to execute, copy wanted
% simulations from their designated directories), and displays long-time
% population properties.
% It will optionally plot in addition the intra-batch dynamics of each of
% the batches in 'batches_to_plt', if given.

%%% NOTE: Use inter-batch data files created by 'sim__serial__interbatch'!

d = dir(['..', filesep, 'Data', filesep, 'Raw', filesep, 'out*.mat']);

for dd = 1:length(d)
    %% Load
    fullname = d(dd).name;
    fullpath = ['..', filesep, 'Data', filesep, 'Raw', filesep, fullname];
    disp(['Reading ', fullname]);
    try
       load(fullpath, 'output', 'params');
    catch
       warning(['Failed loading ', fullname]);
       continue
    end

    c0 = 10^params.log10c0;
    delta = 10.^params.log10delta;
    if isfield(params, 'ctrl0')
        ctrl0 = params.ctrl0;
    else
        ctrl0 = [];
    end 
    model = params.model;

    %% Inter-batch and following intra-batch
    if size(output.rho, 2) > 1
        plot__interbatch(output.rho, output.alpha, c0, params.P, params.E, ...
            delta, ctrl0, model, size(output.rho, 2) - 1)
    end

    if model == 1
        ctrl_b = output.ctrl(:, end);
    elseif model == 2 || model == 3
        ctrl_b = ctrl0;
    end
    [~, ~, ~, ~] = sim__batch(c0, params.K, params.P, params.rho0, ...
        output.rho(:, end), params.E, delta, output.alpha(:, :, end), ...
        ctrl_b, params.is_adaptor, model, 1, size(output.rho,2)); 

    %% Intra-batch of chosen batches
    if exist('batches_to_plt', 'var')
        for i = 1:length(batches_to_plt)
            if model == 1
                ctrl_b = output.ctrl(:, batches_to_plt(i));
            else
                ctrl_b = ctrl0;
            end
            [~, ~, ~, ~] = sim__batch(c0, params.K, params.P, params.rho0, ...
                output.rho(:, batches_to_plt(i)), params.E, delta, ...
                output.alpha(:, :, batches_to_plt(i)), ctrl_b, ...
                params.is_adaptor, model, 1, batches_to_plt(i)); 
        end
    end

    %% Steady-state population
    disp('End population profile is:')
    disp(string(output.pop_frac))
    
    if isfield(output,'ss')
        disp(['Steady-state type: ', int2str(output.ss), ...
            " (see 'app__simulations.m' for details)"])
    else
        disp("'output.ss' does not exist")
    end

    S = - sum(output.pop_frac .* log(output.pop_frac)); % Shannon entropy
    disp(['Effective species no.: ', num2str(exp(S))])
    
end

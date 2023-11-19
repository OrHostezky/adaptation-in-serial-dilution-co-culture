function collect_data__interbatch(dir_name)
% Collects wanted parameters and outputs of all finished inter-batch, 
% serial-dilution simulations in folder 'dir_name' and organizes data in a
% table, where each row corresponds to a single simulation. Saves the table
% in paths: '../Data/Raw/{dir_name}/collected.csv'
%           '../Data/Collected/collected__{dir_name}.csv'

%% Initialize
path = ['..', filesep, 'Data', filesep];
d = dir([path, 'Raw', filesep, dir_name, filesep, 'out*.mat']);
alltab = table;
row = 0;

%% Collect data to a table
for dd = 1:length(d)
    % Load
    fullname = [dir_name, filesep, d(dd).name];
    fullpath = ['..', filesep, 'Data', filesep, 'Raw', filesep, fullname];
    disp(['Reading ', fullname])
    try
        load(fullpath, 'output', 'params')
    catch
        warning(['Failed loading ', fullname])
        continue
    end
    
    if isfield(output, 'ss')
        % Collect
        tempstruct = struct;
        tempstruct.model = params.model;
        tempstruct.b0 = params.b0(1);
        tempstruct.log10c0 = params.log10c0;
        tempstruct.p1 = params.P(1);
        
        for i = 1:params.m
            name = ['E', num2str(i)];
            tempstruct.(name) = params.E(i);
        end
        if params.model == 2 || params.model == 3
            for i = 1:params.m
                name = ['c_preferred', num2str(i)];
                tempstruct.(name) = 2 - params.ctrl0(i);
            end
        end
        for i = 1:params.m
            name = ['log10delta', num2str(i)];
            tempstruct.(name) = params.log10delta(i);
        end
        for i = 1:params.m
            name = ['pop_frac', num2str(i)];
            tempstruct.(name) = output.pop_frac(i);
        end
        for i = 1:params.m
            name = ['strategy', num2str(i)];
            tempstruct.(name) = mean(output.alpha(i, 1, round(0.9 * ...
                end) : end)) / params.E(i); % Mean long-time strategy
        end
        
        tempstruct.batches = size(output.rho, 2) - 1;
        tempstruct.steady_state = output.ss;
        tempstruct.S = - sum((output.pop_frac) .* ...
            log(output.pop_frac)); % Shannon entropy
        %tempstruct.growth_int = output.growth_int;

        % Add to table
        temptab = struct2table(tempstruct);
        temptab.fullname{1} = fullname;
        row = row + 1;
        if row == 1
            alltab = temptab;
        else
            alltab = [alltab; temptab];
        end
    end
end

%% Write data tables
disp('--------------------------------------------')
outfile1 = [path, 'Raw', filesep, dir_name, filesep, 'collected.csv'];
writetable(sortrows(alltab, 'p1'), outfile1);
disp(['Finished collecting to ', outfile1]);

outfile2 = [path, 'Collected', filesep, 'collected__', dir_name, '.csv'];
writetable(sortrows(alltab, 'p1'), outfile2);
disp(['Finished collecting to ', outfile2]);
disp('--------------------------------------------')

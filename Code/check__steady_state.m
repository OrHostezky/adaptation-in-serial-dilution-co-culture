function unfinished = check__steady_state(dir_name)
% Returns rows in data table '../Data/Collected/collected_{dir_name}.csv'
% of inter-batch simulations that have not reached a steady-state of any
% sort (from the ones described in 'app_simulations.m'). If none exist,
% returns an empty array.

col_tab = readtable(['..', filesep, 'Data', filesep, 'Collected', ...
    filesep, 'collected__', dir_name, '.csv']);

unfinished = find(col_tab.steady_state == 0);
if isempty(unfinished)
    disp('No files found')
    unfinished = [];
end
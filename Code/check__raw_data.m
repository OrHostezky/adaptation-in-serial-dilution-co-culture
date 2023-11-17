function missing = check__raw_data(dir_name)
% Returns an array of missing raw-data files ('out_*.mat') in directory
% '../Data/Raw/{dir_name}', represented by their .mat-file number (*). If
% none are missing, returns an empty array.

%%% NOTE: There is an offset of +1 between the .mat-file number and the
%%% simulation number in associated 'to_run.csv'.

dir_path = ['..', filesep, 'Data', filesep, 'Raw', filesep, dir_name];
to_run = readtable([dir_path, filesep, 'to_run.csv']);
runs = height(to_run);
missing = [];

for i = 1:runs
    if ~ isfile([dir_path, filesep, 'out_', num2str(i + 1), '.mat'], 'file')
        missing = [missing, i + 1];
    end
end
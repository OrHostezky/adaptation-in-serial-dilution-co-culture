function split_runs__serial
% Defines a serial-dilution simulation set with varying parameter values
% between simulations in the set. Creates a directory in '../Data/Raw/' for
% the data, in which it puts a table 'to_run.csv' of moving-parameter
% values (each row corresponds to a single simulation), and a parameter
% template 'params.mat', into which moving parameters are assigned for each
% simulation accordingly.
% These simulations are to be applied by the wanted 'app__slurm__serial*'
% SLURM function, which allocates resources and executes the corresponding
% 'sim__serial*' function for all simulations in parallel.

%% Set parameters (use defaults for moving parameters)
%%% NOTE: For parameter description, see 'app__simulations'.

% Directory name (recommended form: '{SIMULATION_TYPE}__{param1}_{val1}__{...}')
dir_name = ['interbatch__2a__model_1__log10c0_[-2_0_2]__p1_[0.55_0.75_0.95]', ...
    '__E_[1_1]__ctrl0_[1_1]__log10(deltas_to_c0)_-2.5to0'];

% Nutrients
params.log10c0 = -2;
params.K = 1;
params.p = 2;
params.P = ones(params.p, 1) / params.p;

% Species
params.rho0 = 1;
params.m = 2;
params.b0 = zeros(params.m, 1) + 1 / params.m;

% Enzymatics
params.E = ones(1,params.m);
%params.alpha0 = params.E' / params.p;
%params.alpha0 = [params.alpha0, params.E' - params.alpha0];
params.ctrl0 = [0; 1];
params.log10delta = ones(1, params.m);
params.isAdaptor = ones(1, params.m);
params.model = 1;

% Control 
%%% NOTE: 'params.pltE', 'params.ss', and 'output.ss' are relevant only
%%% when executing the 'sim__serial__interbatch' function.
params.max_batches = 3e6;
params.plt = 0;
params.pltE = 0;
params.ss = 0;

%% Create a moving-parameters table
% Determine moving-parameter values
log10c0s = [2; 0; -2]; % Start from large c0, typically faster simulation
p1s = [0.55, 0.75, 0.95];
Es = [1, 1];
log10deltas = allcomb(-2.5:0.5:0, -2.5:0.5:0); % Using file exchange

% Create table
row = 0;
for c = 1:length(log10c0s)
    for p = 1:length(p1s)
        for e = 1:size(Es, 1)
            for d = 1:size(log10deltas, 1)
                % Plug values into a structure (vectors as column vectors)
                tempstruct = struct;
                tempstruct.log10c0 = log10c0s(c);
                tempstruct.p1 = p1s(p);
                tempstruct.E = string(mat2str(Es(e, :)'));
                tempstruct.log10delta = string(mat2str(log10c0s(c) + ...
                    log10deltas(d, :)')); % Delta / c0 dimensionless
                
                % Add to table
                temptab = struct2table(tempstruct);
                row = row + 1;
                if row == 1 
                    alltab = temptab;
                else
                    alltab = [alltab; temptab];
                end
            end
        end
    end
end

%% Output
disp('-------------------------------------------------------------')
outdir = ['..', filesep, 'Data', filesep, 'Raw', filesep, dir_name];

% Create parameter template
outparams = [outdir, filesep, 'params.mat'];
disp(['Saving template to ', outparams])
if ~ exist(outdir, 'dir')
    mkdir(outdir)
end
save(outparams, 'params')

% Write table
outtab = [outdir, filesep, 'to_run.csv'];
disp(['Saving table to ', outtab])
writetable(alltab, outtab)

disp('=============================================================')
disp('Done')
disp('=============================================================')

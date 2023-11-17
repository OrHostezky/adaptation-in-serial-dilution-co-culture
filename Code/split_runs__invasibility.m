function split_runs__invasibility
% Defines a set of invasibilty-map simulations with varying parameter
% values between simulations in the set. Creates a directory in
% '../Data/Raw/' for the data, in which it puts a table 'to_run.csv' of
% moving-parameter values (each row corresponds to a single simulation),
% and a parameter template 'params.mat', into which moving parameters are
% assigned for each simulation accordingly.
% These simulations are to be applied by the 'app__slurm__invasibility'
% SLURM function, which allocates resources and executes the
% 'sim__invasibility' function for all simulations in parallel.

%%% NOTE: Map ranges ('p1s', 'log10deltas') and initial pre-existing
%%% species population ('rho_pre') are defined inside the SLURM script.

%% Set parameters (use defaults for moving parameters)
%%% NOTE: For parameter description, see 'app__simulations'.

% Directory name (recommended form: '{SIMULATION_TYPE}__{param1}_{val1}__{...}')
dir_name = ['invasibility__1na_1a__model_1__log10c0_-2to2__p1_0.5to1__E_[1_1]__', ...
    'alpha0_[1_0;0.5_0.5]__ctrl0_[1_1]__log10(delta_to_c0)_-3to0__rho_pre_0.99999'];

% Nutrients
params.log10c0 = 2;
params.K = 1;
params.p = 2;
params.P = ones(params.p, 1) / params.p;

% Species
params.rho0 = 1;
params.m = 2;
params.b0 = zeros(params.m, 1) + 1 / params.m;

% Enzymatics
params.E = ones(1, params.m);
params.alpha0 = params.E' / params.p;
params.alpha0 = [params.alpha0, params.E' - params.alpha0];
params.ctrl0 = [1; 1];
params.log10delta = ones(1, params.m);
params.isAdaptor = [0, 1];
params.model = 1;

% Control
params.plt = 1;

%% Create a moving-parameters table
% Determine moving-parameter values
ctrl0s = [1; 1];
log10c0s = 2:-1:-2; % Start from large c0, typically faster simulation
Es = [1, 1];
alpha0s = [1, 0; 0.5, 0.5];

% Create table
rows = 0;
for c = 1:size(ctrl0s, 2)
    for cc = 1:length(log10c0s)
        for e = 1:size(Es, 1)
            for a = 1:size(alpha0s, 3)
                % Plug values into a structure (vectors as column vectors)
                tempstruct = struct;
                tempstruct.ctrl0 = string(mat2str(ctrl0s(:, c)));
                tempstruct.log10c0 = log10c0s(cc);
                tempstruct.E = string(mat2str(Es(e, :)'));
                tempstruct.alpha0 = string(mat2str(alpha0s(:, :, a)));

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


% This script directly applies the two main simulation types in this
% repository, and is mainly for execution of a limited number of simulations
% or gathering specific, immediate results. To apply large simulation
% sets, it is recommended to follow the Split-Apply-Combine protocol,
% using the .cmd (SLURM) execution files (see '../README.md').
% The simulations are based on a system of species that compete for
% resources (i.e. nutrients) in batches, where growth follows consumption.
% The framework here includes the possibility for species to adapt to
% nutrient concentrations (however, adaptation is strictly limited to the
% 2-nutrient case!), and a few models of adaptation (see below) are built
% in to the system, with the capacity to add more with relative ease (for
% fuller system details and explicit dynamics, see '../README.md').

% The two simulation types are:
% (1) A serial-dilution process: consecutive batches, where species are
% diluted at the end of each batch to start anew, while keeping the
% population profile at dilution. All batches start with the same supply
% (amount and profile), as provided by the user.
% The two sim__serial*() functions simulate this process: 'interbatch'
% suffix means that the state of the system is recorded only at dilutions,
% with the aim of characterizing the long-time behavior, whereas
% 'full_dynamics' suffix means that the whole time-course is recorded.
% (2) Invasibility tests: introducing a miniscule population of a species
% into a much larger population of another, and vice versa, and simulating
% a batch in both cases, with the aim of characterizing invasibility (see
% sim__invasibility_map()). Invasibility character is mapped across a
% domain of parameters. Currently, takes into account a non-adaptor vs
% adaptor co-culture.

% Each simulation will optionally save data and / or plot the recorded
% dynamics.

%%% NOTE: For a desription of the general structure of this repository's
%%% 'Code' section, see Script Index section in '../README.md' file. For a
%%% more specific script description, look at each script specifically.

clear; clc;

%% Set simulation parameters (+ parameter description)
% Data-saving path (recommended form: '../Data/Raw/{SIMULATION_TYPE}__{param1}_{val1}__{...}.mat')
outfile = ['../Data/Raw/interbatch__2a__model_1__log10c0_2__p1_0.6__E_[1_1]__', ...
    'log10(delta_to_c0)_[-1_-2]__ctrl0_[1_0].mat']; % 2a - 2 adaptors
%outfile = []; % Use for not saving data

% Nutrients
p1 = 0.75;                                     % Nutrient-1 fraction
params.log10c0 = 2;                            % Total amount (log)
params.K = 1;                                  % Monod constant [c0]
params.p = 2;                                  % Nutrient no. (MUST be 2 for adaptation!)
params.P = [p1; 1 - p1];                       % Profile

% Species
params.rho0 = 1;                               % Total amount [c0]
params.m = 2;                                  % Species no.
params.b0 = zeros(params.m, 1) + 1 / params.m; % Profile

% Enzymatics
params.E = [1, 1];                             % Budgets
params.alpha0 = [[1; 0] .* params.E', ...
    params.E' - [1; 0] .* params.E'];          % Initial strategies (row - #Sp., col - #Nut.) 
                                               % MUST hold:   sum(params.alpha0, 2) = params.E'

params.log10delta = params.log10c0 - [1, 2];   % Sensing tolerances [c0] (log)
params.ctrl0 = [1; 0];                         % Initial enzyme productions / nutrient preferences:
                                               % 1 (0) - Nutrient-1 (2))

params.is_adaptor = ones(1, params.m);         % 'Identities': 1 (0) for adaptor (non-adaptor)
params.model = 1;                              % The 3 adaptation regimes (i.e. 'models'): 
                                               % 1 - Relative-difference "bang-bang" adaptors
                                               % 2 - Preffered-nutrient "bang-bang" adaptors
                                               % 3 - Tolerance "Monod" adaptors (preferred nutrient)
                  
% Control
%%% NOTE: 'params.pltE', 'params.ss', and 'output.ss' are relevant only
%%% when executing the 'sim__serial__interbatch' function.
params.max_batches = 120; % Batches cut-off

params.plt = 1;  % 1 (0) to (not) plot results***
params.pltE = 0; % 1 (0) to (not) plot enzyme budgets
params.ss = 0;   % 1 (0) to (not) compute and save growth integrals
                 % Output (output.ss):
                 % ss = 0 - no steady-state has been reached / unclear (simulate more batches)
                 % ss = 1 - "exact" steady-state has been reached
                 % ss = 2 - convergence of running logarithmic averages has been reached, small
                 %          fluctuations (X < 0.02)
                 % ss = 3 - same as ss = 2, but moderate fluctuations (0.02 < X < 0.05)
                 % ss = 4 - same as ss = 2, but large fluctuations (X > 0.05)
                 % ss = - 1 - condition for ss = 2 mistakenly fulfilled, but at least one
                 %            population changing monotonically
                 % 'output.ss' regards the inter-batch steady-state

% *** 'sim__serial__interbatch' - inter-batch and following intra-batch dynamics
% *** 'sim__serial__full_dynamics' - full dynamics
% *** 'sim__invasibility_map' - invasibility map

%% Simulate
%%% NOTE: For several simulations, consider creating a string vector for
%%% 'outfile' with varying names (if needed).
sim_type = 1; 

% Serial dilution
if sim_type == 1 || sim_type == 2
    % Determine moving-parameter values
    log10c0s = 2; % Start from large c0, typically faster simulation
    p1s = .6;
    log10deltas_to_c0 = - [1, 2]; % Convenient to use Delta / c0

    for i = 1:size(log10deltas_to_c0, 1)
        for j = 1:length(log10c0s)
            for k = 1:length(p1s)
                % Plug moving-parameter values
                params.log10c0 = log10c0s(j);
                params.P = [p1s(k); 1 - p1s(k)];
                params.log10delta = log10c0s(j) + log10deltas_to_c0(i, :);
                

                % Record inter-batch dynamics
                if sim_type == 1
                    tic; output = sim__serial__interbatch(params, outfile); toc
                % Record full dynamics
                elseif sim_type == 2
                    tic; output = sim__serial__full_dynamics(params, outfile); toc
                end
            end
        end
    end

% Invasibility tests
elseif sim_type == 3
    % Map ranges
    p1s = 0.5:0.1:1;
    log10deltas = params.log10c0 + (-3:0);
    
    % Pre-existing species initial relative population
    rho_pre = 0.99999;
    
    is_invaded = sim__invasibility_map(params, outfile, p1s, ...
        log10deltas, rho_pre);
end

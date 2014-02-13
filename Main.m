
clear;

%% Model Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main data file name
% spec.dataName   = 'Data\data_fullsample.txt';
% spec.dataName   = 'Data\data_sh_20stations.csv';
% spec.dataName   = 'Data\logit_data_salvohuse.csv';
% spec.dataName   = 'Data\data_sh_full_cons.csv';
spec.dataName   = 'Data\data_spec1_sub_cleaned.csv';

% Share data file name
spec.shareName  = 'Data\data_share_full.txt';

% Log file name
[pathstr,name,ext] = fileparts(spec.dataName);
spec.logName    = ['Log\' name '.' datestr(now,'yyyymmdd.HHMM') '.log'];

% Number of consumer groups ( R )
n.conGroup  = 2;

% Number of product characteristic variables ( x_jl )
n.prodChar  = 0;

% Number of consumer characteristic variables ( x_i )
n.conChar   = 6;

% Allow for unobserved product heterogeneity ( xi_jl ) 
%   0 = no
%   1 = yes ( not implemented in this version !! )
spec.unobs  = 0;

%% Estimation Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Impose constraints to bound the diagonal elemets of the Choleski 
%   factor of the covariance matrix away from zeros. This makes sure the 
%   covariance remains positive definite during optimization.
%   1 = yes
%   2 = no
spec.constraint = 2;
spec.boundSize  = 1e-6;

% Base alternative
spec.base       = 1;
spec.scale      = 1 + (spec.base == 1); % not ready to change to other scale yet

% Number of random draws 
n.draw          = 10;

% Random draw type
%   1 = use pseudo-random draws
%   2 = use quasi-random draws - Halton sequence
spec.drawType   = 2;

% Halton sequnce specification
spec.halton.skip        = 1000;
spec.halton.leap        = 100;
spec.halton.scramble    = 'RR2';

% Seed of random draws ( matters only when pseudo-random draws are used )
rand( 'seed', 12345 );

% Solver ( NOTE: do not use KNITRO if frequency simulator is used )
%   1 = use KNITRO 
%   2 = use MATLAB fminsearch
%   3 = use MATLAB fmincon
spec.solver     = 1;

% Common optimization options
opt.maxIter     = 50000;
opt.maxFunEvals = 1e+10;
opt.display     = 'iter';
opt.tolFun      = 1e-8;
opt.tolCon      = 1e-8;
opt.tolX        = 1e-15;
opt.gradObj     = 'on';
opt.gradConstr  = 'off';

% KNITRO/fmincon specific optimization options
opt.algorithm   = 'interior-point';   % 'active-set' or 'interior-point'

%% Logging and timming
diary(spec.logName);
diary on;

display(spec);
display(n);
display(opt);

start_time = tic;
wall_clock = now;

%% Import Data and Construct Data Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConstructData;

%% Run Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diary( spec.logName );

% Start value
theta_0                         = ones(n.theta, 1 );
theta_0(1)                      = -10;

% Run estimation
RunEstimation;

%% Print Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PrintResults;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n\n   Total time      = %.4f seconds\n', toc(start_time));
fprintf(['   Wall-clock time = ' datestr(now - wall_clock,13) '\n']);
diary off;

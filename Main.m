
% clear all;

%% Model Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main data file name
% spec.dataName   = 'Data\data_fullsample.txt';
spec.dataName   = 'Data\sh_data_full.csv';

% Share data file name
spec.shareName  = 'Data\data_share_full.txt';

% Log file name
spec.logName    = 'Log\spec1_full_200.log';

% Number of consumer groups ( R )
n.conGroup  = 0;

% Number of product characteristic variables ( x_jl )
n.prodChar  = 0;

% Number of consumer characteristic variables ( x_i )
n.conChar   = 1;

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

% Number of random draws 
n.draw          = 100;

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
opt.gradObj     = 'off';
opt.gradConstr  = 'on';

% KNITRO/fmincon specific optimization options
opt.algorithm   = 'interior-point';   % 'active-set' or 'interior-point'

%% Import Data and Construct Data Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConstructData;

%% Run Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diary( spec.logName );

% Start value
theta_0                         = ones( n.theta, 1 );
theta_0(1)                      = -10;
theta_0( n.beta + 1 : n.theta ) = theta_0( n.beta + 1 : n.theta ) + ...
                                  linspace( 2, n.s + 1, n.s )';
% Run estimation
RunEstimation;

%% Print Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PrintResults;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diary off;

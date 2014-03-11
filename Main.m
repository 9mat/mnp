
clear;

%% Model Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main data file name
% spec.dataName   = 'Data\data_fullsample.txt';
% spec.dataName   = 'Data\data_sh_20stations.csv';
% spec.dataName   = 'Data\logit_data_salvohuse.csv';
spec.dataName   = 'Data\data_sh_full_spec1_cons.csv';
%spec.dataName   = 'Data\data_spec1_full.txt';

% Share data file name
spec.shareName  = 'Data\data_share_full.txt';

% Log file name
[~,name,~] = fileparts(spec.dataName);
spec.logName    = ['Log\' name '.' datestr(now,'yyyymmdd.HHMM') '.log'];

% Number of consumer groups ( R )
n.conGroup  = 2;

% Number of product characteristic variables ( x_jl )
n.prodChar  = 0;

% Number of consumer characteristic variables ( x_i )
n.conChar   = 10;

% for mfx
spec.paramType = [0;0;0;0;0;3;2;2;0;2*ones(n.conChar-2,1);1];

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
spec.boundSize  = 1e-9;

% Base alternative
spec.base       = 3;
spec.scale      = 1 + (spec.base == 1); % not ready to change to other scale yet

% Number of random draws 
n.draw          = 1000;

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
opt.gradConstr  = 'on';

% KNITRO/fmincon specific optimization options
opt.algorithm   = 'interior-point';   % 'active-set' or 'interior-point'

%% Log

%% Import Data and Construct Data Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConstructData;

%% Run Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare the log file and display the specifications
diary( spec.logName );
diary on;
display(n);
display(spec);
display(opt);
start_time = now;

% Start value
theta_0                         = rand(n.theta, 1 );
theta_0(1)                      = -10;

%% Run estimation
RunEstimation;

%% Print Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PrintResults;

fprintf(['\n\n\n Wall-clock running time = ' datestr(now - start_time,13) '\n']);

%% Save the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except thetaHat MLE spec n opt meanData paramType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% kk
[mfx, P, se_mfx, se_P] = marginalEffect(thetaHat, meanData, n, spec, MLE.cov);
fprintf('\n\nMarginal effects at means and SE\n');
display([mfx sqrt(diag(se_mfx))]);
fprintf('\n\nChoice Probability at means and Se\n'); 
display([P sqrt(diag(se_P))]);

[~,name,~] = fileparts(spec.logName);
save(['Results/' name '.mat']);

diary off;

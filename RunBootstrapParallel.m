function RunBootstrapParallel(bstrind)


%% Model Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main data file name
% spec.dataName   = 'Data\data_fullsample.txt';
% spec.dataName   = 'Data\data_sh_20stations.csv';
% spec.dataName   = 'Data\logit_data_salvohuse.csv';
%spec.dataName   = 'Data\data_sh_full_cleaned_keepchoice1.csv';
spec.dataName   = 'Data\data_new_spec4_2.csv';
spec.include_emidgrade = true;
spec.keep_treattype = -1;

% Log file name
[~,name,ext] = fileparts(spec.dataName);
name2 = [name '.' datestr(now,'yyyymmdd.HHMM')];
spec.logName    = ['Log\' name2 '.log'];

% Share data file name
spec.shareName  = ['Data\share_' name ext];

% Number of consumer groups ( R )
n.conGroup  = 0;

% Number of product characteristic variables ( x_jl )
n.prodChar  = 0;

% Number of consumer characteristic variables ( x_i )
n.conChar   = 8;

% for mfx
% spec.paramType  = [0;0;0;0;0;3;2;2;2;2;2;2;2;2];
% spec.paramGroup = [0;0;0;0;0;1;2;3;3;3;4;4;5;6];
% spec.paramID    = [0;0;0;0;0;1;2;3;4;5;6;7;8;9];

% data_new spec6
% spec.paramType  = [0;0;0;0;0;3;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2];
% spec.paramGroup = [0;0;0;0;0;1;2;3;3;3;4;4;5;6;2;3;3;3;4;4;5;6];
% spec.paramID    = [0;0;0;0;0;1;2;3;4;5;6;7;8;9;2;3;4;5;6;7;8;9];

% data_new spec4_2
spec.paramType  = [0;0;0;0;0;3;2;2;2;2;2;2;2;2];
spec.paramGroup = [0;0;0;0;0;1;2;3;3;3;4;4;5;6];
spec.paramID    = [0;0;0;0;0;1;2;3;4;5;6;7;8;9];

% data_new treatspec1
% spec.paramType = [0;0;0;0;0;3;2;2;2;2;2;2;2;2];
% spec.paramGroup = [0;0;0;0;0;1;3;4;4;4;5;5;6;7];
% spec.paramID = [0;0;0;0;0;1;6;7;8;9;10;11;12;13];

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
spec.constraint = 1;
spec.boundSize  = 1e-6;

% Base alternative
spec.base       = 3;
spec.scale      = 1 + (spec.base == 1); % not ready to change to other scale yet

% Number of random draws 
n.draw          = 10;
n.mfxdraw       = n.draw;

% Random draw type
%   1 = use pseudo-random draws
%   2 = use quasi-random draws - Halton sequence
spec.drawType   = 2;

% Halton sequnce specification
spec.halton.skip        = 1000;
spec.halton.leap        = 100;
spec.halton.scramble    = 'RR2';

% Seed of random draws ( matters only when pseudo-random draws are used )
rng('default');

% Solver ( NOTE: do not use KNITRO if frequency simulator is used )
%   1 = use KNITRO 
%   2 = use MATLAB fminsearch
%   3 = use MATLAB fmincon
spec.solver     = 3;

% Common optimization options
opt.maxIter     = 3;
opt.maxFunEvals = 1e+10;
opt.display     = 'iter';
opt.tolFun      = 1e-8;
opt.tolCon      = 1e-8;
opt.tolX        = 1e-18;
opt.gradObj     = 'on';
opt.gradConstr  = 'on';

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

thetaHats = bootstrapparallel( n, spec, opt, bstrind );
result = [bstrind; thetaHats];

save(['Results/`thetaHats ' num2str(bstrind) '.mat'],'result');

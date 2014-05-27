clear;

%% Model Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main data file name
% spec.dataName   = 'Data/data_fullsample.txt';
% spec.dataName   = 'Data/data_sh_20stations.csv';
% spec.dataName   = 'Data/logit_data_salvohuse.csv';
spec.dataName   = 'Data/data_new.csv';
%spec.dataName   = 'Data/data_spec1_full.txt';

% Share data file name
spec.shareName  = 'Data/data_share_full.txt';

% Log file name
[~,name,~] = fileparts(spec.dataName);
spec.logName    = ['Log/' name '.' datestr(now,'yyyymmdd.HHMM') '.log'];

% Number of consumer groups ( R )
n.conGroup  = 0;

% Number of product characteristic variables ( x_jl )
n.prodChar  = 0;

% Number of consumer characteristic variables ( x_i )
n.conChar   = 14;

% for mfx

% data_sh_full_cons.csv - 10X
% spec.paramType  = [0;0;0;0;0;3;0;2;2;2;2;2;2;2;2;1];
% spec.paramGroup = [0;0;0;0;0;1;0;2;3;3;3;4;4;5;6;7];
% spec.paramID    = [0;0;0;0;0;1;0;2;3;4;5;6;7;8;9;10];

% data_sh_full_23X.csv
% spec.paramType  = [0;0;0;0;0;3;0;2;2;2;2;2;2;2;2;1;2;2;2;2;2;2;2;2;2;2;2;2];
% spec.paramGroup = [0;0;0;0;0;1;0;2;3;3;3;4;4;5;6;7;8;8;8;8;8;9;9;9;9;10;10;10];
% spec.paramID    = [0;0;0;0;0;1;0;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22];
% 
% nX = n.conGroup + n.conChar;
% spec.paramType  = spec.paramType(1:nX+6);
% spec.paramGroup = spec.paramGroup(1:nX+6);
% spec.paramID    = spec.paramID(1:nX+6);

% data_sh_full_23X_dup.csv - 8 Groups, 22X
% spec.paramType  = [0;0;0;0;0;3;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;1;2;2;2;2;2;2;2;2;2;2;2;2;0];
% spec.paramGroup = [0;0;0;0;0;1;2;3;3;3;4;4;5;6;2;3;3;3;4;4;5;6;7;8;8;8;8;8;9;9;9;9;10;10;10;0];
% spec.paramID    = [0;0;0;0;0;1;2;3;4;5;6;7;8;9;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;0];

% data_sh_full_cons_double.csv
% spec.paramType  = [0;0;0;0;0;3;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;1;0];
% spec.paramGroup = [1;0;0;0;0;1;2;3;3;3;4;4;5;6;2;3;3;3;4;4;5;6;7;0];
% spec.paramID    = [2;0;0;0;0;1;2;3;4;5;6;7;8;9;2;3;4;5;6;7;8;9;10;0];

% data_new.csv - 14X
spec.paramType = [0;0;0;0;0;3;2;2;2;2;2;2;2;2;1;2;2;2;2;0];
spec.paramGroup = [0;0;0;0;0;1;2;3;3;3;4;4;5;6;7;8;8;8;8;0];
spec.paramID = [0;0;0;0;0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;0];

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
rand( 'seed', 12345 );

% Solver ( NOTE: do not use KNITRO if frequency simulator is used )
%   1 = use KNITRO 
%   2 = use MATLAB fminsearch
%   3 = use MATLAB fmincon
spec.solver     = 4;

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

[~,name,~] = fileparts(spec.logName);
save(['Settings/setting_' name '.mat']);

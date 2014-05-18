
clear;

%% Model Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main data file name
% spec.dataName   = 'Data\data_fullsample.txt';
% spec.dataName   = 'Data\data_sh_20stations.csv';
% spec.dataName   = 'Data\logit_data_salvohuse.csv';
spec.dataName   = 'Data\data_sh_full_cleaned.csv';
%spec.dataName   = 'Data\data_spec1_full_cleaned.csv';


% Log file name
[~,name,ext] = fileparts(spec.dataName);
spec.logName    = ['Log\' name '.' datestr(now,'yyyymmdd.HHMM') '.log'];

% Share data file name
spec.shareName  = ['Data\share_' name ext];

% Number of consumer groups ( R )
n.conGroup  = 0;

% Number of product characteristic variables ( x_jl )
n.prodChar  = 0;

% Number of consumer characteristic variables ( x_i )
n.conChar   = 9;

% for mfx
spec.paramType = [0;0;0;0;0;3;2*ones(n.conChar-1,1);1];

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
rng('default');

% Solver ( NOTE: do not use KNITRO if frequency simulator is used )
%   1 = use KNITRO 
%   2 = use MATLAB fminsearch
%   3 = use MATLAB fmincon
spec.solver     = 1;

% Common optimization options
opt.maxIter     = 2000;
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

%% Import Data and Construct Data Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConstructData;

%% Run Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diary( spec.logName );

% Start value
% load theta_00.mat;
% theta_0 = thetaHat;
theta_0                         = rand(n.theta, 1 );
theta_0(1)                      = -20;

%% Run estimation
RunEstimation;

%% Print Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PrintResults;

fprintf('\n\n   Total time      = %.4f seconds\n', toc(start_time));
fprintf(['   Wall-clock time = ' datestr(now - wall_clock,13) '\n']);

%% Save the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except dataHeader dataMatrix dataR choicesetcode thetaHat MLE spec n opt meanData;
[~,name,~] = fileparts(spec.logName);
%tic;[mfx, P, se_mfx, se_P] = marginalEffect(thetaHat, meanData, n, spec, MLE.cov);toc;
tic;[ amfx, se_amfx, aP, se_aP ] = AME( dataMatrix, dataR, choicesetcode, thetaHat, n, spec, MLE.cov ); toc;
%%
mfxHeader = {};
for i = 1:numel(spec.paramType)
    if spec.paramType(i) == 0; continue; end;
    if spec.paramType(i) < 3
        mfxHeader{end+1} = dataHeader{i};
    else
        for j=1:n.maxChoice
            mfxHeader{end+1} = [dataHeader{i} num2str(j)];
        end
    end
end

mfxHeader = repmat(mfxHeader, 1, n.maxChoice);
%%
%printmat([mfx, se_mfx, abs(mfx./se_mfx)], 'Marginal Effects at Means', strjoin(mfxHeader), 'MEM SE t');
printmat([amfx, se_amfx, abs(amfx./se_amfx)], 'Average Marginal Effects', strjoin(mfxHeader), 'AME SE t');
 save(['Results/' name '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diary off;

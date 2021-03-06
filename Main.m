
clear;

%% Model Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main data file name
% spec.dataName   = 'Data\data_fullsample.txt';
% spec.dataName   = 'Data\data_sh_20stations.csv';
% spec.dataName   = 'Data\logit_data_salvohuse.csv';
spec.dataName   = 'Data\data_sh_full_cons.csv';
%spec.dataName   = 'Data\data_spec1_full.txt';

% Share data file name
spec.shareName  = 'Data\data_share_full.txt';

% Log file name
[~,name,~] = fileparts(spec.dataName);
spec.logName    = ['Log\' name '.' datestr(now,'yyyymmdd.HHMM') '.log'];

% Number of consumer groups ( R )
n.conGroup  = 0;

% Number of product characteristic variables ( x_jl )
n.prodChar  = 0;

% Number of consumer characteristic variables ( x_i )
n.conChar   = 10;

% for mfx
spec.paramType  = [0;0;0;0;0;3;0;2;2;2;2;2;2;2;2;1];
spec.paramGroup = [0;0;0;0;0;1;0;2;3;3;3;4;4;5;6;7];
spec.paramID    = [0;0;0;0;0;1;0;2;3;4;5;6;7;8;9;10];

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
opt.algorithm   = 'active-set';   % 'active-set' or 'interior-point'

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
theta_0                         = zeros(n.theta, 1 );
theta_0(1)                      = -10;
theta_0(end) = 1; theta_0(end-1) = 0.5;

%% Run estimation
RunEstimation;

%% Print Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PrintResults;

fprintf(['\n\n\n Wall-clock running time = ' datestr(now - start_time,13) '\n']);

%% Save the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except dataMatrix dataR dataHeader choicesetcode thetaHat MLE spec n opt meanData paramType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% kk
tic;[mfx, P, se_mfx, se_P] = marginalEffect(thetaHat, meanData, n, spec, MLE.cov);toc;
tic;[ amfx, se_amfx, aP, se_aP ] = AME( dataMatrix, dataR, choicesetcode, thetaHat, n, spec, MLE.cov ); toc;
%%
mfxHeader = {};
type = mod(spec.paramType, 10);
for i = 1:numel(type)
    if type(i) == 0; continue; end;
    if type(i) < 3
        mfxHeader{end+1} = dataHeader{i};
    else
        for j=1:n.maxChoice
            mfxHeader{end+1} = [dataHeader{i} num2str(j)];
        end
    end
end

mfxHeader = repmat(mfxHeader, 1, n.maxChoice);
printmat([mfx, se_mfx, abs(mfx./se_mfx)], 'Marginal Effects at Means', strjoin(mfxHeader), 'MEM SE t');
printmat([amfx, se_amfx, abs(amfx./se_amfx)], 'Average Marginal Effects', strjoin(mfxHeader), 'AME SE t');

[~,name,~] = fileparts(spec.logName);
save(['Results/' name '.mat']);

diary off;

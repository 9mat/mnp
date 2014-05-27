function run(settingfile)
%RUN Summary of this function goes here
%   Detailed explanation goes here


load(settingfile);
ConstructData;

%% Start value
theta_0         = zeros(n.theta, 1);
theta_0(1)      = -10;
theta_0(end)    = 1; 
theta_0(end-1)  = 0.5;


%% Run Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare the log file and display the specifications
diary( spec.logName );
diary on;

display(n);
display(spec);
display(opt);

start_time = now;


%% Run estimation

if spec.solver == 4 && exist('ipopt','file') == 0
    addpath('/home/svu/a0034101/ipopt');
end

RunEstimation;

%% Print Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PrintResults;

fprintf(['\n\n\n Wall-clock running time = ' datestr(now - start_time,13) '\n']);

%% Save the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except dataMatrix dataR dataHeader choicesetcode thetaHat MLE spec n opt meanData paramType;
[~,name,~] = fileparts(spec.logName);
save(['Results/' name '.mat']);

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
joinHeader = [];
for i=1:numel(mfxHeader)
    joinHeader = [joinHeader mfxHeader{i} ' '];
end

printmat([mfx, se_mfx, abs(mfx./se_mfx)], 'Marginal Effects at Means', joinHeader, 'MEM SE t');
printmat([P, se_P], 'Choice Probability at Means', num2str(1:n.maxChoice), 'Prob se');

printmat([amfx, se_amfx, abs(amfx./se_amfx)], 'Average Marginal Effects', joinHeader, 'AME SE t');
printmat([aP, se_aP], 'Average Choice Probability', num2str(1:n.maxChoice), 'Prob se');

save(['Results/' name '.mat']);
diary off;

end


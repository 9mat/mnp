function run_AME(resultfile)
load(resultfile);

[~,name,~] = fileparts(spec.logName);
diary(['Log/AME_' name '.log']);


tic;[ amfx, se_amfx, aP, se_aP ] = AME( dataMatrix, dataR, choicesetcode, thetaHat, n, spec, MLE.cov ); toc;
%%
mfxHeader = {};
for i = 1:numel(spec.paramType)
    type = mod(spec.paramType(i), 10);
    if type == 0; continue; end;
    if type < 3
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
printmat([aP, se_aP], 'Average Choice Probability', num2str(1:n.maxChoice), 'Prob se');
save(['Results/AME_' name '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diary off;

end
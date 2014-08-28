function [ thetaHats ] = bootstrap( n, spec, opt, nbstr )
%BOOTSTRAP Summary of this function goes here
%   Detailed explanation goes here

%% Import Data and Construct Data Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Data
tempData        = importdata( spec.dataName );
dataMatrix      = tempData.data;
dataHeader      = tempData.colheaders;
clear tempData

% Share data
tempData        = importdata( spec.shareName );
shareMatrix     = tempData.data;
clear tempData;

[dataR, dataS, shareHat, identifiable, paramsid, n] = ConstructData(dataMatrix, shareMatrix, n, spec);

%% Run Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_0                         = zeros(n.theta, 1 );
theta_0(1)                      = -10;
theta_0(end) = 1;
theta_0(end-4:end) = [0;0;1;0;1];

%% Run estimation
fprintf('Bootstrap: Estimate with original data\n');

tic;
thetaHat = RunEstimation(dataR, dataS, theta_0, n, spec, opt, shareHat, identifiable, paramsid);
toc

keep = [(1:n.conGroup+1)';paramsid.beta_1(identifiable.beta_1);paramsid.beta_2(identifiable.beta_2)];
thetaHats = zeros(numel(keep), nbstr);
thetaHats(:,1) = thetaHat(keep);

clear dataR dataS;

clusters = unique(dataMatrix(:,1));
maxConID = max(dataMatrix(:,2))+1;

for ibstr = 1:nbstr-1
    fprintf('Bootstrap iter %d\n', ibstr);
    index = randi(numel(clusters), 1, numel(clusters));
    bstr_sample = clusters(index);
    bstr_dataMatrix = [];
    bstr_shareMatrix = [];
    for i = 1:numel(bstr_sample)
        bstr_newData = dataMatrix(dataMatrix(:,1) == bstr_sample(i),:);
        bstr_newShare = shareMatrix(shareMatrix(:,1) == bstr_sample(i),:);
        
        bstr_newData(:,1) = i;
        bstr_newData(:,2) = bstr_newData(:,2) + i*maxConID;
        bstr_newShare(:,1) = i;
        
        bstr_dataMatrix = [bstr_dataMatrix; bstr_newData];
        bstr_shareMatrix = [bstr_shareMatrix; bstr_newShare];
    end
    
    [dataR, dataS, shareHat, identifiable, paramsid, n] = ConstructData(bstr_dataMatrix, bstr_shareMatrix, n, spec);
    
    theta_0       = zeros(n.theta, 1 );
    theta_0(keep) = thetaHat(keep);
    theta_0(end-4:end) = thetaHat(end-4:end);
    
    bstr_thetaHat = RunEstimation(dataR, dataS, theta_0, n, spec, opt, shareHat, identifiable, paramsid);
    
    thetaHats(:,ibstr+1) = bstr_thetaHat(keep);
    toc
end

end


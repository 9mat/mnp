function [ amfx, se_amfx, P, se_P ] = AME( dataMatrix, dataR, choicesetcode, thetaHat, n, spec, cov )
%AME Summary of this function goes here
%   Detailed explanation goes here


conID = sort(unique(dataMatrix(:,2)));
mfx = zeros(n.mfx*n.maxChoice, n.con);
P = zeros(n.maxChoice, n.con);
d_mfx_d_theta = zeros(n.mfx*n.maxChoice, numel(thetaHat), n.con);
d_P_d_theta = zeros(n.maxChoice, numel(thetaHat), n.con);

uniquecode = sort(unique(choicesetcode));
nchoiceset = numel(unique(choicesetcode));
for i = 1:numel(conID)
    if mod(i,100) == 0; fprintf('   Calculating mfx for person %4d\n',i); end
    index = conID(i) == dataMatrix(:,2);
    k = find(choicesetcode(find(index,1,'first'))==uniquecode);
    m = find(sort(unique(dataMatrix(choicesetcode == uniquecode(k), 1))) == dataMatrix(find(index,1,'first'),1));
    pick = dataR{k}.pick;
    nn = dataR{k}.n;
    todelete = 1 + nn.conGroup + (nn.maxChoice-1)*(nn.prodChar+nn.conChar) + [1:(nn.maxChoice-1)*(m -1 ), ...
        (nn.maxChoice-1)*m + 1: (nn.maxChoice-1)*nn.market];
    pick(todelete) = [];
    [ mfx_i, P_i, ~, ~, d_mfx_d_theta_i, d_P_d_theta_i] = ...
        marginalEffect(thetaHat(pick), dataMatrix(index,:), n, spec, cov(pick,pick));
    mfx(:, i) = mfx_i;
    P(:, i) = P_i;
    d_mfx_d_theta(:, pick, i)  = d_mfx_d_theta_i;
    d_P_d_theta(:, pick, i) = d_P_d_theta_i;
end


amfx = mean(mfx,2);
d_amfx_d_theta = mean(d_mfx_d_theta, 3);
se_amfx = sqrt(diag(d_amfx_d_theta*cov*d_amfx_d_theta'));

P = mean(P, 2);
d_P_d_theta = mean(d_P_d_theta, 3);
se_P = sqrt(diag(d_P_d_theta*cov*d_P_d_theta'));
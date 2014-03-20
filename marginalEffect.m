function [ mfx_full, P_full, se_mfx, se_P, d_mfx_d_theta_full, d_P_d_theta_full] = ...
    marginalEffect(thetaHat, meanData, n, spec, V_theta)
%MFX Summary of this function goes here
%   Detailed explanation goes here
%   paramType = 0 --> no marginal effect
%   paramType = 1 --> continuous and product independent
%   paramType = 2 --> binary and product independent
%   paramType = 3 --> continuous and product dependent
%   paramType = 4 --> binary and product dependent

choiceset = sort(unique(meanData(:,4)));
nchoice = numel(choiceset);
ncon = numel(unique(meanData(:,2)));
nmfx = sum(spec.paramType > 0) + (nchoice-1)*sum(spec.paramType > 2);

data1 = repmat(meanData,nmfx,1);
data2 = data1;
epsilon = 1e-5;
dx = [];

index = 1:size(meanData,1);
for i = 1:numel(spec.paramType)
    if spec.paramType(i) == 0; continue; end
    
    if spec.paramType(i) == 1 || spec.paramType(i) == 2
        if spec.paramType(i) == 1
            data1(index,i) = meanData(:,i)*(1-epsilon);
            data2(index,i) = meanData(:,i)*(1+epsilon);
        else
            data1(index,i) = 0;
            data2(index,i) = 1;
        end
        
        dx = [dx;data2(index(1:nchoice:end),i) - data1(index(1:nchoice:end),i)];
        index = index + nchoice*ncon;
        continue;
    end
    
    for j=1:nchoice
        if spec.paramType(i) == 3
            data1(index(j:nchoice:end),i) = meanData(j:nchoice:end,i)*(1-epsilon);
            data2(index(j:nchoice:end),i) = meanData(j:nchoice:end,i)*(1+epsilon);
        else
            data1(index(j:nchoice:end),i) = 0;
            data2(index(j:nchoice:end),i) = 1;
        end
        
        dx = [dx;data2(index(j:nchoice:end),i) - data1(index(j:nchoice:end),i)];
        index = index + nchoice*ncon;
    end
end

m = size(data1, 1);
data1 = repmat(data1, nchoice, 1);
data2 = repmat(data2, nchoice, 1);
data3 = repmat(meanData, nchoice, 1);
dx = repmat(dx(:), nchoice, 1);
for j = 1:nchoice
    data1((j-1)*m+1:j*m,5) = choiceset(j);
    data2((j-1)*m+1:j*m,5) = choiceset(j);
    data3((j-1)*nchoice+1:j*nchoice,5) = choiceset(j);
end

data = [data1; data2; data3];
conID = repmat(1:2*m+nchoice*ncon, nchoice, 1);
data(:,2) = conID(:);

n.con = numel(unique(data(:,2)));
n.draw = 10;
[dataR, ~] = ConstructDataGroup(data, n, spec);
dataR.draw.uni = repmat(dataR.draw.uni(1,:,:), [n.con 1 1]);
P = ProbitProb(thetaHat, dataR, dataR.n, dataR.spec);

m = (n.con - nchoice*ncon)/2;
P1 = P(1:m);
P2 = P(m+1:2*m);

mfx = (P2-P1)./dx;
P = P(2*m+1:end);


epsilon = 1e-5;
k = numel(thetaHat);
thetaHatPlus = repmat(thetaHat, 1, k).*(1 + epsilon*eye(k));
thetaHatMinus = repmat(thetaHat, 1, k).*(1 - epsilon*eye(k));
d_theta = diag(thetaHatPlus) - diag(thetaHatMinus);

P_plus = zeros(n.con, k);
P_minus = zeros(n.con, k);
for i = 1:k
    P_plus(:,i) = ProbitProb(thetaHatPlus(:,i), dataR, dataR.n, dataR.spec);
    P_minus(:,i) = ProbitProb(thetaHatMinus(:,i), dataR, dataR.n, dataR.spec);
end

d_P_d_theta = bsxfun(@rdivide,P_plus - P_minus, d_theta');

d_P2_d_theta = d_P_d_theta(1:m,:);
d_P1_d_theta = d_P_d_theta(m+1:2*m,:);
d_P_d_theta = d_P_d_theta(2*m+1:end,:);
d_mfx_d_theta = bsxfun(@rdivide, d_P2_d_theta - d_P1_d_theta, dx);

se_mfx = sqrt(diag(d_mfx_d_theta*V_theta*d_mfx_d_theta'));
se_P = sqrt(diag(d_P_d_theta*V_theta*d_P_d_theta'));


mask = 0;
for i = 1:numel(spec.paramType)
    if spec.paramType(i) == 1 || spec.paramType(i) == 2
        mask(end+1) = mask(end) + 1;
    elseif spec.paramType(i) == 3 || spec.paramType(i) == 4
        mask(end+1:end+nchoice) = mask(end) + choiceset;
    end
end
mask(1) = [];
mask = bsxfun(@plus, mask, (0:nchoice-1)'*n.mfx)';
mask = mask(:);

mfx_full = zeros(n.mfx*n.maxChoice, ncon);
mfx_full(mask, :) = reshape(mfx, ncon, nmfx*nchoice)';
P_full = zeros(n.maxChoice, ncon);
P_full(choiceset,:) = reshape(P, ncon, nchoice)';
d_mfx_d_theta_full = zeros(n.mfx*n.maxChoice, k, ncon);
d_mfx_d_theta_full(mask,:,:) = permute(reshape(d_mfx_d_theta, [ncon, nmfx*nchoice, k]), [2,3,1]);
d_P_d_theta_full = zeros(n.maxChoice, k, ncon);
d_P_d_theta_full(choiceset,:,:) = permute(reshape(d_P_d_theta, [ncon, nchoice, k]), [2,3,1]);
end


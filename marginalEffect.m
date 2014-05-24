function [ mfx_full, P_full, se_mfx, se_P, dmfx_full, dP_full] = ...
    marginalEffect(thetaHat, meanData, n, spec, cov)
%MFX Summary of this function goes here
%   Detailed explanation goes here
%   paramType = 0 --> no marginal effect
%   paramType = 1 --> continuous and product independent
%   paramType = 2 --> binary and product independent
%   paramType = 3 --> continuous and product dependent
%   paramType = 4 --> binary and product dependent
% TODO
%   paramType = xy --> a set of dummies of type y 

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
    type = spec.paramType(i);
    sameGroup = spec.paramGroup == spec.paramGroup(i);
    sameID = spec.paramID == spec.paramID(i);
    
    if type == 0; continue; end
    
    if type <= 2 % fuel invariant
        if type == 1
            data1(index,sameID) = meanData(:,sameID)*(1-epsilon);
            data2(index,sameID) = meanData(:,sameID)*(1+epsilon);
        else
            data1(index, sameGroup) = 0;
            data2(index, sameGroup) = 0;
            data2(index, sameID) = 1;
        end
        
        dx = [dx;data2(index(1:nchoice:end),i) - data1(index(1:nchoice:end),i)];
        index = index + nchoice*ncon;
        continue;
    end
    
    for j=1:nchoice
        if type == 3
            data1(index(j:nchoice:end),sameID) = meanData(j:nchoice:end,sameID)*(1-epsilon);
            data2(index(j:nchoice:end),sameID) = meanData(j:nchoice:end,sameID)*(1+epsilon);
        else
            data1(index(j:nchoice:end), sameGroup) = 0;
            data2(index(j:nchoice:end), sameGroup) = 0;
            data2(index, sameID) = 1;
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
n.draw = 1000;
[dataR, ~] = ConstructDataGroup(data, n, spec);
dataR.draw.uni = repmat(dataR.draw.uni(1,:,:), [n.con 1 1]);
[P, dP] = ProbitProb(thetaHat, dataR);

m = (n.con - nchoice*ncon)/2;
P1 = P(1:m);
P2 = P(m+1:2*m);
dP1 = dP(:,1:m);
dP2 = dP(:,m+1:2*m);

mfx = (P2-P1)./dx;
dmfx = bsxfun(@rdivide, dP2 - dP1, dx');

P = P(2*m+1:end);
dP = dP(:,2*m+1:end);



se_mfx = sqrt(diag(dmfx'*cov*dmfx));
se_P = sqrt(diag(dP'*cov*dP));

k = numel(thetaHat);

mask = 0;
for i = 1:numel(spec.paramType)
    if spec.paramType(i) == 1 || spec.paramType(i) == 2
        mask(end+1) = mask(end) + 1;
    elseif spec.paramType(i) == 3 || spec.paramType(i) == 4
        mask(end+1:end+nchoice) = mask(end) + choiceset;
    end
end
mask(1) = [];
mask = bsxfun(@plus, mask, (choiceset-1)*n.mfx)';
mask = mask(:);

mfx_full = zeros(n.mfx*n.maxChoice, ncon);
mfx_full(mask, :) = reshape(mfx, ncon, nmfx*nchoice)';
P_full = zeros(n.maxChoice, ncon);
P_full(choiceset,:) = reshape(P, ncon, nchoice)';
dmfx_full = zeros(n.mfx*n.maxChoice, k, ncon);
dmfx_full(mask,:,:) = permute(reshape(dmfx', [ncon, nmfx*nchoice, k]), [2,3,1]);
dP_full = zeros(n.maxChoice, k, ncon);
dP_full(choiceset,:,:) = permute(reshape(dP', [ncon, nchoice, k]), [2,3,1]);
end


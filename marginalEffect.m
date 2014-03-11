function [ mfx, P, se_mfx, se_P] = marginalEffect(thetaHat, meanData, n, spec, V_theta)
%MFX Summary of this function goes here
%   Detailed explanation goes here
%   paramType = 0 --> no marginal effect
%   paramType = 1 --> continuous and product independent
%   paramType = 2 --> binary and product independent
%   paramType = 3 --> continuous and product dependent
%   paramType = 4 --> binary and product dependent

data1 = [];
data2 = [];
epsilon = 1e-3;
dx = [];

for i = 1:numel(spec.paramType)
    if spec.paramType(i) == 0; continue; end
    
    if spec.paramType(i) == 1 || spec.paramType(i) == 2
        x1 = meanData;
        x2 = meanData;
        
        if spec.paramType(i) == 1
            x1(:,i) = x1(:,i)*(1-epsilon);
            x2(:,i) = x2(:,i)*(1+epsilon);
        else
            x1(:,i) = 0;
            x2(:,i) = 1;
        end
        
        data1(end+1:end+n.maxChoice,:) = x1;
        data2(end+1:end+n.maxChoice,:) = x2;
        
        dx(end+1) = mean(x2(:,i) - x1(:,i));
        continue;
    end
    
    for j=1:n.maxChoice
        x1 = meanData;
        x2 = meanData;
        
        if spec.paramType(i) == 3
            x1(j,i) = x1(j,i)*(1-epsilon);
            x2(j,i) = x2(j,i)*(1+epsilon);
        else
            x1(j,i) = 0;
            x2(j,i) = 1;
        end
        
        data1(end+1:end+n.maxChoice,:) = x1;
        data2(end+1:end+n.maxChoice,:) = x2;
        
        dx(end+1) = x2(j,i) - x1(j,i);
    end
end

m = size(data1, 1);
data1 = repmat(data1, n.maxChoice, 1);
data2 = repmat(data2, n.maxChoice, 1);
data3 = repmat(meanData, n.maxChoice, 1);
dx = repmat(dx(:), n.maxChoice, 1);
for j = 1:n.maxChoice
    data1((j-1)*m+1:j*m,5) = j;
    data2((j-1)*m+1:j*m,5) = j;
    data3((j-1)*n.maxChoice+1:j*n.maxChoice,5) = j;
end

data = [data1; data2; data3];
conID = repmat(1:2*m+n.maxChoice,n.maxChoice,1);
data(:,2) = conID(:);

n.con = numel(unique(data(:,2)));
n.draw = 10000;
[dataR, ~] = ConstructDataGroup(data, n, spec);
dataR.draw.uni = repmat(dataR.draw.uni(1,:,:), [n.con 1 1]);
P = ProbitProb(thetaHat, dataR, n, spec);

m = (n.con - n.maxChoice)/2;
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
    P_plus(:,i) = ProbitProb(thetaHatPlus(:,i), dataR, n, spec);
    P_minus(:,i) = ProbitProb(thetaHatMinus(:,i), dataR, n, spec);
end

d_P_d_theta = bsxfun(@rdivide,P_plus - P_minus, d_theta');

d_P2_d_theta = d_P_d_theta(1:m,:);
d_P1_d_theta = d_P_d_theta(m+1:2*m,:);
d_P_d_theta = d_P_d_theta(2*m+1:end,:);
d_mfx_d_theta = bsxfun(@rdivide, d_P2_d_theta - d_P1_d_theta, dx);

se_mfx = d_mfx_d_theta*V_theta*d_mfx_d_theta';
se_P = d_P_d_theta*V_theta*d_P_d_theta';
end


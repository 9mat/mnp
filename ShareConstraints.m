function [ eq, neq, deq, dneq ] = ShareConstraints( theta, dataS, marketIdByCon, n, shareHat, mask)
%SHARECONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

probCon = zeros(n.con, n.maxChoice); 
d_probCon = zeros(numel(theta), n.con, n.maxChoice); 

for k = 1:size(dataS,1)
    for j = 1:size(dataS,2)
        if ~isempty(dataS{k,j})
             [prob, d_prob]= ProbitProb(theta(dataS{k,j}.pick), dataS{k,j});
             probCon(dataS{k,j}.allConID, j) = prob;
             d_probCon(dataS{k,j}.pick, dataS{k,j}.allConID, j) = d_prob;
        end
    end
end

allmarkets = sort(unique(marketIdByCon));
share = nan(numel(allmarkets), n.maxChoice);
d_share = nan(n.theta, n.market, n.maxChoice);
for k = 1:numel(allmarkets)
    share(k,:) = mean(probCon(marketIdByCon == allmarkets(k),:));
    d_share(:,k,:) = mean(d_probCon(:,marketIdByCon == allmarkets(k),:),2);   
end

share = share(mask.delta==1);
for i = 1:numel(theta)
    d = d_share(i,:,:);
    deq(i,:) = d(mask.delta==1);
end
eq = share(:) - shareHat(:);
neq = [];
dneq = [];

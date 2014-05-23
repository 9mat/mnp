function [ neq, eq, dneq, deq ] = ShareConstraints( theta, dataS, identifiable, n, shareHat)
%SHARECONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

eq = nan(size(identifiable.delta));
deq = zeros([numel(theta) size(eq)]);

shareHat = shareHat';
for m = 1:size(dataS,2)
    %fprintf('Market #%d\n', m);
    for j = 1:size(dataS,1);
        if ~isempty(dataS{j,m});
            dataS_0 = dataS{j,m};
            break;
        end;
    end;
    
    probCon = zeros(n.maxChoice, dataS_0.n.con);
    d_probCon = zeros(numel(theta), n.maxChoice, dataS_0.n.con);
    
    for j = 1:size(dataS,1)
        if ~isempty(dataS{j,m})
            [prob, d_prob]= ProbitProb(theta(dataS{j,m}.pick), dataS{j,m});
            probCon(j, :) = prob;
            d_probCon(dataS{j,m}.pick, j, :) = d_prob;
        end
    end
        
    eq(:,m) = mean(probCon,2) - shareHat(:,m);
    deq(:, :, m) = mean(d_probCon,3);
end

eq = eq(identifiable.delta);
deq = reshape(deq, n.theta, numel(identifiable.delta));
deq = deq(:, identifiable.delta(:));

neq = [];
dneq = [];

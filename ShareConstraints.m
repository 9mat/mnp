function [ neq, eq, dneq, deq ] = ShareConstraints( theta, dataS, identifiable, n, shareHat)
%SHARECONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

eq = nan(size(identifiable.delta));
deq = zeros([numel(theta) size(eq)]);

shareHat = shareHat';

n.marketcon = zeros(size(dataS,2),1);
for m = 1:size(dataS,2)
    for j = 1:size(dataS,1);
        if ~isempty(dataS{j,m});
            n.marketcon(m) = dataS{j,m}.n.con;
            break;
        end;
    end;
end

if nargout > 2
    for m = 1:size(dataS,2)
        probCon = zeros(n.maxChoice, n.marketcon(m));
        d_probCon = zeros(numel(theta), n.maxChoice, n.marketcon(m));
        
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
    
else
    for m = 1:size(dataS,2)
        probCon = zeros(n.maxChoice, n.marketcon(m));
        
        for j = 1:size(dataS,1)
            if ~isempty(dataS{j,m})
                probCon(j, :) = ProbitProb(theta(dataS{j,m}.pick), dataS{j,m});
            end
        end
        
        eq(:,m) = mean(probCon,2) - shareHat(:,m);
    end
    
    
    eq = eq(identifiable.delta);
end

neq = [];
dneq = [];

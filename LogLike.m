function [ nLogLike, d_nLogLike ] = LogLike( theta, dataR, spec )

if nargout == 1
    
    if exist('delta_0.mat', 'file') == 2
        load('delta_0.mat');
    end
    
    nLogLike = 0;
    for k = 1:numel(dataR)
        fprintf('    Solving for group %d\n', k);
        if exist('delta_0', 'var')
            delta0 = delta_0{k};
        else
            delta0 = zeros(dataR(k).n.maxChoice-1, dataR(k).n.market);
        end
        delta = solve_delta(theta(dataR(k).pick), dataR(k), spec, delta0);
        probChosen  = ProbitProb(theta(dataR(k).pick), delta, dataR(k), dataR(k).n, spec);
        probChosen( probChosen == 0 ) = eps;
        nLogLike = nLogLike - sum(log(probChosen));
        delta_000{k} = delta;
    end
    delta_0 = delta_000;
    save('delta_0.mat', 'delta_0');
    
elseif nargout > 1
    nLogLike = 0;
    d_nLogLike = zeros(size(theta));
    
    for k = 1:numel(dataR)
        [ probChosen, d_probChosen ]    = ProbitProb( theta(dataR(k).pick), dataR(k), dataR(k).n, spec );
        probChosen( abs(probChosen) < eps )   = eps;
        nLogLike    = nLogLike-sum( log( probChosen ) );
        d_nLogLike(dataR(k).pick)  = d_nLogLike(dataR(k).pick) ...
            - sum( bsxfun( @rdivide, d_probChosen, probChosen' ), 2 );
    end
    
    
    
end

    function sse = SSE(theta, delta, dataR, spec)
        prob = zeros(dataR.n.con, dataR.n.maxChoice);
        for j = 1:dataR.n.maxChoice
            if j ~= spec.base
                prob(:,j) = ProbitProb(theta,delta,dataR.dataS(j),dataR.n, spec);
            end
        end
        prob(:, spec.base) = 1 - sum(prob,2);
        
        share = zeros(dataR.n.market, dataR.n.maxChoice);
        for i = 1:dataR.n.con
            share(dataR.marketID(i),:) = share(dataR.marketID(i),:) + prob(i,:);
        end
        share = bsxfun(@times, share, 1./sum(share,2));
        
        sse = sum((share(:) - dataR.share(:)).^2);
    end

    function delta = solve_delta(theta, dataR, spec, delta0)
        obj = @(x) SSE(theta, x, dataR, spec );
        ub  = ones(size(delta0)) * 200;
        lb  = ones(size(delta0)) * -200;
        
        global opt;
        optOption   = optimset( 'MaxIter', 100, ...
            'MaxFunEvals', opt.maxFunEvals, ...
            'Display', 'off', ...
            'TolFun', opt.tolFun, ...
            'TolCon', opt.tolCon, ...
            'TolX', opt.tolX, ...
            'Algorithm', opt.algorithm, ...
            'GradObj', 'off', ...
            'GradConstr', 'off', ...
            'DerivativeCheck', 'on' );
        
        
        delta = ...
            ktrlink( obj, delta0, [], [], [], [], lb, ub, ...
            [], optOption );
    end
end

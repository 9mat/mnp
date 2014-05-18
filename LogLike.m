function [ nLogLike, d_nLogLike ] = LogLike( theta, dataR, spec )

if nargout == 1
    
    nLogLike = 0;
    for k = 1:numel(dataR)
        probChosen  = ProbitProb( theta(dataR{k}.pick), dataR{k}, dataR{k}.n, spec );
        probChosen( probChosen == 0 ) = eps;
        nLogLike = nLogLike - sum(log(probChosen));
    end
    
elseif nargout > 1
    nLogLike = 0;
    d_nLogLike = zeros(size(theta));
    
    for k = 1:numel(dataR)
        [ probChosen, d_probChosen ]    = ProbitProb( theta(dataR{k}.pick), dataR{k}, dataR{k}.n, spec );
        d_nLogLike(dataR{k}.pick)  = d_nLogLike(dataR{k}.pick) ...
            - sum( bsxfun( @rdivide, d_probChosen, probChosen' ), 2 );
        probChosen( probChosen < eps )   = eps;
        nLogLike    = nLogLike-sum( log( probChosen ) );
    end
    
    
    
end


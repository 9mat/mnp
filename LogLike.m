function [ nLogLike, d_nLogLike ] = LogLike( theta, dataR, spec )

n = dataR.n;

if nargout == 1
    
    nLogLike = 0;
    for k = 1:numel(dataR)
        probChosen  = ProbitProb( theta, dataR(k), n, spec );
        probChosen( probChosen == 0 ) = eps;
        nLogLike = nLogLike - sum(log(probChosen));
    end
    
elseif nargout > 1
    nLogLike = 0;
    d_nLogLike = zeros(size(theta));
    
    for k = 1:numel(dataR)
        [ probChosen, d_probChosen ]    = ProbitProb( theta, dataR(1), n, spec );
        probChosen( probChosen == 0 )   = eps;
        nLogLike    = nLogLike-sum( log( probChosen ) );
        d_nLogLike(dataR.pick)  = d_nLogLike(dataR.pick) ...
            - sum( bsxfun( @rdivide, d_probChosen, probChosen' ), 2 );
    end
    
end


function [ nLogLike, d_nLogLike ] = LogLike( theta, dataR )

if nargout == 1
    
    nLogLike = 0;
    for k = 1:numel(dataR)
        probChosen  = ProbitProb( theta(dataR{k}.pick), dataR{k} );
        probChosen( probChosen == 0 ) = eps;
        nLogLike = nLogLike - sum(log(probChosen));
    end
    
elseif nargout > 1
    nLogLike = 0;
    d_nLogLike = zeros(size(theta));
    
    for k = 1:numel(dataR)
        [ probChosen, d_probChosen ]    = ProbitProb( theta(dataR{k}.pick), dataR{k} );
        probChosen( abs(probChosen) < eps )   = eps;
        nLogLike    = nLogLike-sum( log( probChosen ) );
        d_nLogLike(dataR{k}.pick)  = d_nLogLike(dataR{k}.pick) ...
            - sum( bsxfun( @rdivide, d_probChosen, probChosen' ), 2 );
    end
    
end


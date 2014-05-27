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
%         probChosen(probChosen == 0) = eps;
        d_ln_prob = bsxfun( @rdivide, d_probChosen, probChosen' );
        d_ln_prob(isnan(d_ln_prob)) = 0;
        d_nLogLike(dataR{k}.pick)  = d_nLogLike(dataR{k}.pick) ...
            - sum( d_ln_prob, 2 );
        probChosen( abs(probChosen) < eps )   = eps;
        nLogLike    = nLogLike - sum( log( probChosen ) );
    end
    
    d_nLogLike(d_nLogLike == -Inf) = -1e10;
    d_nLogLike(d_nLogLike == Inf) = 1e10;    
    
end


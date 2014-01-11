function [ nLogLike, d_nLogLike ] = LogLike( theta, dataR, spec )

n = dataR.n;

if nargout == 1
    
    probChosen  = ProbitProb( theta, dataR(1), n, spec );
    
    % Modify 0 values in 'probChosen' to machine epsilon to avoid log( 0 )
    %   or division by 0  
    probChosen( probChosen == 0 )   = eps;

    nLogLike    = -sum( log( probChosen ) );

elseif nargout > 1
   
    [ probChosen, d_probChosen ]    = ProbitProb( theta, dataR(1), n, spec );
    
    % Modify 0 values in 'probChosen' to machine epsilon to avoid log( 0 )
    %   or division by 0        
    probChosen( probChosen == 0 )   = eps;
        
    nLogLike    = -sum( log( probChosen ) );
    d_nLogLike  = -sum( bsxfun( @rdivide, d_probChosen, probChosen' ), 2 ); 
    
end


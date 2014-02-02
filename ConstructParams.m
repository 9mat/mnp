function [ params ] = ConstructParams( theta, n, spec )

%% Construct Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size( theta, 2 ) > size( theta, 1 )
    theta = theta';
end

% Common price parameter ( alpha_0 )
%   alpha_0 is [ 1 x 1 ]
temp    = 0;
if ( 1 - spec.unobs ) > 0
    params.alpha_0  = theta( temp + 1 );    
end

% Consumer group specific price parameters ( alpha_r )
%   alpha_r is [ 1 x n.conGroup ]
temp    = 1 * ( 1 - spec.unobs );
if n.conGroup > 0
    params.alpha_r 	= theta( temp + 1 : temp + n.conGroup )'; 
end

% Product characteristic parameters ( beta_1 )
%   beta_1 is [ n.maxChoice x n.prodChar ]
temp    = temp + n.conGroup;
if ( 1 - spec.unobs ) * n.prodChar > 0
    params.beta_1 	= theta( temp + 1 : temp + n.prodChar * n.maxChoice );
    params.beta_1 	= reshape( params.beta_1, [ n.maxChoice n.prodChar ] );  
end

% Consumer characteristic parameters ( beta_2 )
%   beta_2 is [ n.maxChoice x n.conChar ]
temp    = temp + n.prodChar * n.maxChoice * ( 1 - spec.unobs );
if n.conChar > 0
    temp_beta_2     = theta( temp + 1 : ...
                             temp + n.conChar * ( n.maxChoice - 1 ) );
    temp_beta_2 	= reshape( temp_beta_2, ...
                               [ ( n.maxChoice - 1 ) n.conChar ] );
    params.beta_2   = [ temp_beta_2( 1 : spec.base - 1, : );
                        zeros( 1, n.conChar );
                        temp_beta_2( spec.base : end, : ) ];
end

% FE
% delta is [ (n.maxChoice - 1) x n.market ]
temp        = temp + n.conChar * ( n.maxChoice - 1 );
params.delta  = zeros(n.maxChoice - 1, n.market);
params.delta(:) = theta(temp + 1 : temp + numel(params.delta));

% Choleski factor of the ( differenced ) covariance matrix 
%   - the base alternative is spec.base 
%   - S is a lower-triangular matrix
%   - S is [ n.maxChoice x n.maxChoice ]
temp        = temp + numel(params.delta);

mask.S = tril(true(n.maxChoice - 1));
if numel(theta) - temp < sum(mask.S(:))
    mask.S(1,1) = false;
end
params.S = eye( n.maxChoice - 1, n.maxChoice - 1 );
params.S(mask.S) = theta(temp+1:end);

clear temp

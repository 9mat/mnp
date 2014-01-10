function [ cIneq, cEq, d_cIneq, d_cEq ] = NonLCon( theta, n, bound )

%% Compute Non-Linear Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inequality constraints 
cIneq       = zeros( n.maxChoice - 2, 1 );
k           = n.beta;

for j = 1 : n.maxChoice - 2    
    k           = k + ( n.maxChoice - j );        
    cIneq(j)    = bound - theta(k) .^ 2;
end
clear k

% Equality constraints 
cEq     = [];

%% Compute Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 2
    
    % Inequality constraints 
    d_cIneq     = zeros( n.theta, n.maxChoice - 2 );    
    k           = n.beta;
    
    for j = 1 : n.maxChoice - 2
        k               = k + ( n.maxChoice - j );
        d_cIneq( k, j ) = -2 * theta(k);
    end
        
    % Equality constraints 
    d_cEq       = [];
end


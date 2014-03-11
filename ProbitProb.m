function [ probChosen, d_probChosen ] = ProbitProb( theta, dataR, n, spec )

%% Construct Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta = [ theta; 2; 3; 4; 5; 6;];
% theta = [ theta(1:end-2); theta(end-1); theta(end); 7; 4; 1;];
params 	= ConstructParams( theta, dataR.n, dataR.spec );
base    = dataR.spec.base;

% Common price parameter ( alpha_0 )
%   alpha_0 is [ 1 x 1 ]
if ( 1 - spec.unobs ) > 0
    alpha_0     = params.alpha_0;  
end

% Consumer group specific price parameters ( alpha_r )
%   alpha_r is [ 1 x n.conGroup ]
if n.conGroup > 0
    alpha_r     = params.alpha_r;   
end

% Product characteristic parameters ( beta_1 )
%   beta_1 is [ n.maxChoice x n.prodChar ]
if ( 1 - spec.unobs ) * n.prodChar > 0
    beta_1      = params.beta_1;
end

% Consumer characteristic parameters ( beta_2 )
%   beta_2 is [ n.maxChoice x n.conChar ]
if n.conChar > 0
    beta_2      = params.beta_2;  
end

% Choleski factor of the ( differenced ) covariance matrix 
%   - the base alternative is dataR.base 
%   - S is a lower-triangular matrix
%   - S is [ n.maxChoice x n.maxChoice ]
S	= params.S;

%% Compute ( Differenced ) Deterministic Utilities %%%%%%%%%%%%%%%%%%%%%%%%

params.delta = reshape(params.delta, n.maxChoice - 1, n.market);
V   = params.delta(:, dataR.marketID);
V   = [V(1:spec.base-1, :); zeros(1, n.con); V(spec.base:end, :)];
ind = sub2ind(size(V), dataR.choice, 1:n.con);
V   = bsxfun(@minus, V, V(ind));
V(ind) = [];
V = reshape(V, n.maxChoice-1, n.con);

if ( 1 - spec.unobs ) > 0
    %   alpha_0 is [ 1 x 1 ]
    %   dataR.diff.price  is [ n.maxChoice - 1 x 1 x n.con ]
    V   = V + alpha_0 * permute(dataR.diff.price, [1 3 2]);
end

if n.conGroup > 0 
    %   alpha_r is [ 1 x n.conGroup ]
    %   dataR.diff.conGroupP is [ 1 x n.conGroup x n.con ]
    temp    = bsxfun( @times, alpha_r, dataR.diff.conGroupP );   

    V       = V + permute(sum( temp, 2 ), [1 3 2] );
    clear temp
end
       
if ( 1 - spec.unobs ) * n.prodChar > 0            
    %   beta_1 is [ n.maxChoice x n.prodChar ]
    %   dataR.diff.prodChar is [ (n.maxChoice - 1) x (n.maxChoice * n.prodChar) x n.con ]
    beta_1  = reshape( beta_1, [ 1 (n.maxChoice * n.prodChar) 1 ] );
    temp    = bsxfun( @times, beta_1, dataR.diff.prodChar );

    V       = V + permute( sum( temp, 2 ), [1 3 2] );
    clear temp
end

if n.conChar > 0
    %   beta_2 is [ n.maxChoice x n.conChar ]
    %   dataR.diff.conChar is [ (n.maxChoice - 1) x (n.maxChoice * n.conChar) x n.con ]
    beta_2  = reshape( beta_2, [ 1 (n.maxChoice * n.conChar) 1 ] );
    temp    = bsxfun( @times, beta_2, dataR.diff.conChar );

    V       = V + permute( sum( temp, 2 ), [1 3 2] );
    clear temp
end


%% Difference the Choleski Factor of the Covariance Matrix %%%%%%%%%%%%%%%%
 
S_j         = zeros( n.maxChoice - 1, n.maxChoice - 1, n.maxChoice );                

for j = 1 : n.maxChoice
    MS = dataR.M(:,:,j,base) * S;
    S_j(:,:,j) = mychol(MS*MS'); % use a more stable decomposition
end
     
S_i = S_j(:,:,dataR.choice);

%% Calculate Probabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
a           = repmat(permute(V, [2 3 1]), [1 n.draw 1]);    
w           = zeros( n.con, n.draw, n.maxChoice - 2 );
ub          = zeros( n.con, n.draw, n.maxChoice - 1 ); 

% In the loops, use temporary matrices aj and ubj to avoid repeated 
% indexing to matrices a and ub, which can be very costly
for j = 1 : n.maxChoice - 1        
    aj = repmat(permute(V(j, :), [2 3 1]), [1 n.draw]);
    if j > 1
        w(:,:,j-1) = norminv(ubj .* dataR.draw.uni(:,:,j-1));
        for h = 1:j-1
            aj = aj + bsxfun( @times, w(:,:,h), squeeze(S_i(j,h,:)));
        end
    end
    aj = bsxfun(@rdivide, aj, -squeeze(S_i(j,j,:)));
    ubj = 0.5*erfc(-aj/sqrt(2)); % 3x faster than normcdf
    
    a(:,:,j) = aj;
    ub(:,:,j) = ubj;
end
clear aj ubj;

ub(ub==0) = eps;
probChosen = prod(ub,3);

if nargout == 1
    probChosen 	= mean( probChosen, 2 );  
end

%% Compute Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mn = permute(dataR.M(:, :, dataR.choice, dataR.spec.base), [2 3 1]);

if nargout > 1     
            
    % Pre-compute normpdf(a) and normpdf(w) to speed up computation
    normpdf_a                   = normpdf(a);      
    normpdf_w                   = normpdf(w);
    
    % Modify 0 values in 'ub' to machine epsilon to avoid division 
    %   by 0  
    normpdf_w( normpdf_w == 0 ) = eps;
    
    u_normpdf_a_w   = dataR.draw.uni .* normpdf_a(:,:,1:end-1) ./ normpdf_w; 
    normpdf_a_ub    = reshape( normpdf_a ./ ub, ...
                               [ 1 n.con n.draw ( n.maxChoice - 1 ) ] );
    clear normpdf_w ub;      
    
    % Derivatives of V wrt to beta
    if ( 1 - spec.unobs ) > 0
        d_V_beta    = permute( dataR.diff.price, [ 2 3 1 ] ); 
    else
        d_V_beta    = [];
    end
    
    if n.conGroup > 0
        d_V_beta    = [ d_V_beta; ...
                        permute( dataR.diff.conGroupP, [ 2 3 1 ] ) ];
    end
    
    if ( 1 - spec.unobs ) * n.prodChar > 0        
        d_V_beta    = [ d_V_beta; ....
                        permute( dataR.diff.prodChar, [ 2 3 1 ] ) ];        
    end
    
    if n.conChar > 0        
        d_V_beta    = [ d_V_beta; ....
                        permute( dataR.diff.conChar, [ 2 3 1 ] ) ];        
    end
    
    d_V_beta    = reshape( d_V_beta, ...
                           [ n.beta_all n.con 1 ( n.maxChoice - 1 ) ] );   
    
    % Derivatives of a wrt to beta
    d_a_beta    = zeros( n.beta_all, n.con, n.draw, n.maxChoice - 1 );    
    
    d_a_delta   = zeros( n.maxChoice-1, n.con, n.draw, n.maxChoice - 1 );    
    
    % Derivatives of w wrt to beta
    d_w_beta    = reshape( u_normpdf_a_w, ...
                           [ 1 n.con n.draw ( n.maxChoice - 2 ) ] );             
    d_w_beta    = repmat( d_w_beta, [ n.beta_all 1 1 1 ] ); 
 
    d_w_delta   = reshape( u_normpdf_a_w, ...
                           [ 1 n.con n.draw ( n.maxChoice - 2 ) ] );             
    d_w_delta   = repmat( d_w_delta, [ n.maxChoice-1 1 1 1 ] ); 

    % Derivtatives of a wrt to s_i
    d_a_s_i     = zeros( n.maxChoice - 1, n.maxChoice - 1, n.con, ...
                         n.draw, n.maxChoice - 1 );
                     
    % Derivatives of w wrt to s_i              
    d_w_s_i     = reshape( u_normpdf_a_w( :, :, : ) , ...
                           [ 1 1 n.con n.draw ( n.maxChoice - 2 ) ] );
    
    for l = 1 : n.maxChoice - 1 
        
        % Derivatives of a wrt beta
        if l==1
            d_a_beta_l = repmat( d_V_beta(:,:,:,l), [1 1 n.draw]);
        elseif l > 1            
                    
            d_w_beta( :, :, :, l - 1 )  = d_w_beta( :, :, :, l - 1 ) .* d_a_beta_l;
            d_a_beta_l = repmat( d_V_beta(:,:,:,l), [1 1 n.draw]);
            for h = 1 : l - 1            
                d_a_beta_l  = d_a_beta_l + bsxfun( @times, d_w_beta( :, :, :, h ), ...
                                            squeeze( S_i( l, h, : ) )' );
            end     
        end
        d_a_beta_l     = bsxfun( @rdivide, d_a_beta_l, -squeeze(S_i( l, l, : ) )' );
        d_a_beta( :, :, :, l )      = d_a_beta_l;
            
        % Derivatives of a wrt s_i                                      
        d_a_s_i( l, l, :, :, l ) 	= -bsxfun( @rdivide, ...
                    reshape( a( :, :, l ), [ 1 1 n.con n.draw 1 ] ), ...
                    reshape( S_i( l, l, : ), [ 1 1 n.con 1 1 ] ) );
                
        for i = 1 : n.maxChoice - 1
            for j = 1 : n.maxChoice - 1
                if ( i < l ) && ( i >= j )                    
                    for h = 1 : ( l - 1 )
                        temp    = d_a_s_i( i, j, :, :, h );
                        temp    = bsxfun( @times, temp, ...
                                          d_w_s_i( :, :, :, :, h ) );
                        temp    = bsxfun( @times, temp, ...
                                          reshape( S_i( l, h, : ) ./ ...
                                                   S_i( l, l, : ), ...
                                                   [ 1 1 n.con 1 1 ] ) );
                        
                        d_a_s_i( i, j, :, :, l )    = -temp + ...
                                                d_a_s_i( i, j, :, :, l );
                    end
                    clear temp;
                    
                elseif ( i == l ) && ( i > j )
                    d_a_s_i( i, j, :, :, l )    = -bsxfun( @rdivide, ...
                        reshape( w( :, :, j ), [ 1 1 n.con n.draw 1 ] ), ...
                        reshape( S_i( l, l, : ), [ 1 1 n.con 1 1 ] ) );
                end                
            end
        end
        
        % derivative of a wrt to delta (i.e. fixed effects)      
        if l==1
            d_a_delta_l = repmat( Mn(:,:,l), [1 1 n.draw]);
        elseif l > 1            
            d_w_delta( :, :, :, l - 1 )  = d_w_delta( :, :, :, l - 1 ) .* d_a_delta_l;
            d_a_delta_l = repmat( Mn(:,:,l), [1 1 n.draw]);
            for h = 1 : l - 1
                d_a_delta_l  = d_a_delta_l + bsxfun( @times, d_w_delta( :, :, :, h ), ...
                    squeeze( S_i( l, h, : ) )' );
            end
        end
        d_a_delta_l     = bsxfun( @rdivide, d_a_delta_l, -squeeze(S_i( l, l, : ) )' );
        d_a_delta( :, :, :, l )      = d_a_delta_l;
        
    end 
    clear d_w_beta d_w_s_i d_V_beta S_i d_a a w;  
    
    % Derivatives of a wrt s
    d_a_s       = zeros( n.s + 1, n.con, n.draw, n.maxChoice - 1 );
    
    temp1        = 0;
    for j = 1 : n.maxChoice - 1
        temp2   = n.maxChoice - j;
        
        d_a_s( temp1 + 1 : temp1 + temp2, :, :, : ) = ...
                reshape( d_a_s_i( j : n.maxChoice - 1, j, :, :, : ), ...
                         [ temp2 1 n.con n.draw ( n.maxChoice - 1 ) ] );
                 
        temp1    = temp1 + temp2;         
    end
    clear temp1 temp2 d_a_s_i;
    
    % Derivatives of L wrt beta
    d_L_beta    = bsxfun( @times, normpdf_a_ub, d_a_beta );
    d_L_beta    = sum(d_L_beta , 4 );
    d_L_beta    = ...
        mean( bsxfun( @times, d_L_beta, ...
                      reshape( probChosen, [ 1 n.con n.draw ] ) ), 3 );  
    
    clear d_a_beta;
    
    % Derivatives of L wrt delta
    d_L_delta   = sum( bsxfun( @times, normpdf_a_ub, d_a_delta ), 4 );
    d_L_delta   = ...
        mean( bsxfun( @times, d_L_delta, ...
                      reshape( probChosen, [ 1 n.con n.draw ] ) ), 3 );  
    
    clear d_a_delta;
    
    d_L_delta_full = zeros((n.maxChoice - 1)*n.market, n.con);
    for k = 1:n.market
        marketindex = dataR.marketID == k;
        feindex = (n.maxChoice - 1)*(k-1)+1:(n.maxChoice - 1)*k;
        d_L_delta_full(feindex, marketindex) = d_L_delta(:,marketindex);
    end

                  
    % Derivatives of L wrt s_i
    d_L_s_i     = sum( bsxfun( @times, normpdf_a_ub, d_a_s ), 4 );
    d_L_s_i     = ...
        mean( bsxfun( @times, d_L_s_i, ...
                      reshape( probChosen, [ 1 n.con n.draw ] ) ), 3 ); 
    d_L_s_i( d_L_s_i == Inf )   = 1e+10;
    
    clear normpdf_a_ub d_a_s;

    % Derivatives of s_n wrt s
    
    % Because matrices A and B depends on choices but not on individuals,
    % first calculate A and B for each choice and replicate them to
    % all the consumers according to their choices 
    % (note: in Bolduc, equation (B2) to (B4), A and B depend on consumers
    % because of individual choice sets, which here assumed away)
    for j = 1 : n.maxChoice
       A_j  = pinv( ...
                        kron( S_j(:,:,j), eye( n.maxChoice - 1 ) ) * dataR.L + ...
                        kron( eye( n.maxChoice - 1 ), S_j(:,:,j) ) * dataR.K )';
       B_j  = ( kron( dataR.M( :, :, j, base ) * S, ...
                              dataR.M( :, :, j, base ) ) * ...
                        dataR.L + ...
                        kron( dataR.M( :, :, j, base ), ...
                              dataR.M( :, :, j, base ) * S ) * ...
                        dataR.K )';
       BA_j(:,:,j)  = B_j * A_j;
    end
    
    % Derivatives of s_n wrt s
    d_s_n_s = BA_j(:,:,dataR.choice);
    
    clear A_j B_j BA_j S_j;
    
    % Derivatives of L wrt s
    d_L_s       = zeros( n.s_all, n.con );
    for i = 1 : n.con        
        d_L_s( :, i )        = d_s_n_s( :, :, i ) * d_L_s_i( :, i );
    end    
    
    clear d_s_n_s d_L_s_i;
        
    probChosen      = mean( probChosen, 2 );
    d_probChosen    = [ d_L_beta( dataR.betaIndex, : ); 
                        d_L_delta_full;
                        d_L_s( dataR.sIndex, : ) ];
    
end















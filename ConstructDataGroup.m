
function [dataR, data] = ConstructDataGroup(dataMatrix, n, spec)

% check that data is relatively balanced (i.e. choice sets are the same for
% all consumers)
alternative = dataMatrix(:, 4);
allchoices = unique(alternative);
count = histc(alternative, allchoices);
assert(all(count(2:end)==count(1)));

% check if the choice set includes the scale alternative, otherwise we
% need to estimate the first element of the covariance matrix
dataR.missingscale = ~any(allchoices == spec.scale);

% sort data by alternative
% !!!important for later reshaping
[~, sortindex] = sort(alternative);
dataMatrix = dataMatrix(sortindex,:);

data.marketID       = dataMatrix( :, 1 );
data.conID          = dataMatrix( :, 2 );
data.choiceSetSize  = dataMatrix( :, 3 );
data.alternative    = dataMatrix( :, 4 );
data.choice         = dataMatrix( :, 5 );
data.price          = dataMatrix( :, 6 );

% recode choices and alternatives in running order (simplify later code??)
dataR.mapChoice = sort(allchoices);
for i = 1:numel(allchoices)
    % use negative sign to differentiate old and new codes
    data.alternative(data.alternative == allchoices(i)) = -i;
    data.choice(data.choice == allchoices(i)) = -i;
    
    if allchoices(i) == spec.base
        spec.base = i;
    end
end

data.alternative = - data.alternative;
data.choice = - data.choice;

dataR.allConID = sort(unique(data.conID));

%recode market ID
dataR.allmarkets = sort(unique(data.marketID));
for i = 1:numel(dataR.allmarkets)
    data.marketID(data.marketID == dataR.allmarkets(i)) = -i;
end
data.marketID = -data.marketID;

% Matrix of consumer group indicators ( i in r )
temp            = 6;
if n.conGroup > 0
    data.conGroup   = dataMatrix( :, temp + 1 : temp + n.conGroup );
end

% Matrix of product characteric variables ( x_jl )
temp            = 6 + n.conGroup;
if n.prodChar > 0
    data.prodChar   = dataMatrix( :, temp + 1 : temp + n.prodChar );
end

% Matrix of consumer characteric variables ( x_i )
temp            = 6 + n.conGroup + n.prodChar;
if n.conChar > 0
    data.conChar    = dataMatrix( :, temp + 1 : temp + n.conChar );
end
clear temp

% Number of markets
n.market    = length( unique( data.marketID ) );

% Number of consumers
n.con       = length( unique( data.conID ) );

% Largest number of choices
n.maxChoice = max( data.choiceSetSize );

%% Number of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% beta = [ alpha_0; alpha_r; beta_1; beta_2 ]
n.beta      = 1 * ( 1 - spec.unobs ) + ...
              n.conGroup + ...
              n.prodChar * n.maxChoice * ( 1 - spec.unobs ) + ...
              n.conChar * ( n.maxChoice - 1 );
n.beta_all  = n.beta + n.conChar;

% s = elements of S, where Omega = S * S'
n.s         = ( n.maxChoice - 1 ) * n.maxChoice / 2 - 1;
n.s_all     = n.s + 1;

n.delta     = (n.maxChoice - 1)*n.market;

% theta = [ beta; s ];
n.theta     = n.beta + n.s;          
n.theta_all = n.beta_all + n.s_all;

% The following indexing variables indicate which parameters are NOT
%   normalized and are to be estimated
if n.conChar > 0
    dataR.betaIndex     = true( n.beta_all, 1 );
    temp                = n.beta - n.conChar * ( n.maxChoice - 1 );
    for k = 1 : n.conChar
        dataR.betaIndex( temp + spec.base )  = 0;
           
        temp     = temp + ( n.maxChoice );  
    end
else
    dataR.betaIndex     = true( n.beta, 1 );
end
clear temp
    
dataR.sIndex    = [ false; true( n.s, 1 ) ];
if dataR.missingscale
    dataR.sIndex(1) = true;
end

%% Reshape Data Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a data matrix "dataR" by reshaping "data"
dataR.conID         = zeros( n.maxChoice, 1 );
dataR.choiceSetSize = zeros( n.maxChoice, 1 );
dataR.alternative   = zeros( n.maxChoice, 1 );

dataR.choice        = zeros( 1, n.con );
dataR.marketID      = zeros( 1, n.con );
dataR.price         = zeros( n.maxChoice, 1, n.con );
dataR.conGroup      = zeros( n.maxChoice, n.conGroup, n.con );
dataR.conGroupP     = zeros( n.maxChoice, n.conGroup, n.con );
dataR.prodChar     	= zeros( n.maxChoice, n.prodChar, n.con );
dataR.conChar       = zeros( n.maxChoice, n.conChar, n.con );

for i = 1 : n.con
    
    % NOTE: this reshaping depends crucial on the sorting in the beginning,
    % i.e. data must be sorted by alternative first
    index1          = ( data.conID == data.conID(i) );
    
    dataR.choice(i) = unique( data.choice( index1 ) );
    dataR.marketID(i) = unique( data.marketID( index1 ) );
        
    dataR.price( :, :, i )          = data.price( index1, : );
    dataR.alternative( :, :, i )    = data.alternative( index1, : );
    
    if n.conGroup > 0
        dataR.conGroup( :, :, i )   = data.conGroup( index1, : );
        dataR.conGroupP( :, :, i )  = bsxfun( @times, ...
            data.conGroup( index1, : ), ...
            data.price( index1, : ) );
    end
    
    if n.prodChar > 0
        dataR.prodChar( :, :, i )   = data.prodChar( index1, : );
    end
    
    if n.conChar > 0
        dataR.conChar( :, :, i )    = data.conChar( index1, : );
    end
end

%% Construct Differencing Matrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dataR.N( :, :, j ) is a linear operator that differences against 
%   alternative j
dataR.N     = zeros( n.maxChoice - 1, n.maxChoice, n.maxChoice );

temp        = eye( n.maxChoice - 1 );
for j = 1 : n.maxChoice
    dataR.N( :, :, j )  = [ temp( :, 1 : j - 1 ), ...
                            -ones( n.maxChoice - 1, 1 ), ...
                            temp( : ,j : end ) ];    
end
clear temp

% dataR.M( :, :, j, k ) is a linear operator that changes the
%    base alternative from k to j
dataR.M     = zeros( n.maxChoice - 1, n.maxChoice - 1, ...
                     n.maxChoice , n.maxChoice );

temp1   = eye( n.maxChoice - 2 );
for j = 1 : n.maxChoice
    for k = 1 : n.maxChoice
        
        if j == k
            dataR.M( :, :, j, k )   = eye( n.maxChoice - 1 );
        else
            c                       = ( j < k ) * j + ...
                                      ( j > k ) * ( j - 1 );
            r                       = ( j < k ) * ( k - 1 ) + ...
                                      ( j > k ) * k;
                          
            temp2                   = [ temp1( :, 1 : c - 1 ), ...
                                        -ones( n.maxChoice - 2, 1) ...
                                        temp1( :, c : end ) ];
            temp3                   = zeros( 1, n.maxChoice - 1 );
            temp3(c)                = -1;            
            
            dataR.M( :, :, j, k )   = [ temp2( 1 : r - 1, : ); ...
                                        temp3; ...
                                        temp2( r : end , : ) ]; 
            clear c r temp2 temp3
        end
    end
end
clear temp1

%% The Next Two Matrices Are Useful For Computing Derivatives %%%%%%%%%%%%%

% dataR.L is a linear operator such that S(:) = L * [ 1; s ]
dataR.L = zeros( ( n.maxChoice - 1 )^ 2, n.s + 1 );
temp1   = tril( ones( n.maxChoice - 1, n.maxChoice - 1 ), 0 );
temp2   = find( temp1 == 1 );

for j = 1 : n.s + 1
    dataR.L( temp2(j), j ) = 1;
end
clear temp1 temp2

% dataR.K is a linear operator such that SS(:) = K * [ 1; s ], where 
%   SS = S'
dataR.K = zeros( ( n.maxChoice - 1 )^ 2, n.s + 1 );
temp1   = linspace( 1, 1 + n.s, 1 + n.s )';
temp2   = dataR.L * temp1;
temp2   = reshape( temp2, n.maxChoice - 1, n.maxChoice - 1 )';
for j = 1 : n.s + 1
    find( temp2 == temp1(j) );
    dataR.K( temp2 == temp1(j), j ) = 1;
end
clear temp1 temp2

%% Difference the Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( 1 - spec.unobs ) > 0
    dataR.diff.price        = zeros( n.maxChoice - 1, 1, n.con );         
end  

if n.conGroup > 0  
    dataR.diff.conGroupP    = zeros( n.maxChoice - 1, n.conGroup, n.con );  
end

if ( 1 - spec.unobs ) * n.prodChar > 0
   dataR.diff.prodChar      = ...
                zeros( n.maxChoice - 1, n.prodChar * n.maxChoice, n.con );
end

if n.conChar > 0
   dataR.diff.conChar       = ...
                zeros( n.maxChoice - 1, n.conChar * n.maxChoice, n.con );
end
            
for i = 1 : n.con        
    if ( 1 - spec.unobs ) > 0
        dataR.diff.price( :, :, i )     = ...
            dataR.N( :, :, dataR.choice(i) ) * dataR.price( :, :, i );            
    end  

    if n.conGroup > 0  
        dataR.diff.conGroupP( :, :, i ) = ...
            dataR.N( :, :, dataR.choice(i) ) * dataR.conGroupP( :, :, i );
    end
end

if ( 1 - spec.unobs ) * n.prodChar > 0 
    temp    = zeros( n.prodChar * n.maxChoice, n.maxChoice - 1, n.con );

    for i = 1 : n.con
        for j = 1 : n.maxChoice
            if j == dataR.choice(i)
                for l = 1 : n.prodChar
                       temp( ( l -  1 ) * n.maxChoice + j, : ,i ) = ...
                           -repmat( dataR.prodChar( j, l, i ), ...
                                    [ 1 ( n.maxChoice - 1 ) 1 ] );   
                end
            else
                r   = ( j < dataR.choice(i) ) * j + ...
                      ( j > dataR.choice(i) ) * ( j - 1 );

                for l = 1 : n.prodChar
                    temp( ( l - 1 ) * n.maxChoice + j, r ,i ) =  ...
                                            dataR.prodChar( j, l, i );   
                end
            end
        end
    end

    dataR.diff.prodChar = permute( temp, [ 2 1 3 ] ); 
    clear temp
end

if n.conChar > 0
    temp    = zeros( n.conChar * n.maxChoice, n.maxChoice - 1, n.con );

    for i = 1 : n.con
        for j = 1 : n.maxChoice
            if j == dataR.choice(i)
                for l = 1 : n.conChar
                       temp( ( l -  1 ) * n.maxChoice + j, : ,i ) = ...
                           -repmat( dataR.conChar( j, l, i ), ...
                                    [ 1 ( n.maxChoice - 1 ) 1 ] );   
                end
            else
                r   = ( j < dataR.choice(i) ) * j + ...
                      ( j > dataR.choice(i) ) * ( j - 1 );

                for l = 1 : n.conChar
                    temp( ( l - 1 ) * n.maxChoice + j, r ,i ) =  ...
                                            dataR.conChar( j, l, i );   
                end
            end
        end
    end

    dataR.diff.conChar = permute( temp, [ 2 1 3 ] ); 
    clear temp
end

%% Make Random Draws %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if spec.drawType == 1
    
    % Make uniform draws
    dataR.draw.uni  = rand( n.maxChoice - 2, n.con, n.draw );
    
elseif spec.drawType == 2 && n.maxChoice > 2
    
    tempDraw        = haltonset( ( n.maxChoice - 2 ) * n.con, ...
                                 'Skip', spec.halton.skip, ...
                                 'Leap', spec.halton.leap );
    tempDraw        = scramble( tempDraw, spec.halton.scramble );
    
    % Make uniform draws
    dataR.draw.uni  = net( tempDraw, n.draw );
    dataR.draw.uni  = reshape( dataR.draw.uni, ...
                        [ ( n.maxChoice - 2 ) n.con n.draw ] );
    dataR.draw.uni  = permute(dataR.draw.uni, [2 3 1]);  
    
elseif n.maxChoice == 2
    dataR.draw.uni = zeros(n.con, n.draw, 0);
end

dataR.n = n;
dataR.spec = spec;


        

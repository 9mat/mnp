
%% Construct Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% thetaHat = [ thetaHath; 2; 3; 4; 5; 6;];
% thetaHat = [ thetaHat(1:end-1); 0; thetaHat(end); 0; 0; 1;];
params  = ConstructParams( thetaHat, n, spec );
params.se = ConstructParams( MLE.se, n, spec );

% Common price parameter ( alpha_0 )
%   alpha_0 is [ 1 x 1 ]
if ( 1 - spec.unobs ) > 0
    alpha_0 = params.alpha_0;    
end

% Consumer group specific price parameters ( alpha_r )
%   alpha_r is [ 1 x n.conGroup ]
if n.conGroup > 0
    alpha_r = params.alpha_r;   
end

% Product characteristic parameters ( beta_1 )
%   beta_1 is [ n.maxChoice x n.prodChar ]
if ( 1 - spec.unobs ) * n.prodChar > 0
    beta_1  = params.beta_1;
end

% Consumer characteristic parameters ( beta_2 )
%   beta_2 is [ n.maxChoice x n.conChar ]
if n.conChar > 0
    beta_2 	= params.beta_2;  
end

% Choleski factor of the (differenced) covariance matrix 
%   - the base alternative is spec.base 
%   - L is a lower-triangular matrix
%   - L is [ n.maxChoice x n.maxChoice ]
S	= params.S;
                                
% Covariance matrix
omega                       = zeros( n.maxChoice, n.maxChoice );
omega( 2 : end, 2 : end )   = S * S';

delta = NaN(n.maxChoice, n.market);
delta(mask.delta==1) = params.delta;
delta(spec.base,:) = 0;

se_delta = NaN(size(delta));
se_delta(mask.delta==1) = params.se.delta;
se_delta(spec.base,:) = 0;
clear temp

%% Print Some Prelimiary Informations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf( '\n\n ******************************************************\n' );
fprintf( ' ******************************************************\n' );
fprintf( ' \n\t Main Data File:\t%s\n', spec.dataName );
fprintf( ' \n\t Number of consumers:\t\t\t\t\t%5i\n', n.con );
fprintf( ' \n\t Number of alternatives:\t\t\t\t%2i\n', n.maxChoice );
fprintf( ' \n\t Base alternative:\t\t\t\t\t\t%2i\n', spec.base );
fprintf( ' \n\t Number of consumer groups:\t\t\t\t%2i\n', n.conGroup );
fprintf( ' \n\t Number of product characteristics:\t\t%2i\n', n.prodChar );
fprintf( ' \n\t Number of consumer characteristics:\t%2i\n', n.conChar );
fprintf( ' \n\t Number of estimated parameters:\t\t%2i\n', n.theta );
fprintf( ' \n ' );
fprintf( ' \n\t   Log-Likelihood: \t\t\t %10.4f ', MLE.value );
fprintf( ' \n\t  Run-time (sec.): \t\t\t %10.4f ', MLE.time );
fprintf( ' \n\t\t\t\t  AIC: \t\t\t %10.4f ', MLE.AIC );
fprintf( ' \n\t\t\t\t  BIC: \t\t\t %10.4f ', MLE.BIC );


fprintf( '\n\n ******************************************************\n' );
fprintf( ' ******************************************************\n' );

%% Print Estimated Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

headerIndex = 5;
thetaIndex  = 0;

fprintf( '\n\n\n   ********** Estimated Parameters ********** \n\n' );
    fprintf( '\t|\t\t\t\t\t\t|\t   Coef. |  Std. Err. | t-stat ' ) ;
    fprintf( '| p > |t| |   95 %% Conf. Interval |\n\n' );    
    

if ( 1 - spec.unobs ) > 0    
    
    fprintf( '\n\t ***** alpha_0 ***** \n\n\n' );
    headerIndex     = headerIndex + 1;
    thetaIndex      = thetaIndex  + 1;
    fprintf( ...
        '\t| %18s \t| %10.4f | %10.4f | %6.2f | %7.3f | %10.4f %10.4f |\n', ...
         dataHeader{ headerIndex }, thetaHat(thetaIndex), ...
         MLE.se(thetaIndex), MLE.tStat(thetaIndex), MLE.p(thetaIndex), ...
         MLE.ci_lb(thetaIndex), MLE.ci_ub(thetaIndex) );    
     
end

if n.conGroup > 0
    
    fprintf( '\n\n\t ***** alpha_r ***** \n\n\n' );   
    
    for r = 1 : n.conGroup
        headerIndex     = headerIndex + 1;
        thetaIndex      = thetaIndex  + 1;
    fprintf( ...
        '\t| %18s \t| %10.4f | %10.4f | %6.2f | %7.3f | %10.4f %10.4f |\n', ...
         dataHeader{ headerIndex }, thetaHat(thetaIndex), ...
         MLE.se(thetaIndex), MLE.tStat(thetaIndex), MLE.p(thetaIndex), ...
         MLE.ci_lb(thetaIndex), MLE.ci_ub(thetaIndex) );    
    end
    
end

if ( 1 - spec.unobs ) * n.prodChar > 0 
    fprintf( '\n\n\t ***** beta_1 ***** \n\n' );
    
    for j = 1 : n.maxChoice
        thetaIndex      = thetaIndex  + 1;
        fprintf( '\n\t| \t ** product %2i ** \t|\t\t', j );        
        fprintf( '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t |\n\n' );

        for s = 1 : n.prodChar
            headerIndex     = headerIndex + 1;
            thetaIndex      = thetaIndex + ( s - 1 ) * n.maxChoice;
            
            fprintf( '\t| %18s \t| %10.4f | %10.4f | %6.2f ', ...
                     dataHeader{ headerIndex }, ...
                     thetaHat(thetaIndex), MLE.se(thetaIndex), ...
                     MLE.tStat(thetaIndex) );
            fprintf( '| %7.3f | %10.4f %10.4f |\n', ...
                     MLE.p(thetaIndex), MLE.ci_lb(thetaIndex), ...
                     MLE.ci_ub(thetaIndex) );    

            thetaIndex      = thetaIndex - ( s - 1 ) * n.maxChoice;     
        end  
         
        headerIndex    = headerIndex - n.prodChar;
    end
    
    headerIndex     = headerIndex + n.prodChar;
    thetaIndex      = thetaIndex + ( n.prodChar - 1 ) * n.maxChoice;
end

if n.conChar > 0    
    fprintf( '\n\n\t ***** beta_2 ***** \n\n' );
    
    for j = 1 : n.maxChoice
        if j == spec.base
            fprintf( '\n\t| \t ** product %2i ** \t|\t\t', j );
            fprintf( '\t\t\t\t *** base alternative *** \t\t\t\t\t |\n\n' );
            
            for s = 1 : n.conChar
                headerIndex     = headerIndex + 1;
                fprintf( '\t| %18s \t| %10.4f | %10.4f | %6.2f ', ...
                        dataHeader{ headerIndex }, 0, 0, 0 );    
                fprintf( '| %7.3f | %10.4f %10.4f |\n', 0, 0, 0 );             
            end
        else
            thetaIndex      = thetaIndex  + 1;
            fprintf( '\n\t| \t ** product %2i ** \t|\t\t', j );        
            fprintf( '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t |\n\n' );
         
            for s = 1 : n.conChar
                headerIndex     = headerIndex + 1;
                thetaIndex      = thetaIndex + ...
                                  ( s - 1 ) * ( n.maxChoice - 1 );
                
                fprintf( '\t| %18s \t| %10.4f | %10.4f | %6.2f ', ...
                         dataHeader{ headerIndex }, ...
                         thetaHat(thetaIndex), MLE.se(thetaIndex), ...
                         MLE.tStat(thetaIndex) );
                fprintf( '| %7.3f | %10.4f %10.4f |\n', ...
                         MLE.p(thetaIndex), MLE.ci_lb(thetaIndex), ...
                         MLE.ci_ub(thetaIndex) );    
                         
                thetaIndex      = thetaIndex - ...
                                  ( s - 1 ) * ( n.maxChoice - 1 );     
            end
        end        
         
        headerIndex    = headerIndex - n.conChar;
    end
    
    headerIndex     = headerIndex + n.conChar;
    thetaIndex      = thetaIndex + ( n.conChar - 1 ) * ( n.maxChoice - 1 );
end

fprintf( '\n\n\t ****** Choleski Factor of The Differenced Covariance' );
fprintf( ' Matrix ( base alternative = %1i ) ****** \n\n\n', spec.base );

thetaIndex 	= thetaIndex + n.delta;
for i = 1 : n.maxChoice - 1
    for j = 1 : n.maxChoice - i
        if ( i == 1 ) && ( j == 1 )
            fprintf( '\t|\t\t\t\t s_%1i%1i \t| %10.4f | %10.4f ', ...
                     i, j, 0, 0 );    
            fprintf( '| %6.2f | %7.3f | %10.4f %10.4f |\n', 0, 0, 0, 0 );     
        else            
            thetaIndex 	= thetaIndex + 1;
            fprintf( '\t|\t\t\t\t s_%1i%1i \t| %10.4f | %10.4f ', ...
                     i, j, thetaHat(thetaIndex), MLE.se(thetaIndex) );    
            fprintf( '| %6.2f | %7.3f | %10.4f %10.4f |\n', ...
                     MLE.tStat(thetaIndex), MLE.p(thetaIndex), ...
                     MLE.ci_lb(thetaIndex), MLE.ci_ub(thetaIndex) );   
        end
    end
end

%% Print Covariance Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf( '\n\n   ********** Covariance Matrices ********** \n\n' );

fprintf( '\n\t ***** Full Covariance Matrices ***** \n\n' );
for j = 1 : n.maxChoice
    fprintf( '\t' );
    for k = 1 : n.maxChoice
        fprintf( '%12.4f ', omega( j, k ) );
    end
    fprintf( '\n' );
end

% Find a data group with all the alternatives; this is to obtain the
% correct M matrix used to get Differenced Covariance Matrices below
% !assumption: there exists such group (with all the alternatives)
for k = 1:numel(dataR)
    if dataR{k}.n.maxChoice == n.maxChoice
        kmax = k;
    end
end

for h = 1 : n.maxChoice  
        
    fprintf( '\n\t ***** Differenced Covariance Matrices');
    fprintf( ' ( base alternative = %1i ) ***** \n\n', h );

    tempOmega   = dataR{kmax}.M( :, :, h ,spec.base ) * S;
    tempOmega   = tempOmega * tempOmega';

    for j = 1 : ( n.maxChoice - 1 )
        fprintf( '\t' );
        for k = 1 : ( n.maxChoice - 1 )
            fprintf( '%12.4f ', tempOmega( j, k ) );
        end
        fprintf( '\n' );
    end
    
    clear tempOmega
end

%%
delta_table = reshape([delta(:) se_delta(:)]', [2*n.maxChoice n.market])';
fprintf( '\n***** Fixed Effect estimates');
printmat(delta_table, 'Fixed Effects', num2str(allmarkets'), 'choice_1  se_1 choice_2 se_2 choice_3 se_3 choice_4 se_4');


%% Run Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bounds on all parameters
ub  = ones( n.theta, 1 ) * 200;    
lb  = ones( n.theta, 1 ) * -200;

% Define objective function
obj = @(x) LogLike( x, dataR, spec );

% Run Optimization
tic
if spec.solver == 1
        
    optOption   = optimset( 'MaxIter', opt.maxIter, ...
                            'MaxFunEvals', opt.maxFunEvals, ...
                            'Display', opt.display, ...
                            'TolFun', opt.tolFun, ...
                            'TolCon', opt.tolCon, ...
                            'TolX', opt.tolX, ...
                            'Algorithm', opt.algorithm, ...
                            'GradObj', opt.gradObj, ...
                            'GradConstr', opt.gradConstr, ...
                            'DerivativeCheck', 'on' );
                        
    if spec.constraint == 1
        
        [ thetaHat, MLE.value ] = ...
                ktrlink( obj, theta_0, [], [], [], [], lb, ub, ...
                         @(x) NonLCon( x, n, spec.boundSize ), optOption );
                        
    elseif spec.constraint == 2
        
        [ thetaHat, MLE.value ] = ...
                ktrlink( obj, theta_0, [], [], [], [], lb, ub, ...
                         [], optOption );
    end
         

%     testOption  = optimset( 'MaxIter', 1, 'GradObj', 'off', 'Display','iter');    
%     [thetaHat,a,b,c,grad_fd]   = fminunc( @(x) LogLike( x, dataR, n, spec ), theta_0, testOption );
%     [f, grad] = LogLike( thetaHat, dataR, n, spec );
%     [grad_fd grad]
%     max(abs(grad_fd(2:end)-grad(2:end)))
%         
%     thetaHat    = fminunc( @(x) LogLike( x, dataR, n, spec ), theta_0, optOption );

elseif spec.solver == 2
    
    optOption   = optimset( 'MaxIter', opt.maxIter, ...
                            'MaxFunEvals', opt.maxFunEvals, ...
                            'Display', opt.display, ...
                            'TolFun', opt.tolFun, ...
                            'TolX', opt.tolX );
    
    [ thetaHat, MLE.value ] = fminsearch( obj, theta_0, optOption ); 
    
elseif spec.solver == 3
        
    optOption   = optimset( 'MaxIter', opt.maxIter, ...
                            'MaxFunEvals', opt.maxFunEvals, ...
                            'Display', opt.display, ...
                            'TolFun', opt.tolFun, ...
                            'TolCon', opt.tolCon, ...
                            'TolX', opt.tolX, ...
                            'Algorithm', opt.algorithm, ...
                            'GradObj', opt.gradObj, ...
                            'GradConstr', opt.gradConstr, ...
                            'DerivativeCheck', 'off' );
                        
    if spec.constraint == 1
        
        [ thetaHat, MLE.value ] = ...
                fmincon( obj, theta_0, [], [], [], [], lb, ub, ...
                         @(x) NonLCon( x, n, spec.boundSize ), optOption );
                        
    elseif spec.constraint == 2
        
        [ thetaHat, MLE.value ] = ...
                fmincon( obj, theta_0, [], [], [], [], lb, ub, ...
                         [], optOption );
    end
end

% Estimation time in seconds
MLE.time    = toc;

%% Compute MLE Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define negative log likelihood function
nLoglike    = @( x, data, cens, freq ) obj( x );

% Compute MLE covariance matrix by finite difference
fprintf('\n Computing MLE Covariance Matrix...');
MLE.covOpt  = statset( 'GradObj', opt.gradObj );
MLE.cov     = mlecov( thetaHat, data, 'nloglf', nLoglike, ...
                      'options', MLE.covOpt );

% Standard errors
MLE.se      = sqrt( diag( MLE.cov ) );

% t-statistics
MLE.tStat   = thetaHat ./ MLE.se;

% p > |t|
MLE.p       = 2 * ( 1 - normcdf( abs( MLE.tStat ) ) );

% 95% Confidence Interval
MLE.ci_lb   = thetaHat - 1.96 * MLE.se;
MLE.ci_ub   = thetaHat + 1.96 * MLE.se;

% MLE value
MLE.value   = -MLE.value;
                  
% AIC
MLE.AIC     = 2 * ( n.theta - MLE.value );

% BIC
MLE.BIC     = -2 * MLE.value + n.theta * log( n.con );

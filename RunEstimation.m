
%% Run Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bounds on all parameters
ub  = ones( n.theta, 1 ) * 2;    
lb  = ones( n.theta, 1 ) * -2;
ub(end-1:end) = 20;
lb(end-1:end) = -20;
ub(1) = 50;
lb(1) = -50;

optOption   = optimset( 'MaxIter', opt.maxIter, ...
    'MaxFunEvals', opt.maxFunEvals, ...
    'GradObj', 'on', ... opt.gradObj, ...
    'GradConstr', 'on', ... opt.gradConstr, ...
    'Display', opt.display, ...
    'Hessian', 'user-supplied', ...
    'HessFcn', @(x,y) SSE_hss(x, dataS, marketIdByCon, mapConID, n, shareHat, mask), ... ); %, ...
    'TolFun', 1e-5, ...
    'DerivativeCheck', 'on' );
%         'TolFun', opt.tolFun, ...
%         'TolCon', opt.tolCon, ...
%         'TolX', opt.tolX, ...
%         'Algorithm', opt.algorithm, ...
%         'DerivativeCheck', 'on' );


% Find a feasible solution by minizing sse
% [ theta_0, sse ] = ...
%     ktrlink( @(x)SSE(x, dataS, marketIdByCon, mapConID, n, shareHat, mask ), ...
%     theta_0, [], [], [], [], lb, ub, ...
%     [], optOption );



% Bounds on all parameters
ub  = ones( n.theta, 1 ) * 100;    
lb  = ones( n.theta, 1 ) * -100;

% Define objective function
obj = @(x) LogLike( x, dataR );
%obj = @(x) NFXLogLike(x, dataR, dataS, marketIdByCon, n, shareHat, mask, opt);
constraints = @(x) ShareConstraints(x,dataS,marketIdByCon,mapConID,n,shareHat,mask);

index = 1+n.conGroup+(n.maxChoice-1)*(n.conChar + n.prodChar);

% Run Optimization
tic
if spec.solver == 1
    
    optOption   = optimset( 'MaxIter', opt.maxIter, ...
        'MaxFunEvals', opt.maxFunEvals, ...
        'GradObj', opt.gradObj, ...
        'GradConstr', opt.gradConstr, ...
        'Display', opt.display,... ); %, ...
        'TolFun', opt.tolFun, ...
        'TolCon', opt.tolCon, ...
        'TolX', opt.tolX); %, ...
%         'Algorithm', opt.algorithm, ...
%         'DerivativeCheck', 'on' );
    
    optOption2   = optimset( 'MaxIter', 20000, ...
        'MaxFunEvals', opt.maxFunEvals, ...
        'Display', opt.display, ...
        'TolFun', opt.tolFun, ...
        'TolCon', opt.tolCon, ...
        'TolX', opt.tolX, ...
        'Algorithm', opt.algorithm, ...
        'SubproblemAlgorithm', 'cg', ...
        'GradObj', opt.gradObj, ...
        'GradConstr', opt.gradConstr );
                       
    if spec.constraint == 1
        
        [ thetaHat, MLE.value ] = ...
                knitromatlab( obj, theta_0, [], [], [], [], lb, ub, ...
                         @(x) NonLCon( x, n, spec.boundSize ), optOption );
                     
    elseif spec.constraint == 2
        
        [ thetaHat, MLE.value ] = ...
                knitromatlab( obj, theta_0, [], [], [], [], lb, ub, ...
                         [], [], optOption2 );
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
elseif spec.solver == 4
    funcs.objective = obj;
    funcs.gradient = @(x) LogLikeGrad(x, dataR);
    
    options.lb = lb;
    options.ub = ub;
    
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.derivative_test = 'first-order';
    options.ipopt.max_iter = opt.maxIter;
%     options.ipopt.derivative_test_perturbation = 1e-5;
    
    [ thetaHat, info ] = ipopt(theta_0, funcs, options);
    MLE.value = obj(thetaHat);
end

% Estimation time in seconds
MLE.time    = toc;

%% Compute MLE Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define negative log likelihood function
nLoglike    = @( x, data, cens, freq ) obj( x );

% Compute MLE covariance matrix by finite difference
fprintf('\n Computing MLE Covariance Matrix...');
MLE.covOpt  = statset( 'GradObj', opt.gradObj );
%MLE.cov     = mlecov( thetaHat, data, 'nloglf', nLoglike, ...
%                      'options', MLE.covOpt );
MLE.cov = bhhh(thetaHat, dataR);

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
